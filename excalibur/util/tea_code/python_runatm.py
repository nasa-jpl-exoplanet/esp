"""
python_runatm.py
----------------
Thermochemical-equilibrium driver that consumes the in-memory pre-atmosphere
dictionary produced by `python_makeatm.build_pre_atm` and returns a pandas
DataFrame instead of a `.tea` file.

Mathematical core (iterate, makeheader, etc.) is unchanged from original TEA framework.
Only:
  • removed all disk reads of *.atm / abundances.txt;
  • removed final file write of results;
  • refactored into a callable `run_tea(...)` function that returns a DataFrame.
"""
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import ctypes
from pathlib import Path

from excalibur.util.tea_code import readconf
from excalibur.util.tea_code import iterate
from excalibur.util.tea_code import makeheader
from excalibur.util.tea_code import updated_balance


# -----------------------------------------------------------------------------


def _multiproc_worker(
    pres_arr,
    temp_arr,
    atom_arr,
    free_energy,
    heat,
    stoich_arr,
    guess,
    maxiter,
    verb,
    times,
    xtol,
    shared_abn,
):
    """Exactly the same per-layer loop as original runatm, minus file I/O."""
    abn = np.ctypeslib.as_array(shared_abn.get_obj()).reshape(
        (len(pres_arr), len(stoich_arr))
    )
    n_pres = pres_arr.size
    for q in range(n_pres):
        if verb > 1:
            print(f"\nLayer {q + 1:d}:")
        g_RT = makeheader.calc_gRT(free_energy, heat, temp_arr[q])
        try:

            max_retry = 5
            for n in range(max_retry):
                result = iterate.iterate(
                    pres_arr[q],
                    stoich_arr,
                    atom_arr[q],
                    g_RT,
                    maxiter,
                    verb,
                    times,
                    guess,
                    xtol,
                )
                if isinstance(result, tuple):
                    break
                xtol *= 0.3

            else:
                continue

            if not isinstance(result, tuple) or len(result) < 6:
                raise ValueError("iterate returned scalar")

            y, x, delta, y_bar, x_bar, delta_bar = result
            guess = (x, x_bar)
            abn[q, :] = x / x_bar

        except Exception as exc:
            print(
                f"[Layer {q + 1}] iterate failed: {exc} — using balanced guess"
            )
            x, x_bar = guess
            abn[q, :] = x / x_bar
        else:
            y, x, delta, y_bar, x_bar, delta_bar = result
        continue


def run_tea(pre_atm, cfg_file, desc="tea_output"):
    """
    Parameters
    ----------
    pre_atm : dict      – output from `build_pre_atm`
    cfg_file: str       – TEA.cfg path (defaults to current dir)
    desc    : str       – label used only for log/print messages

    Returns
    -------
    pandas.DataFrame  with columns ['Pressure', 'Temp', *species ]
    """
    TEApars, _ = readconf.readcfg(cfg_file)
    maxiter, savefiles, verb, times, _, location_out, xtol, _ = TEApars
    thermo_dir = os.path.join(os.path.dirname(cfg_file), "gdata")

    pres_arr = np.asarray(pre_atm["pressure"], dtype=float)
    temp_arr = np.asarray(pre_atm["temperature"], dtype=float)
    atom_arr = np.asarray(pre_atm["atom_abundances"], dtype=float)
    atom_name = np.asarray(pre_atm["atom_name"])
    speclist = np.asarray(pre_atm["output_species"])

    n_pres = pres_arr.size
    n_spec = speclist.size

    free_energy, heat = makeheader.read_gdata(speclist, thermo_dir)

    thermo_dir = Path(cfg_file).with_name("gdata")

    missing = [
        sp
        for sp in speclist
        if not any(
            (thermo_dir / f"{sp}{ext}").exists() for ext in (".dat", ".txt")
        )
    ]
    if missing:
        raise FileNotFoundError(
            f"Thermo file not found for species: {', '.join(missing)} "
            f"in {thermo_dir}. Add the file(s) or remove those species."
        )

    stoich_arr, elements = makeheader.read_stoich(
        speclist, stoich_file=Path(cfg_file).with_name("stoich.txt")
    )
    elem_idx = [np.where(atom_name == el)[0][0] for el in elements]
    atom_arr = atom_arr[:, elem_idx]  # (n_layers, n_elements)

    guess = updated_balance.balance(stoich_arr, atom_arr[0], verb)

    shared_abn = mp.Array(ctypes.c_double, int(n_pres * n_spec))
    p = mp.Process(
        target=_multiproc_worker,
        args=(
            pres_arr,
            temp_arr,
            atom_arr,
            free_energy,
            heat,
            stoich_arr,
            guess,
            maxiter,
            verb,
            times,
            xtol,
            shared_abn,
        ),
    )
    p.start()
    p.join()  # mandatory

    abn = np.ctypeslib.as_array(shared_abn.get_obj()).reshape((n_pres, n_spec))

    cols = ["Pressure", "Temp"] + speclist.tolist()
    df = pd.DataFrame(np.column_stack((pres_arr, temp_arr, abn)), columns=cols)
    return df
