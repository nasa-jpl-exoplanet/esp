"""
-----------------
Generates the TEA pre-atmospheric structure as a python object rather than a file

Uses the same abundance-conversion logic used in TEA from Asplund-2009 (https://arxiv.org/pdf/0909.0948 (Table 1, pg 42.)) solar data.
"""
import numpy as np
from excalibur.util.tea_code import readconf

__all__ = ["build_pre_atm"]


def _dex_to_fraction(dex_values):
    """Convert log10(dex) to number-fraction relative to hydrogen."""
    num_dens = 10.0 ** np.asarray(dex_values, dtype=float)
    H_num = 10.0**12.0  # Using Asplund table
    return (num_dens / H_num).tolist()


def build_pre_atm(
    pressure,
    temperature,
    abundances_path=None,
    input_elem=None,
    output_species=None,
    cfg_file="TEA.cfg",
):
    """
    Parameters
    ----------
    pressure, temperature : 1-D arrays (same length)
    abundances_path       : path to Asplund-style abundances.txt
    input_elem            : iterable of element symbols   (H C N …)
    output_species        : iterable of species strings   (H2O_g …)
    cfg_file              : configuration file

    Returns
    -------
    dict  {
        "pressure"        : np.ndarray [n_layers],
        "temperature"     : np.ndarray [n_layers],
        "atom_name"       : list[str]     (order kept from input_elem),
        "atom_abundances" : np.ndarray [n_layers, n_elem]  (fractions),
        "output_species"  : list[str]
    }
    """

    TEApars, PREATpars = readconf.readcfg(cfg_file)
    _, _, _, _, abun_file_cfg, _, _, _ = TEApars
    _, _, input_elem_cfg, output_spec_cfg = PREATpars

    input_elem = input_elem or input_elem_cfg.split()
    output_species = output_species or output_spec_cfg.split()
    abundances_path = abundances_path or abun_file_cfg

    data = np.genfromtxt(abundances_path, comments='#', dtype=str)

    symbols = data[:, 1]
    dex = data[:, 2].astype(float)
    sel = np.isin(symbols, input_elem)
    frac = _dex_to_fraction(dex[sel])

    atom_name = symbols[sel].tolist()
    n_layers = len(pressure)
    atom_abund = np.tile(frac, (n_layers, 1))

    return {
        "pressure": np.asarray(pressure, dtype=float),
        "temperature": np.asarray(temperature, dtype=float),
        "atom_name": atom_name,
        "atom_abundances": atom_abund.astype(float),
        "output_species": output_species,
    }


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    raise RuntimeError("Import and call build_pre_atm(...) from your driver.")
