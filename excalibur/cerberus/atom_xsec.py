'''Computation of the grids for atomic cross sections Ca, K, Na
with VALD line lists ds'''

# pylint: disable=invalid-name
# pylint: disable=no-member

import numpy as np
import json
import warnings

from scipy import special
from scipy.interpolate import interp1d
from scipy.constants import c, h, k, m_e, physical_constants, eV

# reference physical parameters for
# colisional broadening
TREF = 296e0
PREF = 1e0

# physical constants
e_cgs = 4.80320427e-10  # Electron charge in cgs units (statC)
c_cgs = 1e2 * c  # Speed of light in cgs units (cm)
me_cgs = 1e3 * m_e  # Electron mass in cgs units (g)
c2 = h * c_cgs / k
amu = physical_constants["atomic mass constant"][0]

# Element masses
ATOMIC_MASS = {
    "H": 1.008,
    "He": 4.0026,
    "Na": 22.990,
    "K": 39.098,
    "Ca": 40.078,
}


# --------- ----------------------------------------------------------
# -- LINE LIST LOADING -- --------------------------------------------
def load_lines():
    with open(
        "/proj/sdp/data/CERBERUS/ATOM_XSEC/lines_infos/lines.json",
        "r",
        encoding="utf-8",
    ) as f:
        data = json.load(f)
    for element in data:
        if element == "units":
            continue
        data[element] = {
            float(wl): values for wl, values in data[element].items()
        }
    return data


# --------- ----------------------------------------------------------
# -- PARTITION FUNCTION LOADING -- -----------------------------------
def load_part_func(specie):
    with open(
        "/proj/sdp/data/CERBERUS/ATOM_XSEC/lines_infos/partition_function.json",
        "r",
        encoding="utf-8",
    ) as f:
        partition_func = json.load(f)
        for element in partition_func:
            partition_func[element]["T"] = np.array(
                [float(x) for x in partition_func[element]["T"]]
            )
            partition_func[element]["Q"] = np.array(
                [float(x) for x in partition_func[element]["Q"]]
            )
    temp = partition_func[specie]["T"]
    Q = partition_func[specie]["Q"]
    return temp, Q


# --------- ----------------------------------------------------------
# -- COLISIONAL BROADENING COEFFICIENTS -- ---------------------------
def gammavld(gamma_vdw, ms, broadener):

    alphah = 0.666793  # Polarisability of atomic hydrogen (A^-3)
    mh = 1.007825  # Mass of atomic hydrogen (u)
    mp = 0.0
    alphap = 0.0

    if broadener == 'H2':
        alphap = 0.805000  # Polarisability of molecular hydrogen (A^-3)
        mp = 2.01565  # Mass of molecular hydrogen (u)

    elif broadener == 'He':
        alphap = 0.205052  # Polarisability of helium (A^-3)
        mp = 4.002603  # Mass of helium (u)

    # Compute Lorentzian HWHM
    gammal0 = (
        2.2593427e7
        * gamma_vdw
        * np.power(((mh * (ms + mp)) / (mp * (ms + mh))), (3e0 / 1e1))
        * np.power((alphap / alphah), (2e0 / 5e0))
    )

    return gammal0


# --------- ----------------------------------------------------------
# -- COLISIONAL BROADENING -- ----------------------------------------
def h2hebroadening(
    gamma0h2,
    T,
    pressure,
    xh2,
    gamma0he,
):
    gamma = gamma0h2 * np.power((TREF / T), 7e-1) * (
        pressure / PREF
    ) * xh2 + gamma0he * np.power(  # H2+He Lorentzian HWHM for given T, pressure, and J (ang. mom.)
        (TREF / T), 7e-1
    ) * (
        pressure / PREF
    ) * (
        1e0 - xh2
    )

    return gamma


# --------- ----------------------------------------------------------
# -- PARTITION FUNCTION -- -------------------------------------------
def partition_function(specie, T):
    temp, Q = load_part_func(specie)
    temp_min, temp_max = np.min(temp), np.max(temp)
    f = interp1d(temp, Q, kind='linear', fill_value="extrapolate")

    if np.any((T < temp_min) | (T > temp_max)):
        warnings.warn(
            f"Extrapolation used for T = {T}. "
            f"Valid range is [{temp_min}, {temp_max}]",
            RuntimeWarning,
        )

    return f(T)


# --------- ----------------------------------------------------------
# -- LINE INTENSITY -- -----------------------------------------------
def vald_intensity(gf, elow, nu0, T, Q):
    S = (
        ((gf * np.pi * e_cgs**2) / (me_cgs * c_cgs**2))
        * (1e0 / Q)
        * np.exp(-1e0 * c2 * elow / T)
    for line in data[specie].values():
        sigma_total += single_line_sigma(w_grid, specie, line, Q, parameters)
    return sigma_total


# --------- ----------------------------------------------------------
# -- GRID STORING IN /proj/sdp/data/CERBERUS/ATOM_XSEC -- ------------
def grid_generation(elem, temp, pressure, xh2, wgrid):

    temp_grid, pressure_grid, xh2_grid = np.meshgrid(
        temp, pressure, xh2, indexing='ij'
    )

    temp_comb = temp_grid.flatten()
    pressure_comb = pressure_grid.flatten()
    xh2_comb = xh2_grid.flatten()

    parameters = np.column_stack((temp_comb, pressure_comb, xh2_comb))

    grid = atom_xsec(w_grid=wgrid, specie=elem, parameters=parameters)

    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/temp.npy", temp)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/pressure.npy", pressure)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/X_H2.npy", xh2)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/wgrid.npy", wgrid)

    grid_4d = grid.reshape((len(temp), len(pressure), len(xh2), -1))

    np.save(
        "/proj/sdp/data/CERBERUS/ATOM_XSEC/" + elem + "/grid_4d.npy", grid_4d
    )

    return
