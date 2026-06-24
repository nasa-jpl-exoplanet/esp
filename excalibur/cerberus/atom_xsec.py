'''Computation of the grids for atomic cross sections Ca, K, Na
with VALD line lists'''

import numpy as np
import json
import warnings

from scipy.special import wofz
from scipy.interpolate import interp1d
from scipy.constants import c, h, k, m_e, physical_constants, eV

# reference physical parameters for
# colisionnal broadening
TREF = 296e0
PREF = 1e0

# physical constants
ecgs = 4.80320427e-10  # Electron charge in cgs units (statC)
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
        "/proj/sdp/data/CERBERUS/ATOM_XSEC/lines_infos/lines.json", "r"
    ) as f:
        data = json.load(f)
    # reconvertir les clés wavelength en float
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
# -- COLISIONNAL BROADENING COEFFICIENTS COMPUTATION -- --------------
def gammavld(gamma_vdw, m_s, broadener):

    alphaH = 0.666793  # Polarisability of atomic hydrogen (A^-3)
    mH = 1.007825  # Mass of atomic hydrogen (u)

    if broadener == 'H2':
        alphap = 0.805000  # Polarisability of molecular hydrogen (A^-3)
        mp = 2.01565  # Mass of molecular hydrogen (u)

    elif broadener == 'He':
        alphap = 0.205052  # Polarisability of helium (A^-3)
        mp = 4.002603  # Mass of helium (u)

    # Compute Lorentzian HWHM
    gammaL0 = (
        2.2593427e7
        * gamma_vdw
        * np.power(((mH * (m_s + mp)) / (mp * (m_s + mH))), (3e0 / 1e1))
        * np.power((alphap / alphaH), (2e0 / 5e0))
    )

    # Temperature exponent
    nL = 7e-1

    return gammaL0, nL


# --------- ----------------------------------------------------------
# -- COLISIONNAL BROADENING COMPUTATION -- ---------------------------
def H2Hebroadening(
    gamma0H2,
    T,
    nLH2,
    pressure,
    XH2,
    gamma0He,
    nLHe,
):
    gamma = gamma0H2 * np.power((TREF / T), nLH2) * (
        pressure / PREF
    ) * XH2 + gamma0He * np.power(  # H2+He Lorentzian HWHM for given T, pressure, and J (ang. mom.)
        (TREF / T), nLHe
    ) * (
        pressure / PREF
    ) * (
        1e0 - XH2
    )

    return gamma


# --------- ----------------------------------------------------------
# -- PARTITION FUNCTION COMPUTATION -- -------------------------------
def partition_function(specie, T):
    temp, Q = load_part_func(specie)
    Tmin, Tmax = np.min(temp), np.max(temp)
    f = interp1d(temp, Q, kind='linear', fill_value="extrapolate")

    def Qwarn(T):
        if np.any((T < Tmin) | (T > Tmax)):
            warnings.warn(
                f"Extrapolation used for T = {T}. "
                f"Valid range is [{Tmin}, {Tmax}]",
                RuntimeWarning,
            )
        return f(T)

    return Qwarn(T)


# --------- ----------------------------------------------------------
# -- LINE INTENSITY COMPUTATION -- -----------------------------------
def Svald(gf, Elow, nu0, T, Q):
    S = (
        ((gf * np.pi * ecgs**2) / (me_cgs * c_cgs**2))
        * (1e0 / Q)
        * np.exp(-1e0 * c2 * Elow / T)
        * (1e0 - np.exp(-1e0 * c2 * nu0 / T))
    )
    return S


# --------- ----------------------------------------------------------
# -- VOIGT PROFILE COMPUTATION -- ------------------------------------
def Voigt(nu_grid, nu0, T, mass, gamma):
    alpha = np.sqrt(2e0 * k * T * np.log(2e0) / mass) * nu0 / c
    x = np.sqrt(np.log(2e0)) * (np.outer(1.0 / alpha, (nu_grid - nu0)))
    y = np.sqrt(np.log(2e0)) * (np.outer(gamma / alpha, np.ones(len(nu_grid))))
    z = x + 1j * y
    coeff = np.sqrt(np.log(2e0) / np.pi) * np.outer(
        1e0 / alpha, np.ones(len(nu_grid))
    )
    return coeff * np.real(wofz(z))


# --------- ----------------------------------------------------------
# -- CROSS SECTION COMPUTATION FOR 1 LINE-- --------------------------
def single_line_sigma(w_grid, specie, line, Q, T, pressure, XH2):

    # line infomations
    nu0 = 1.0 / line["WL_vacuum"] * 1e4  # conversion from micron to cm-1
    nu_grid = 1.0 / w_grid * 1e4  # conversion from micron to cm-1
    gf = line["gf"]
    Elow = (
        line["E_low"] * eV / h / c_cgs
    )  # conversion from eV to cm-1 for the computation of the line intensity
    gamma_nat = line["Rad"]
    gamma_vdw = line["Waals"]

    # line intensity computation
    S = Svald(gf, Elow, nu0, T, Q)

    # broadening
    gamma0H2, nLH2 = gammavld(gamma_vdw, ATOMIC_MASS[specie], 'H2')
    gamma0He, nLHe = gammavld(gamma_vdw, ATOMIC_MASS[specie], 'He')
    gamma_brd = H2Hebroadening(
        gamma0H2, T, nLH2, pressure, XH2, gamma_0_He, nLHe
    )
    gamma = gamma_brd + 1e0 / (4e0 * np.pi * (c_cgs)) * gamma_nat

    # line cutoff
    wing_cutoff = 3e1

    # voigt profile construction
    idx = np.where(np.abs(nu_grid - nu0) <= wing_cutoff)[0]
    sigma = np.zeros((len(T), len(nu_grid)))
    if len(idx):
        phi = Voigt(nu_grid[idx], nu0, T, mass, gamma)
        sigma[:, idx] = S[:, None] * phi

    return sigma


# --------- ----------------------------------------------------------
# -- CROSS SECTION COMPUTATION WITH LINE BY LINE -- ------------------
def atom_xsec(
    w_grid,
    specie,
    T,
    pressure,
    XH2,
):
    # partition function computation
    Q = partition_function(specie, T)

    # lines loading
    data = load_lines()

    sigma_total = np.zeros((len(T), len(wavelength_grid)))
    for line in data[specie].values():
        sigma_total += single_line_sigma(
            w_grid, specie, line, Q, T, pressure, XH2
        )
    return sigma_total


# --------- ----------------------------------------------------------
# -- GRID STORING IN /proj/sdp/data/CERBERUS/ATOM_XSEC -- ------------
def grid_generation(elem, temp, pressure, XH2, wgrid):

    Tgrid, Pgrid, XH2grid = np.meshgrid(temp, pressure, XH2, indexing='ij')

    Tcomb = T_grid.flatten()
    Pcomb = P_grid.flatten()
    XH2comb = XH2grid.flatten()

    grid = atom_xsec(
        w_grid=wgrid,
        specie=elem,
        T=Tcomb,
        pressure=Pcomb,
        XH2=XH2comb,
    )

    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/wgrid.npy", wgrid)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/temp.npy", temp)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/pressure.npy", pressure)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/X_H2.npy", pressure)

    grid_4d = grid.reshape((len(temp), len(pressure), len(XH2), -1))

    np.save(
        "/proj/sdp/data/CERBERUS/ATOM_XSEC/" + elem + "/grid_4d.npy", grid_4d
    )

    return
