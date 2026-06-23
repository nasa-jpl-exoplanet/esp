''' Computation of the grids for atomic cross sections Ca, K, Na '''

import numpy as np
import json
import warnings

from scipy.special import wofz
from scipy.interpolate import interp1d
from scipy.constants import (c,h,k,m_e,e,pi,physical_constants,eV)
T_ref = 296.0
P_ref = 1.0
e_cgs = 4.80320427e-10  # Electron charge in cgs units (statC)
c_cgs = 100.0*c         # Speed of light in cgs units (cm)
m_e_cgs = 1000.0*m_e    # Electron mass in cgs units (g)
c2 = h*c_cgs/k

amu = physical_constants["atomic mass constant"][0]

ATOMIC_MASS = {
    "H": 1.008,
    "He": 4.0026,
    "Na": 22.990,
    "K": 39.098,
    "Ca": 40.078
}

def get_mass(specie):
    elem = specie.split()[0]
    return ATOMIC_MASS[elem]*amu

def load_lines():
    with open("/proj/sdp/data/CERBERUS/ATOM_XSEC/lines_infos/lines.json", "r") as f:
        data = json.load(f)
    # reconvertir les clés wavelength en float
    for element in data:
        if element == "units":
            continue
        data[element] = {
            float(wl): values
            for wl, values in data[element].items()
    }
    return data

def load_part_func(specie):
    with open("/proj/sdp/data/CERBERUS/ATOM_XSEC/lines_infos/partition_function.json", "r") as f:
        partition_func = json.load(f)
        for element in partition_func:
            partition_func[element]["T"] = np.array([float(x) for x in partition_func[element]["T"]])
            partition_func[element]["Q"] = np.array([float(x) for x in partition_func[element]["Q"]])
    T_raw = partition_func[specie]["T"]
    Q_T_raw = partition_func[specie]["Q"]
    return T_raw,Q_T_raw

def gamma_L_VALD(gamma_vdw, m_s, broadener):

    alpha_H = 0.666793  # Polarisability of atomic hydrogen (A^-3)
    m_H = 1.007825      # Mass of atomic hydrogen (u)

    if (broadener == 'H2'):
        alpha_p = 0.805000  # Polarisability of molecular hydrogen (A^-3)
        m_p = 2.01565       # Mass of molecular hydrogen (u)

    elif (broadener == 'He'):
        alpha_p = 0.205052  # Polarisability of helium (A^-3)
        m_p = 4.002603      # Mass of helium (u)

    # Compute Lorentzian HWHM
    gamma_L_0 = (2.2593427e7 * gamma_vdw * np.power(((m_H*(m_s+m_p))/(m_p*(m_s+m_H))), (3.0/10.0)) *
                                           np.power((alpha_p/alpha_H), (2.0/5.0)))

    # Temperature exponent
    n_L = 0.7

    return gamma_L_0, n_L

def read_atom(gamma_vdw, m):
    gamma_0_H2, n_L_H2 = gamma_L_VALD(gamma_vdw, (m/amu), 'H2')
    gamma_0_He, n_L_He = gamma_L_VALD(gamma_vdw, (m/amu), 'He')
    return gamma_0_H2, gamma_0_He, n_L_H2, n_L_He

def compute_H2_He_broadening(gamma_0_H2, T_ref, T, n_L_H2, P, P_ref, X_H2, gamma_0_He, n_L_He, X_He):
    gamma = (gamma_0_H2 * np.power((T_ref/T), n_L_H2) * (P/P_ref) * X_H2 +   # H2+He Lorentzian HWHM for given T, P, and J (ang. mom.)
             gamma_0_He * np.power((T_ref/T), n_L_He) * (P/P_ref) * X_He)    # Note that these are only a function of J''

    return gamma

def partition_function(specie, T):
    T_raw, Q_T_raw = load_part_func(specie)
    T_min, T_max = np.min(T_raw), np.max(T_raw)
    f = interp1d(T_raw, Q_T_raw, kind='linear', fill_value="extrapolate")
    def Q_with_warning(T):
        if np.any((T < T_min) | (T > T_max)):
            warnings.warn(
                f"Extrapolation used for T = {T}. "
                f"Valid range is [{T_min}, {T_max}]",
                RuntimeWarning
            )
        return f(T)
    
    return Q_with_warning(T)

    
def compute_line_intensity_VALD(gf, E_low, nu0, T, Q_T):
    S = ((gf*pi*e_cgs**2)/(m_e_cgs*c_cgs**2)) * (1.0/Q_T) * np.exp(-1.0*c2*E_low/T) * (1.0 - np.exp(-1.0*c2*nu0/T))
    return S

def Voigt_profile(nu_grid,nu0,T,mass,gamma):
    alpha = np.sqrt(2.0*k*T*np.log(2)/mass)*nu0/c
    x = np.sqrt(np.log(2.0)) * ((nu_grid-nu0)/alpha)
    y = np.sqrt(np.log(2.0)) * (gamma/alpha)
    z = x + 1j * y
    coeff = np.sqrt(np.log(2.0)/np.pi) * (1.0/alpha)
    return (coeff * np.real(wofz(z)))

def single_line_sigma(wavelength_grid,specie,line,Q_T,T,P,T_ref,P_ref,X_H2,X_He):
    
    #information for each line
    nu0 = 1./line["WL_vacuum"] * 1e4
    nu_grid = 1./wavelength_grid * 1e4
    gf = line["gf"]
    E_low = line["E_low"] * eV / h /c_cgs #conversion from eV to cm-1 for the computation of the line intensity
    gamma_nat = line["Rad"]
    gamma_vdw = line["Waals"]
    mass = get_mass(specie)

    # collisional broadening
    gamma_0_H2, gamma_0_He, n_L_H2, n_L_He = read_atom(gamma_vdw, mass)

    # line cutoff at 30 cm-1 away from the line 
    wing_cutoff = 30.

    sigma = np.zeros((len(T), len(wavelength_grid)))
    for i in range(len(T)):
        #line intensity computation
        S = compute_line_intensity_VALD(gf, E_low, nu0, T[i], Q_T[i])

        #broadening computation
        gamma_brd = compute_H2_He_broadening(gamma_0_H2, T_ref, T[i], n_L_H2, P[i], P_ref, X_H2, gamma_0_He, n_L_He, X_He)
        gamma = gamma_brd  + 1.0/(4.0*np.pi*(c_cgs)) * gamma_nat

        idx = np.where(np.abs(nu_grid - nu0) <= wing_cutoff)[0]
        # line profile computation
        if len(idx):
            phi = Voigt_profile(nu_grid[idx],nu0,T[i],mass,gamma)
            sigma[i,idx] = S * phi

    return sigma

def atom_xsec(wavelength_grid,specie,T,P,T_ref=T_ref,P_ref=P_ref,X_H2=0.85,X_He=0.15):
    
    #partition function computation
    Q_T = partition_function(specie, T)
    
    #lines loading
    data = load_lines()
    sigma_total = np.zeros((len(T), len(wavelength_grid)))

    #line by line computation
    for line in data[specie].values():
        sigma_total += single_line_sigma(wavelength_grid,specie,line,Q_T,T,P,T_ref,P_ref,X_H2,X_He)
    return sigma_total

def grid_generation(elem, temp, pressure, wgrid):

    T_grid, P_grid = np.meshgrid(temp, pressure, indexing = 'ij')

    T_comb = T_grid.flatten()
    P_comb = P_grid.flatten()

    grid = atom_xsec(wavelength_grid=wgrid,specie=elem,T=T_comb,P=P_comb)

    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/wgrid.npy",wgrid)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/temp.npy",temp)
    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/pressure.npy",pressure)

    grid_3d = grid.reshape((len(temp), len(pressure), -1))

    np.save("/proj/sdp/data/CERBERUS/ATOM_XSEC/" + elem + "/grid_3d.npy", grid_3d)

    return