'''utilites used by many tasks'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-locals,too-many-statements,too-many-arguments,too-many-positional-arguments

# -- IMPORTS -- ------------------------------------------------------
import numpy as np


# ------------ -------------------------------------------------------
# -- CHEMICAL EQUILIBRIUM -- -----------------------------------------
def TEA(p, temp, C2Or=0.0, X2Hr=0.0, N2Or=0.0):
    '''
    TEA chemical equilibrium abundances
    '''

    mixratio, nH2, nHe = crbce(
        p, temp, X2Hr=X2Hr * 1.1, C2Or=C2Or * 0.9, N2Or=N2Or
    )

    return mixratio, nH2, nHe


# ------------ -------------------------------------------------------
# -- CHEMICAL EQUILIBRIUM -- -----------------------------------------
def crbce(p, temp, C2Or=0.0, X2Hr=0.0, N2Or=0.0):
    '''
    G. ROUDIER: BURROWS AND SHARP 1998 + ANDERS & GREVESSE 1989
    '''
    solar = {
        'nH': 9.10e-1,
        'nHe': 8.87e-2,
        'nO': 7.76e-4,
        'nC': 3.29e-4,
        'nNE': 1.12e-4,
        'nN': 1.02e-4,
        'nMg': 3.49e-5,
        'nSi': 3.26e-5,
        'nFe': 2.94e-5,
        'nS': 1.68e-5,
        'nAr': 3.29e-6,
        'nAl': 2.77e-6,
        'nCa': 1.99e-6,
        'nNa': 1.87e-6,
        'nNi': 1.61e-6,
        'nCr': 4.40e-7,
        'nP': 3.39e-7,
        'nMn': 3.11e-7,
        'nCl': 1.71e-7,
        'nK': 1.23e-7,
        'nTi': 7.83e-8,
        'nCo': 7.34e-8,
        'nF': 2.75e-8,
        'nV': 9.56e-9,
        'nLi': 1.86e-9,
        'nRb': 2.31e-10,
        'nCs': 1.21e-11,
    }
    a1 = 1.106131e6
    b1 = -5.6895e4
    c1 = 62.565
    d1 = -5.81396e-4
    e1 = 2.346515e-8
    RcalpmolpK = 1.9872036  # cal/mol/K

    solCtO = solar['nC'] / solar['nO']
    solNtO = solar['nN'] / solar['nO']

    # define metal as all elements except H and He
    metal = solar.copy()
    metal['nH'] = 0.0
    metal['nHe'] = 0.0
    # sum up the metal part.  It's 10^-2.84 dex
    nXsolar = np.sum(np.array([x[1] for x in metal.items()]))

    # (C2O,N2O are limited by the prior bounds, so limits not really needed)
    C2Or = max(C2Or, -10.0)
    C2Or = min(C2Or, 10.0)
    N2Or = max(N2Or, -10.0)
    N2Or = min(N2Or, 10.0)

    Xfactor = 10.0**X2Hr
    Cfactor = 10.0**C2Or
    Nfactor = 10.0**C2Or

    # OLD linear method, with max limit at 2.84 dex
    # if Xfactor >= 1.0 / nXsolar:
    #    nHplusHeold = 1e-16
    #    nXold = 1.0
    # else:
    #    nHplusHeold = 1.0 - Xfactor * nXsolar
    #    nXold = Xfactor * nXsolar

    # NEW method defined such that X2Hr is the X-to-H ratio
    #  for both methods, nHplusHe + nX totals to 1
    nHplusHe = (1.0 - nXsolar) / (1.0 - nXsolar + Xfactor * nXsolar)
    nX = Xfactor * nXsolar / (1.0 - nXsolar + Xfactor * nXsolar)

    # breakdown nHplusHe into nH and nHe
    # nHold = nHplusHe
    # nHeold = nHplusHe * solar['nHe'] / solar['nH']
    nH = nHplusHe * solar['nH'] / (solar['nHe'] + solar['nH'])
    nH2 = nHplusHe / 2.0
    nHe = nHplusHe * solar['nHe'] / (solar['nHe'] + solar['nH'])

    # print('           nHHe-old, nHHe-new', nHplusHeold, nHplusHe)
    # print('           nH-old, nH-new', nHold, nH)
    # print('           nHe-old, nHe-new', nHeold, nHe)
    # print('           X-old, X-new', nXold, nX)

    pH2 = nH2 * p  # array
    K1 = np.exp(
        (a1 / temp + b1 + c1 * temp + d1 * temp**2 + e1 * temp**3)
        / (RcalpmolpK * temp)
    )
    AH2 = (pH2**2.0) / (2.0 * K1)
    ACpAO = nX / nXsolar / nH * solar['nO'] * (1.0 + Cfactor * solCtO)
    ACtAO = Cfactor * solCtO * (solar['nO'] ** 2) * ((nX / nXsolar / nH) ** 2)
    BCO = ACpAO + AH2 - np.sqrt((ACpAO + AH2) ** 2 - 4.0 * ACtAO)
    # <--
    # GMR: We should be prepared for T-P profile change that later
    # -->
    nCO = np.mean(BCO * pH2 / p)
    nCH4 = np.mean((2.0 * nX / nXsolar / nH * solar['nC'] - BCO) * (pH2 / p))
    nH2O = np.mean((2.0 * nX / nXsolar / nH * solar['nO'] - BCO) * (pH2 / p))

    a2 = 8.16413e5
    b2 = -2.9109e4
    c2 = 58.5878
    d2 = -7.8284e-4
    e2 = 4.729048e-8
    K2 = np.exp(
        (a2 / temp + b2 + c2 * temp + d2 * temp**2 + e2 * temp**3)
        / (RcalpmolpK * temp)
    )
    AN = nX / nXsolar * Nfactor * solNtO * solar['nO'] / nH  # solar['nN']/nH
    AH2 = (pH2**2.0) / (8.0 * K2)
    # take absolute value, just to be sure that there's never a negative value
    BN2 = AN + AH2 - np.sqrt(np.abs((AN + AH2) ** 2.0 - (AN) ** 2.0))
    BNH3 = 2.0 * (AN - BN2)
    nN2 = np.nanmean(BN2 * pH2 / p)
    nNH3 = np.nanmean(BNH3 * pH2 / p)

    nCO = np.max([nCO, 1e-16])
    nCH4 = np.max([nCH4, 1e-16])
    nH2O = np.max([nH2O, 1e-16])
    nN2 = np.max([nN2, 1e-16])
    nNH3 = np.max([nNH3, 1e-16])

    # make sure that the sum of all mixing ratios does not exceed 1
    #  ok actually it's working great. no problems here (molecules are 60-75%)
    # nMolecules = nH2O + nCH4 + nNH3 + nN2 + nCO
    # print('nSummed vs nMax',nMolecules,1 - nHplusHe)
    # print('nSummed vs nMax',nMolecules / (1 - nHplusHe))

    mixratio = {
        'H2O': np.log10(nH2O) + 6.0,
        'CH4': np.log10(nCH4) + 6.0,
        'NH3': np.log10(nNH3) + 6.0,
        'N2': np.log10(nN2) + 6.0,
        'CO': np.log10(nCO) + 6.0,
    }
    # print('mixratio', mixratio)
    return mixratio, nH2, nHe


# -------------------- -----------------------------------------------
# -- MEAN MOLECULAR WEIGHT -- ----------------------------------------
def getmmw(mixratio, protosolar=True, fH2=None, fHe=None, verbose=False):
    '''
    G. ROUDIER: Mean molecular weight estimate assuming proton mass dominated nucleus
    '''
    molsum = 0.0
    mmw = 0.0
    weights = {
        'H2': 2.0,
        'He': 4.0,
        'CH4': 16.0,
        'NH3': 17.0,
        'H2O': 18.0,
        'H2CO': 30.0,
        'TIO': 64,
        'HCN': 27.0,
        'N2': 28.0,
        'C2H2': 26.0,
        'NO2': 46.0,
        'N2O': 44.0,
        'O3': 48.0,
        'HNO3': 63.0,
        'O2': 32.0,
        'CO': 28.0,
        'CO2': 44.0,
        'NO': 30.0,
        'OH': 17.0,
    }
    for elem in mixratio:
        molsum = molsum + 10.0 ** (mixratio[elem] - 6.0)
        mmw = mmw + 10.0 ** (mixratio[elem] - 6.0) * weights[elem]
        if verbose:
            print('molsum,mmw', molsum, mmw)
    mrH2He = 1.0 - molsum
    # Lodders 2010
    if protosolar:
        HEoH2 = 2.0 * 2.343 * 1e9 / (2.431 * 1e10)
    else:
        HEoH2 = fHe / fH2
    if verbose:
        print('fHe,fH2', fHe, fH2)
    if verbose:
        print('HEoH2', HEoH2)
    mrH2 = mrH2He / (1.0 + HEoH2)
    mrHe = HEoH2 * mrH2
    mmw = mrH2 * weights['H2'] + mrHe * weights['He'] + mmw
    if verbose:
        print('*** mmw ***', mmw)
    return mmw, mrH2, mrHe
