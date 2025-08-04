'''system consistency ds'''

# Heritage code shame:
# pylint: disable=invalid-name

# -- IMPORTS -- -----------------------------------------------------
import numpy as np
import excalibur.system.core as syscore
from excalibur.util.logg import calculate_logg
import logging

log = logging.getLogger(__name__)

# ----------------- -------------------------------------------------
# --       Consistency checks between related parameters           --
# -------------------------------------------------------------------


def consistency_checks(priors, ignoredPlanets):
    '''
    Run through a series of data validation checks, for internal self-consistency
    '''

    verbose = False
    # useful option for checking consistency in test target values
    # verbose = True

    inconsistencies = []

    ok = consistency_check_M_R_LOGG_star(priors, verbose=verbose)
    if not ok:
        inconsistencies.append('LOGG*')

    ok = consistency_check_M_R_RHO_star(priors, verbose=verbose)
    if not ok:
        inconsistencies.append('RHO*')

    ok = consistency_check_R_T_Lstar(priors, verbose=verbose)
    if not ok:
        inconsistencies.append('L*')

    for planet_letter in priors['planets']:
        # special check needed for HAT-P-11 c
        #  (it's a non-transiting RV planet that shouldn't be here, but there's an Archive bug)
        if planet_letter in ignoredPlanets:
            log.info(
                'SKIPPING consistency check for ignored planets: %s',
                planet_letter,
            )
        else:
            ok = consistency_check_sma_P_Mstar(
                priors, planet_letter, verbose=verbose
            )
            if not ok:
                inconsistencies.append(planet_letter + ':semi-major axis')

            # planet density does not exist (so no need to check it)
            # ok = consistency_check_M_R_RHO_planet(priors, planet_letter)
            # if not ok: inconsistencies.append(planet_letter+':rho')

            ok = consistency_check_M_R_LOGG_planet(
                priors, planet_letter, verbose=verbose
            )
            if not ok:
                inconsistencies.append(planet_letter + ':logg')

            ok = consistency_check_Teq_sma_Lstar(
                priors, planet_letter, verbose=verbose
            )
            if not ok:
                inconsistencies.append(planet_letter + ':Teq')

            # impact and inclination should be consistent
            # but impact doesn't exist here; not saved by target
            # ok = consistency_check_inc_impact(priors, planet_letter)
            # if not ok: inconsistencies.append(planet_letter+':impact')

    return inconsistencies


# -------------------------------------------------------------------
def close_to(A, B, eps=1.0e-1):
    '''
    Check whether two values are more or less equal
    (within fractional 'eps' of each other)
    '''
    # print('ccheck  close enough?',A,B,np.abs((A-B)/(A+B)*2))

    if A == 0 and B == 0:
        close_enough = True
    elif np.abs((A - B) / (A + B) * 2) > eps:
        close_enough = False
    else:
        close_enough = True
    return close_enough


# -------------------------------------------------------------------
def consistency_check_M_R_RHO_star(starInfo, verbose=False):
    '''
    Verify that the stellar density is consistent with the stellar mass,radius
    '''
    # get Msun and Rsun definitions, for calculating stellar density from M*,R*
    sscmks = syscore.ssconstants(cgs=True)

    R = starInfo['R*']
    M = starInfo['M*']
    RHO = starInfo['RHO*']
    # print('ccheck R,M,RHO',R,M,RHO)

    consistent = True
    if RHO == '' or R == '' or M == '':
        # print('ccheck PASS: one of the M/R/RHO* fields is missing')
        pass
    else:
        RHOcheck = (
            float(M)
            * sscmks['Msun']
            / (4.0 * np.pi / 3.0 * (float(R) * sscmks['Rsun']) ** 3)
        )
        # print('ccheck RHOcheck',RHOcheck)
        if not close_to(float(RHO), RHOcheck):
            consistent = False
            if verbose:
                print('density value vs should be', RHO, RHOcheck)

    return consistent


# -------------------------------------------------------------------
def consistency_check_M_R_RHO_planet(starInfo, planet_letter, verbose=False):
    '''
    Verify that the stellar density is consistent with the stellar mass,radius
    '''
    # get Mjup and Rjup definitions, for calculating planet density from Mp,Rp
    sscmks = syscore.ssconstants(cgs=True)

    R = starInfo[planet_letter]['rp']
    M = starInfo[planet_letter]['mass']
    RHO = starInfo[planet_letter]['rho']
    # print('ccheck R,M,RHO planet',R,M,RHO)

    consistent = True
    if RHO == '' or R == '' or M == '':
        # print('ccheck PASS: one of the M/R/RHO planet fields is missing')
        pass
    else:
        RHOcheck = (
            float(M)
            * sscmks['Mjup']
            / (4.0 * np.pi / 3.0 * (float(R) * sscmks['Rjup']) ** 3)
        )
        # print('ccheck RHOcheck',RHOcheck)
        if not close_to(float(RHO), RHOcheck):
            consistent = False
            if verbose:
                print('density value vs should be', RHO, RHOcheck)

    return consistent


# -------------------------------------------------------------------
def consistency_check_M_R_LOGG_star(starInfo, verbose=False):
    '''
    Verify that the stellar log(g) is consistent with the stellar mass,radius
    '''

    # get Msun and Rsun definitions, for calculating stellar density from M*,R*
    sscmks = syscore.ssconstants(cgs=True)

    R = starInfo['R*']
    M = starInfo['M*']
    LOGG = starInfo['LOGG*']
    # print('ccheck R,M,LOGG',R,M,LOGG)

    consistent = True
    if LOGG == '' or R == '' or M == '':
        # print('ccheck PASS: one of the M/R/LOGG* fields is missing')
        pass
    else:
        LOGGcheck = calculate_logg(M, R, sscmks, units='solar')
        # print('ccheck  LOGGcheck',LOGGcheck)
        if not close_to(float(LOGG), LOGGcheck):
            consistent = False
            if verbose:
                print('LOGG* vs check', LOGG, LOGGcheck)

    return consistent


# -------------------------------------------------------------------
def consistency_check_M_R_LOGG_planet(starInfo, planet_letter, verbose=False):
    '''
    Verify that the planetary log(g) is consistent with the planet mass,radius
    '''
    # get Mjup and Rjup definitions, for calculating stellar density from Mp,Rp
    sscmks = syscore.ssconstants(cgs=True)

    R = starInfo[planet_letter]['rp']
    M = starInfo[planet_letter]['mass']
    LOGG = starInfo[planet_letter]['logg']
    # print('ccheck R,M,LOGG planet',R,M,LOGG)

    consistent = True
    if LOGG == '' or R == '' or M == '':
        # print('ccheck PASS: one of the M/R/LOGG planet fields is missing')
        pass
    else:
        LOGGcheck = calculate_logg(M, R, sscmks, units='Jupiter')
        # print('ccheck  LOGGcheck',LOGGcheck)
        if not close_to(float(LOGG), LOGGcheck):
            consistent = False

    if verbose:
        print('logg vs check', LOGG, LOGGcheck)

    return consistent


# -------------------------------------------------------------------
def consistency_check_sma_P_Mstar(starInfo, planet_letter, verbose=False):
    '''
    Verify that the semi-major axis and period are self-consistent
    '''
    # get Msun and Rsun definitions, for calculating stellar density from M*,R*
    sscmks = syscore.ssconstants(cgs=True)

    M = starInfo['M*']
    P = starInfo[planet_letter]['period']
    sma = starInfo[planet_letter]['sma']
    # print('ccheck M P sma',M,P,sma)

    consistent = True
    if M == '' or P == '' or sma == '':
        # print('ccheck PASS: one of the M/P/sma fields is missing')
        pass
    else:
        GM = sscmks['G'] * float(M) * sscmks['Msun']
        smaCheck = (GM * (float(P) * sscmks['day'] / 2.0 / np.pi) ** 2) ** (
            1.0 / 3.0
        )
        smaCheck /= sscmks['AU']
        if not close_to(float(sma), smaCheck):
            consistent = False
            if verbose:
                print('ccheck  smaCheck', smaCheck)

    return consistent


# -------------------------------------------------------------------
def consistency_check_R_T_Lstar(starInfo, verbose=False):
    '''
    Verify that the stellar luminosity matches the star radius,temperature
    '''

    # L and R are kept in solar units, but get Tsun definition here
    sscmks = syscore.ssconstants(cgs=True)

    R = starInfo['R*']
    T = starInfo['T*']
    L = starInfo['L*']
    # print('ccheck R T log-L',R,T,L)

    consistent = True
    if R == '' or T == '' or L == '':
        # print('ccheck PASS: one of the R/T/Lstar fields is missing')
        pass
    else:
        Lcheck = R**2 * (T / sscmks['Tsun']) ** 4  # (solar units)
        if not close_to(float(L), Lcheck):
            consistent = False
            if verbose:
                print('ccheck  L,Lcheck', float(L), Lcheck)

    return consistent


# -------------------------------------------------------------------
def consistency_check_Teq_sma_Lstar(starInfo, planet_letter, verbose=False):
    '''
    Verify that the planet equilibrium temperature matches it's stellar radiation
    '''

    L = starInfo['L*']
    Teq = starInfo[planet_letter]['teq']
    sma = starInfo[planet_letter]['sma']
    # print('ccheck M P sma',M,P,sma)

    consistent = True
    if L == '' or Teq == '' or sma == '':
        # print('ccheck PASS: one of the Teq/L/sma fields is missing')
        pass
    else:
        TeqCheck = 278.3 * float(L) ** 0.25 / float(sma) ** 0.5
        if not close_to(float(Teq), TeqCheck):
            consistent = False
            if verbose:
                print('ccheck  Teq', Teq, TeqCheck)

    return consistent
