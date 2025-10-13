'''ariel planet_metallicity ds'''

# Heritage code shame:
# pylint: disable=invalid-name

# -- IMPORTS -- ------------------------------------------------------
import excalibur.system.core as syscore
import numpy as np
import logging

log = logging.getLogger(__name__)

# np.random.seed(123)

# ______________________________________________________


def massMetalRelationDisp(logmetStar, Mp,
                          thorngren=False, chachan=False, dispersion=0.3):
    '''
    Add some realistic scatter to the mass-metallicity relation
    (not that we know reality)
    '''
    logmet = massMetalRelation(logmetStar, Mp, thorngren=thorngren, chachan=chachan)

    logmet += np.random.normal(scale=dispersion)

    return logmet


# ______________________________________________________


def massMetalRelation(logmetStar, Mp, thorngren=False, chachan=False):
    '''
    Calculate an assumed planet metallicity based on its mass
     default option: FINESSE linear relationship
     thorngren option: Thorngren 2016 relationship
     chachan option: Chachan 2025 relationship
    Planet mass (Mp) is in Jupiter masses
    '''

    if logmetStar == '':
        log.warning(
            '--< Star metallicity missing in Ariel-sim : add to overwriter.py >--'
        )
        logmetStar = 0

    if chachan:
        # mass-metallicity relation from Chachan et al 2025
        Mcore = 14.7  # -1.6+1.8 Earth masses
        fZ = 0.09  # +-0.01

        Zsun = 0.014

        sscmks = syscore.ssconstants(cgs=True)
        Mp_Earths = Mp * sscmks['Mjup'] / sscmks['Mearth']
        if Mp_Earths < Mcore:
            M_Z = Mp_Earths
        else:
            M_Z = Mcore + fZ * (Mp_Earths - Mcore)
        # this isn't quite right. should be ratio with hydrogen
        #  (but then it would go to infinite)
        Zplanet = M_Z / Mp_Earths

        # ignores Zstar!?
        logmet = np.log10(Zplanet / Zsun)

    elif thorngren:
        # mass-metallicity relation from Thorngren et al 2016
        slope = -0.45
        # metallicity for Jupiter mass (trend value; Jupiter itself is lower)
        intercept = np.log10(9.7)

        logmet = logmetStar + intercept + slope * np.log10(Mp)

    else:
        #   mass-metallicity relation from FINESSE proposal (Fortney motivated)
        # Assume an inverse-linear relationship between planet mass and metallicity
        # Include a limit on metallicity of +2.0 dex, relative to the parent star
        #  Mp = 1Jup gives met = 0.5 dex
        #  Mp = 10/318 = 10Earths gives met = +2 dex
        #  Mp = 1/318 = 1Earth would give met = +3 dex, but capped at +2 dex

        slope = -1.0
        maxMetal = 2.0
        # Mpivot = 1.   # (earth units)
        # intercept = maxMetal - slope*Mpivot
        # Mpivot = -1.5  # (jupiter units)
        intercept = 0.5  # metallicity for Jupiter mass

        logmet = logmetStar + intercept + slope * np.log10(Mp)

        # change so that it can handle an array of masses (from cerberus/plotting)
        if isinstance(Mp, float):
            logmet = min(maxMetal, logmet)
        else:
            logmet[np.where(logmet > maxMetal)] = maxMetal

    return logmet


# ______________________________________________________


def randomStarMetal():
    '''
    If there's no stellar metallicity, pull from random distribution
    '''

    # this is from my excel check of Hinkel's Hypatia catalog
    #  logmetStar=0.06 + 0.29*random.gauss(0.,1.)
    # from Kepler-detection-based (Buchhave 2011)
    logmetStar = -0.01 + 0.25 * np.random.normal()

    return logmetStar


# ______________________________________________________


def randomCtoO_linear(logCtoOaverage=-0.26, logCtoOdispersion=0.3):
    '''
    Assign a random C-to-O ratio to each system
    Allow a small fraction (~5%) of stars to have more C than O
    Actually that's too much.  Consensus at the May2023 JWST conference was less than that I think
    Let's go with -0.2+-0.1, which gives 2.3% with more C than O
    March 2024 update (allowing planets to have more variety than stars):
     - use solar C/O as median
     - have a stdev of 0.3 dex
    '''

    # this is from my excel check of Hinkel's Hypatia catalog
    # logCtoO=-0.13 + 0.19*random.gauss(0.,1.)

    # this is motivated by Fortney's argument that C>O is rare
    # he gets 1-5% and suggests even lower
    # let's just stick with 5% C/O>1,
    # in order to see how well FINESSE can do with some oddballs
    # this gives 4.8% with C>O
    # logCtoO=-0.2 + 0.12*np.random.normal()
    # logCtoO=-0.2 + 0.1*np.random.normal()

    # logCtoO_solar = -0.26  # solar C/O is 0.55

    logCtoO = logCtoOaverage + logCtoOdispersion * np.random.normal()

    CtoO = 10.0**logCtoO
    return CtoO


# ______________________________________________________
