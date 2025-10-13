'''estimators ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-positional-arguments

# -- IMPORTS -- ------------------------------------------------------
import excalibur.system.core as syscore

import numpy as np
import scipy
import math

from astropy.modeling.models import BlackBody

# from astropy import units as u  # CI problem. strange
import astropy.units

from excalibur.ariel.metallicity import massMetalRelation


# ------------- ------------------------------------------------------
# -- ESTIMATOR PROTOTYPES -- -----------------------------------------
class Estimator:
    '''Estimator ds'''

    def __init__(
        self, name, descr, plot='hist', units=None, scale=None, ref=None
    ):
        '''__init__ ds'''
        self._name = name
        self._descr = descr
        self._plot = plot
        self._units = units
        self._scale = scale
        self._ref = ref
        return

    def name(self):
        '''name ds'''
        return self._name

    def descr(self):
        '''descr ds'''
        return self._descr

    def plot(self):
        '''plot ds'''
        return self._plot

    def units(self):
        '''units ds'''
        return self._units

    def scale(self):
        '''scale ds'''
        return self._scale

    def ref(self):
        '''ref ds'''
        return self._ref


class PlEstimator(Estimator):
    '''PlEstimator ds'''

    def __init__(
        self,
        name,
        descr,
        plot='hist',
        units=None,
        scale=None,
        method=None,
        ref=None,
    ):
        '''__init__ ds'''
        Estimator.__init__(self, name, descr, plot, units, scale, ref)
        self._method = method

    def run(self, priors, ests, pl):
        '''run ds'''
        return self._method(priors, ests, pl)


class StEstimator(Estimator):
    '''StEstimator ds'''

    def __init__(
        self,
        name,
        descr,
        plot='hist',
        units=None,
        scale=None,
        method=None,
        ref=None,
    ):
        '''__init__ ds'''
        Estimator.__init__(self, name, descr, plot, units, scale, ref)
        self._method = method

    def run(self, priors, ests):
        '''run ds'''
        return self._method(priors, ests)


# -------------------------- -----------------------------------------
# -- COLLECTION OF ESTIMATORS ----------------------------------------


class TeqEstimator(PlEstimator):
    '''TeqEstimator ds'''

    def __init__(self):
        '''__init__ ds'''
        PlEstimator.__init__(
            self,
            name='Teq',
            descr='Equilibrium temperature',
            units='K',
            ref='albedo = 0',
        )

    def run(self, priors, ests, pl):
        '''run ds'''
        sscmks = syscore.ssconstants(mks=True)

        if 'sma' not in priors[pl].keys():
            return 'missing semi-major axis'
        if priors[pl]['sma'] == '':
            return 'missing semi-major axis'
        if priors['R*'] == '':
            return 'missing R*'
        if priors['T*'] == '':
            return 'missing T*'

        eqtemp = priors['T*'] * np.sqrt(
            priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
        )

        return round(eqtemp, 3)


class HEstimator(PlEstimator):
    '''HEstimator ds'''

    def __init__(self):
        '''__init__ ds'''
        PlEstimator.__init__(
            self,
            name='H',
            descr='Atmospheric scale height (CBE)',
            units='km',
            ref='Chachan metallicity',
        )

    def run(self, priors, ests, pl):
        '''run ds'''
        sscmks = syscore.ssconstants(cgs=True)

        if 'sma' not in priors[pl].keys():
            return 'missing semi-major axis'
        if priors[pl]['sma'] == '':
            return 'missing semi-major axis'
        if priors['R*'] == '':
            return 'missing R*'
        if priors['T*'] == '':
            return 'missing T*'

        eqtemp = priors['T*'] * np.sqrt(
            priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
        )

        if priors[pl]['mass'] == '':  # abort for targets without mass estimates
            return None
        if priors[pl]['logg'] == '':
            g = (
                sscmks['G']
                * priors[pl]['mass']
                * sscmks['Mjup']
                / (priors[pl]['rp'] * sscmks['Rjup']) ** 2
            )
        else:
            g = 10.0 ** float(priors[pl]['logg'])

        mmw = pl_mmw_chachan(priors, ests, pl)

        H = sscmks['Rgas'] * eqtemp / mmw / g

        # convert cm to km
        H = H / 1.0e5

        return round(H, 3)


class HmaxEstimator(PlEstimator):
    '''HmaxEstimator ds'''

    def __init__(self):
        '''__init__ ds'''
        PlEstimator.__init__(
            self,
            name='H_max',
            descr='Atmospheric scale height (max)',
            units='km',
            ref='solar composition',
        )

    def run(self, priors, ests, pl):
        '''run ds'''
        sscmks = syscore.ssconstants(cgs=True)

        if 'sma' not in priors[pl].keys():
            return 'missing semi-major axis'
        if priors[pl]['sma'] == '':
            return 'missing semi-major axis'
        if priors['R*'] == '':
            return 'missing R*'
        if priors['T*'] == '':
            return 'missing T*'

        eqtemp = priors['T*'] * np.sqrt(
            priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
        )

        if priors[pl]['mass'] == '':  # abort for targets without mass estimates
            return None
        if priors[pl]['logg'] == '':
            g = (
                sscmks['G']
                * priors[pl]['mass']
                * sscmks['Mjup']
                / (priors[pl]['rp'] * sscmks['Rjup']) ** 2
            )
        else:
            g = 10.0 ** float(priors[pl]['logg'])

        mmw = pl_mmwmin(priors, ests, pl)

        H = sscmks['Rgas'] * eqtemp / mmw / g

        # convert cm to km
        H = H / 1.0e5

        return round(H, 3)


def pl_metals_thorngren(priors, _ests, pl):
    '''pl_metals_thorngren ds'''
    if priors[pl]['mass'] == '':  # abort for targets without mass estimates
        # print('no mass for this planet')
        return None

    # logmetStar = 0
    logmetStar = priors['FEH*']

    # 318 Earth masses per Jupiter mass
    # pivot point (where metallicity is max value) is at 10 Earth masses
    # metallicity = 3 - np.log10(318.*priors[pl]['mass'])
    # metallicity = min(metallicity, 2)
    # above is more or less the same as massMetalRelation(thorngren=False)
    # metallicity = massMetalRelation(logmetStar, priors[pl]['mass'], thorngren=False)
    # print('Mp,metallicity old',priors[pl]['mass'],metallicity)

    metallicity = massMetalRelation(
        logmetStar, priors[pl]['mass'],
        thorngren=True,
    )
    # print('Mp,metallicity new',priors[pl]['mass'],metallicity)

    return round(metallicity, 4)


def pl_metals_chachan(priors, _ests, pl):
    '''pl_metals_chachan ds'''
    if priors[pl]['mass'] == '':  # abort for targets without mass estimates
        # print('no mass for this planet')
        return None

    logmetStar = priors['FEH*']

    metallicity = massMetalRelation(
        logmetStar, priors[pl]['mass'],
        chachan=True,
    )
    # print('Mp,metallicity new',priors[pl]['mass'],metallicity)

    return round(metallicity, 4)


def pl_mmw_thorngren(priors, _ests, pl):
    '''pl_mmw ds'''
    metallicity = pl_metals_thorngren(priors, _ests, pl)
    if metallicity is None:
        return None

    (a, b, c) = (2.274, 0.02671737, 2.195719)
    mmw = a + b * np.exp(c * metallicity)
    return round(mmw, 5)


def pl_mmw_chachan(priors, _ests, pl):
    '''pl_mmw ds'''
    metallicity = pl_metals_chachan(priors, _ests, pl)
    if metallicity is None:
        return None

    (a, b, c) = (2.274, 0.02671737, 2.195719)
    mmw = a + b * np.exp(c * metallicity)
    return round(mmw, 5)


def pl_mmwmin(_priors, _ests, _pl):
    '''pl_mmwmin ds'''
    return 2.274


def pl_modulation(priors, _ests, pl):
    '''spectral modulation (2 H R_p / R_*^2); assumed mass-metal relation'''

    sscmks = syscore.ssconstants(cgs=True)

    # abort for targets without mass or semi-major axis values
    if priors[pl]['mass'] == '':
        return 'missing planet mass'
    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'
    # abort if there's no stellar radius or temperature
    if priors['R*'] == '':
        return 'missing R*'
    if priors['T*'] == '':
        return 'missing T*'

    eqtemp = priors['T*'] * np.sqrt(
        priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
    )

    # g = sscmks['G'] * priors[pl]['mass']*sscmks['Mjup'] / \
    #             (priors[pl]['rp']*sscmks['Rjup'])**2
    g = 10.0 ** float(priors[pl]['logg'])

    mmw = pl_mmw_chachan(priors, _ests, pl)

    H = sscmks['Rgas'] * eqtemp / mmw / g

    # Spectral modulation is a dimensionless estimate for changes in the transit depth
    #  It is defined as the change in transit depth when increasing the planet size by 1 scale height
    #   = 2 H R_p / R_*^2
    spectral_modulation = (
        2
        * H
        * priors[pl]['rp']
        / priors['R*'] ** 2
        * sscmks['Rjup']
        / sscmks['Rsun'] ** 2
    )
    # print('spectral modulation',spectral_modulation)
    return spectral_modulation


def pl_modulationmax(priors, _ests, pl):
    '''spectral modulation (2 H R_p / R_*^2); solar composition'''

    sscmks = syscore.ssconstants(cgs=True)

    # abort for targets without mass or semi-major axis values
    if priors[pl]['mass'] == '':
        return 'missing planet mass'
    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'
    # abort if there's no stellar radius or temperature
    if priors['R*'] == '':
        return 'missing R*'
    if priors['T*'] == '':
        return 'missing T*'

    eqtemp = priors['T*'] * np.sqrt(
        priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
    )

    # g = sscmks['G'] * priors[pl]['mass']*sscmks['Mjup'] / \
    #             (priors[pl]['rp']*sscmks['Rjup'])**2
    g = 10.0 ** float(priors[pl]['logg'])

    # H/He-dominant atmosphere (for minimum mmw case)
    mmw_min = pl_mmwmin(priors, _ests, pl)

    H_max = sscmks['Rgas'] * eqtemp / mmw_min / g

    # Spectral modulation is a dimensionless estimate for changes in the transit depth
    #  It is defined as the change in transit depth when increasing the planet size by 1 scale height
    #   = 2 H R_p / R_*^2
    spectral_modulation_max = (
        2
        * H_max
        * priors[pl]['rp']
        / priors['R*'] ** 2
        * sscmks['Rjup']
        / sscmks['Rsun'] ** 2
    )
    # print('spectral modulation max',spectral_modulation_max)
    return spectral_modulation_max


def pl_ZFOM(priors, _ests, pl):
    '''pl_ZFOM'''

    error_message = ''
    if priors[pl]['mass'] == '':
        error_message = 'missing planet mass'
    if 'sma' not in priors[pl].keys():
        error_message = 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        error_message = 'missing semi-major axis'
    if priors['R*'] == '':
        error_message = 'missing R*'
    if priors['T*'] == '':
        error_message = 'missing T*'
    if not error_message == '':
        return error_message

    if 'Hmag' in priors.keys():
        Hmag = priors['Hmag']

        if Hmag == '':
            return 'blank Hmag'
    else:
        return 'no Hmag'

    # calculate Zellem Figure-of-Merit (Zellem 2017, Eq.10)
    modulation = pl_modulation(priors, _ests, pl)
    ZFOM = modulation / 10 ** (Hmag / 5)

    # normalization to make it scaled similar to TSM
    return round(ZFOM * 1.0e8, 5)


def pl_ZFOMmax(priors, _ests, pl):
    '''pl_ZFOMmax ds'''

    error_message = ''
    if priors[pl]['mass'] == '':
        error_message = 'missing planet mass'
    if 'sma' not in priors[pl].keys():
        error_message = 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        error_message = 'missing semi-major axis'
    if priors['R*'] == '':
        error_message = 'missing R*'
    if priors['T*'] == '':
        error_message = 'missing T*'
    if not error_message == '':
        return error_message

    if 'Hmag' in priors.keys():
        Hmag = priors['Hmag']

        if Hmag == '':
            return 'blank Hmag'
    else:
        return 'no Hmag'

    # calculate Zellem Figure-of-Merit (Zellem 2017, Eq.10)
    modulation = pl_modulationmax(priors, _ests, pl)
    ZFOMmax = modulation / 10 ** (Hmag / 5)
    # print('ZFOMmax',ZFOMmax*1.e8)

    # normalization to make it scaled similar to TSM
    return round(ZFOMmax * 1.0e8, 5)


def pl_density(priors, _ests, pl):
    '''pl_density ds'''
    # abort for targets without mass estimates
    if priors[pl]['mass'] == '':
        return 'missing planet mass'
    if priors[pl]['rp'] == '':
        return 'missing planet radius'

    sscmks = syscore.ssconstants(cgs=True)
    volume = (4.0 / 3) * math.pi * priors[pl]['rp'] ** 3
    density = priors[pl]['mass'] / volume
    # now convert from Jupiter masses per Jupiter radii to g/cm^3
    conversion = sscmks['Mjup'] / (sscmks['Rjup'] ** 3)
    density = density * conversion

    return round(density, 5)


def st_luminosity(priors, _ests):
    '''st_luminosity ds'''
    sscmks = syscore.ssconstants(mks=True)
    Tsun = sscmks['Tsun']
    if 'L*' in priors.keys() and priors['L*'] != '':
        # est = 10.**float(priors['L*'])
        # autofill now corrects for luminosity being logged in the Archive
        est = priors['L*']
        # est_ref = priors['L*_ref']
        if priors['R*'] != '' and priors['T*'] != '':
            # print('L* vs T*R* consistency?',
            #      est/(priors['R*']**2*(priors['T*']/Tsun)**4),
            #      est, priors['R*']**2*(priors['T*']/Tsun)**4)
            pass
    else:
        if priors['R*'] == '':
            est = 'missing R*'
        elif priors['T*'] == '':
            est = 'missing T*'
        else:
            est = priors['R*'] ** 2 * (priors['T*'] / Tsun) ** 4
            # est_ref = 'from R_star & T_star'
            est = round(est, 5)
    return est  # ,est_ref


def st_spTyp(priors, _ests):
    '''st_spTyp ds'''
    if 'spTyp' in priors.keys():
        est = priors['spTyp']
    else:
        est = 'N/A'
    return est


def pl_insolation(priors, ests, pl):
    '''pl_insolation ds'''
    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'

    # if there's an error message for luminosity, just pass it along (e.g. 'missing R_*')
    if isinstance(ests['luminosity'], str):
        return ests['luminosity']

    insolation = ests['luminosity'] * (priors[pl]['sma'] ** -2)

    # apply eccentricity correction
    insolation *= math.sqrt(1 / (1 - priors[pl]['ecc'] ** 2))

    return round(insolation, 4)


def st_rotationPeriod(priors, _ests):
    '''st_rotationPeriod ds'''
    if 'AGE*' in priors.keys():
        Age = priors['AGE*']
    else:
        Age = 5
        # *** important assumption - arbitrary guess for the age if it's blank ***
    if Age == '':
        Age = 5  # Gyr

    if 'T*' not in priors.keys():
        return 'missing T*'
    if priors['T*'] == '':
        return 'missing T*'

    if priors['T*'] < 3500:
        # rotation period for M2.5-6.0 V stars
        # applies for GJ 1132, GJ 3053, K2-18, GJ 1214
        if Age < 0.012 + 0.061:
            Prot = 1
        else:
            Prot = (Age - 0.012) / 0.061
    elif priors['T*'] < 4000:
        # rotation period for M0 V stars
        # applies for Kepler 138, K2-3
        if Age < 0.365 + 0.019:
            Prot = 1
        else:
            Prot = ((Age - 0.365) / 0.019) ** (1 / 1.457)
    else:
        # rotation period for K-type stars
        # applies for 55 Cnc, GJ 9827, GJ 97658
        Prot = 25 * (Age / 4.6) ** 0.5

    return round(Prot, 4)


def st_coronalTemp(priors, _ests):
    '''st_coronalTemp ds'''
    Prot = st_rotationPeriod(priors, _ests)

    if isinstance(Prot, str):
        return 'missing P_rot'

    # this formula needs work; it is not continuous at the 18.9-day break
    if Prot > 18.9:
        Tcorona = 1.5e6 * (27 / Prot) ** 1.2
    else:
        Tcorona = 1.98e6 * (27 / Prot) ** 0.37

    return Tcorona


def pl_windVelocity(priors, ests, pl):
    '''pl_windVelocity ds'''
    sscmks = syscore.ssconstants(cgs=True)

    if 'M*' not in priors.keys():
        if priors['LOGG*'] == '':
            return None
        priors['M*'] = (
            10.0 ** float(priors['LOGG*'])
            / sscmks['G']
            / sscmks['Msun']
            * (priors['R*'] * sscmks['Rsun']) ** 2
        )
        pass

    if priors['M*'] == '':
        if priors['LOGG*'] == '':
            return None
        priors['M*'] = (
            10.0 ** float(priors['LOGG*'])
            / sscmks['G']
            / sscmks['Msun']
            * (priors['R*'] * sscmks['Rsun']) ** 2
        )
        pass

    mmw = 0.67

    Tcorona = st_coronalTemp(priors, ests)

    # Tcorona can be undefined, e.g. 'missing P_rot'
    if isinstance(Tcorona, str):
        return Tcorona

    v_crit = np.sqrt(sscmks['Rgas'] * Tcorona / mmw)

    r_crit = sscmks['G'] * sscmks['Msun'] * priors['M*'] / 2 / v_crit**2

    # Parker solar wind solution: v^2 - ln(v^2) = 4(1/r + ln(r)) - 3
    #   comes from  dv/dr(v - 1/v) = 2(1/r - 1/r^2)
    # (r,v are normalized by r_crit,v_crit)

    def parkerSolution(v, r):
        '''parkerSolution ds'''
        return v**2 - np.log(v**2) - 4 / r - 4 * np.log(r) + 3

    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'

    r = priors[pl]['sma'] * sscmks['AU'] / r_crit
    # (the second field is an initial guess for v)
    v = scipy.optimize.root(parkerSolution, r, args=r)['x'][0]

    # convert cm/s to km/s
    windVelocity = v * v_crit / 1.0e5
    return round(windVelocity, 4)


def pl_windDensity(priors, ests, pl):
    '''pl_windDensity ds'''
    sscmks = syscore.ssconstants(cgs=True)

    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'
    r = priors[pl]['sma'] * sscmks['AU'] / sscmks['Rsun']  # RSun units

    numberDensity = 3.3e5 / r**2 + 4.1e6 / r**4 + 8.0e7 / r**6

    mmw = 0.67
    m_H = 1.6735e-24
    windDensity = numberDensity * m_H * mmw

    # adjustment based on rotation.  fast spin fives higher density
    Prot = st_rotationPeriod(priors, ests)
    if isinstance(Prot, str):
        return 'missing P_rot'

    windDensity *= (18.9 / Prot) ** 0.6
    return windDensity


def pl_windMassLoss(priors, ests, pl):
    '''pl_windMassLoss ds'''
    sscmks = syscore.ssconstants(cgs=True)

    entrainmentEfficiency = 0.3

    windVelocity = pl_windVelocity(priors, ests, pl)
    if isinstance(windVelocity, str):
        return 'missing wind velocity'
    if windVelocity is None:
        return 'missing wind velocity'

    # convert velocity from km/s to cm/s
    windVelocity *= 1.0e5
    windDensity = pl_windDensity(priors, ests, pl)

    massLossRate = (
        2
        * np.pi
        * entrainmentEfficiency
        * (priors[pl]['rp'] * sscmks['Rjup']) ** 2
        * windDensity
        * windVelocity
    )

    # convert to Jupiter masses per gigayear
    mdotJupUnits = massLossRate / sscmks['Mjup'] * 3.16e7 * 1.0e9
    return mdotJupUnits


# X-ray flux based on stellar age and spectral type (T_* actually)
def Lxray(priors):
    '''Lxray ds'''
    if 'AGE*' in priors.keys():
        Age = priors['AGE*']
    else:
        Age = 5
        # *** important assumption - arbitrary guess for the age if it's blank ***
    if Age == '':
        Age = 5  # Gyr
    logAge = np.log10(Age * 1.0e9)

    # these approx formulae are from by-eye fits to from Raissa's notebook figures
    st_temp = priors['T*']
    if st_temp in ('', 0):
        return 'missing T*'

    if st_temp < 3900:  # M stars
        L_X = 10.0 ** (26.5 + 1.3 * (10.0 - logAge))
    elif st_temp < 5300:  # K stars
        L_X = 10.0 ** (28.0 + 0.7 * (10.0 - logAge))
    else:  # G and F stars
        L_X = 10.0 ** (28.3 + 0.5 * (10.0 - logAge))

    return round(L_X, 5)


def pl_evapMassLoss(priors, _ests, pl):
    '''pl_evapMassLoss ds'''
    sscmks = syscore.ssconstants(cgs=True)

    heatingEfficiency = 0.1

    # X-ray flux based on stellar age and spectral type (T_* actually)
    L_X = Lxray(priors)
    if isinstance(L_X, str):
        return 'missing L_X'

    # assumed relationship between X-ray and EUV fluxes
    #  (from Raissa email 7/9/21: 'logFeuv = 2.63 + 0.58*logFx)

    # need to convert luminosity to a surface flux
    stellar_surface_area = (
        4 * np.pi * (priors['R*'] * sscmks['Rsun']) ** 2
    )  # cm2
    surfaceFlux_X = L_X / stellar_surface_area  # erg/s/cm2 = mW/m2

    surfaceFlux_EUV = 426.6 * (surfaceFlux_X) ** 0.58

    L_EUV = surfaceFlux_EUV * stellar_surface_area

    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'

    # F_X = L_X / (4 * np.pi * (priors[pl]['sma']*sscmks['AU'])**2)
    F_EUV = L_EUV / (4 * np.pi * (priors[pl]['sma'] * sscmks['AU']) ** 2)

    # assume the base of the photoevaporative flow is just the planet radius
    R_base = priors[pl]['rp'] * sscmks['Rjup']

    massLossRate = (
        np.pi
        * heatingEfficiency
        * F_EUV
        * R_base**3
        / (sscmks['G'] * priors[pl]['mass'] * sscmks['Mjup'])
    )

    # convert to Jupiter masses per gigayear
    mdotJupUnits = massLossRate / sscmks['Mjup'] * 3.16e7 * 1.0e9

    return mdotJupUnits


def st_COratio(priors, _ests):
    '''st_COratio ds'''
    if 'FEH*' in priors.keys() and priors['FEH*'] != '':
        # from Table 2 in da Silva 2024
        est = 0.012 + 0.175 * priors['FEH*']
        est = round(est, 5)
    else:
        est = 'N/A'
    return est


def st_NOratio(priors, _ests):
    '''st_NOratio ds'''
    if 'FEH*' in priors.keys() and priors['FEH*'] != '':
        # from Table 2 in da Silva 2024
        est = -0.068 + 0.578 * priors['FEH*']
        est = round(est, 5)
    else:
        est = 'N/A'
    return est


def st_COratio_Nissen(priors, _ests):
    '''st_COratio_Nissen ds'''
    if 'FEH*' in priors.keys() and priors['FEH*'] != '':
        # this is equation 2 from Nissen 2013
        est = -0.002 + 0.22 * priors['FEH*']
        est = round(est, 5)
    else:
        est = 'N/A'
    return est


def pl_beta_rad(priors, _ests, pl):
    '''Beta = radiation pressure / gravity'''
    sscmks = syscore.ssconstants(cgs=True)

    # K2-3 doesn't have a planet radius! (for all 3 planets)
    if isinstance(priors[pl]['rp'], str):
        return 'missing planet radius'

    # (this is actually acceleration, not force)
    Fgrav = (
        sscmks['G']
        * priors[pl]['mass']
        * sscmks['Mjup']
        / (priors[pl]['rp'] * sscmks['Rjup']) ** 2
    )

    Qrad = 1  # assume perfect absorbtion
    rho_grain = 1  # assume density of water
    grainsize_um = 0.1
    grainsize_cm = grainsize_um / 1.0e4
    mass2area = 4 / 3 * rho_grain * grainsize_cm

    Lstar = st_luminosity(priors, _ests)

    if isinstance(Lstar, str):
        return Lstar
    if 'sma' not in priors[pl].keys():
        return 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        return 'missing semi-major axis'

    # (and this is actually pressure (force per area))
    Frad = (
        Qrad
        * Lstar
        * sscmks['Lsun']
        / sscmks['c']
        / (4 * np.pi * (priors[pl]['sma'] * sscmks['AU']) ** 2)
    )

    # (include the grain mass-to-area ratio here)
    Beta = Frad / Fgrav / mass2area
    # print('Beta = ',Beta)

    # there was a detailed comparison here against the Owens formula
    # but apparently a commented out block of code isn't kosher
    # it won't pass the moronic PEP3/dawgie review,
    # even after several changes to accomodate the requirements
    # so instead there's some comments here talking about stupid PEP3/dawgie crap

    # and don't forget the delete the whitespace on the blank lines here!

    return round(Beta, 6)


def pl_TSM(priors, _ests, pl):
    '''pl_TSM'''

    sscmks = syscore.ssconstants(cgs=True)

    error_message = ''
    if priors[pl]['mass'] == '':
        error_message = 'missing planet mass'
    if 'sma' not in priors[pl].keys():
        error_message = 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        error_message = 'missing semi-major axis'
    if priors['R*'] == '':
        error_message = 'missing R*'
    if priors['T*'] == '':
        error_message = 'missing T*'
    if not error_message == '':
        return error_message
    if 'Jmag' in priors.keys():
        Jmag = priors['Jmag']
        if Jmag == '':
            return 'blank Jmag'
    elif 'Hmag' in priors.keys():
        # print(' NOTE: using Hmag instead of Jmag for TSM (temp fix)')
        Jmag = priors['Hmag']
        if Jmag == '':
            return 'blank Jmag'
    else:
        return 'no Jmag'

    # calculate TSM figure-of-merit (Kempton 2018, Eq.1)

    eqtemp = priors['T*'] * np.sqrt(
        priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
    )

    Rp_earth = priors[pl]['rp'] * sscmks['Rjup'] / sscmks['Rearth']
    Mp_earth = priors[pl]['mass'] * sscmks['Mjup'] / sscmks['Mearth']
    Rstar_sun = priors['R*']

    if Rp_earth < 1.5:
        scaleFactor = 0.190
    elif Rp_earth < 2.75:
        scaleFactor = 1.26
    elif Rp_earth < 4.0:
        scaleFactor = 1.28
    else:
        scaleFactor = 1.15
    # scaleFactor = 1
    # print('scaleFactor',scaleFactor)

    TSM = (
        scaleFactor
        * eqtemp
        * Rp_earth**3
        / Rstar_sun**2
        / Mp_earth
        / 10 ** (Jmag / 5)
    )
    # print('TSM',TSM)
    # zfommax = pl_ZFOMmax(priors, _ests, pl)
    # print('zfommax',zfommax)
    # print('zellem/kempton ratio',zfommax / TSM)

    return round(TSM, 4)


def pl_ESM(priors, _ests, pl):
    '''pl_ESM'''

    sscmks = syscore.ssconstants(cgs=True)

    error_message = ''
    if 'sma' not in priors[pl].keys():
        error_message = 'missing semi-major axis'
    if priors[pl]['sma'] == '':
        error_message = 'missing semi-major axis'
    if priors['R*'] == '':
        error_message = 'missing R*'
    if priors['T*'] == '':
        error_message = 'missing T*'
    if not error_message == '':
        return error_message

    if 'Kmag' in priors.keys():
        Kmag = priors['Kmag']
        if Kmag == '':
            return 'blank Kmag'
    elif 'Hmag' in priors.keys():
        print(' NOTE: using Hmag instead of Kmag for ESM (temp fix)')
        Kmag = priors['Hmag']
        if Kmag == '':
            return 'blank Kmag'
    else:
        return 'no Kmag'

    # calculate ESM figure-of-merit (Kempton 2018, Eq.4)

    Tstar = priors['T*']
    eqtemp = Tstar * np.sqrt(
        priors['R*'] * sscmks['Rsun/AU'] / (2.0 * priors[pl]['sma'])
    )
    Tday = 1.1 * eqtemp

    Rp_cgs = priors[pl]['rp'] * sscmks['Rjup']
    Rstar_cgs = priors['R*'] * sscmks['Rsun']

    scaleFactor = 4.29e6

    # BBplanet = BlackBody(temperature=Tday*u.K)
    # BBstar = BlackBody(temperature=Tstar*u.K)
    BBplanet = BlackBody(temperature=Tday * astropy.units.K)
    BBstar = BlackBody(temperature=Tstar * astropy.units.K)

    # using u for units doesn't work because CI is dumb I guess
    #  it doesn't think that there is such a thing as micron.  strange
    # ESM = scaleFactor * BBplanet(7.5*u.micron) / BBstar(7.5*u.micron) \
    #    * (Rp_cgs / Rstar_cgs)**2 \
    #    / 10**(Kmag/5)
    ESM = (
        scaleFactor
        * BBplanet(7.5 * astropy.units.micron)
        / BBstar(7.5 * astropy.units.micron)
        * (Rp_cgs / Rstar_cgs) ** 2
        / 10 ** (Kmag / 5)
    )
    # print('ESM factor1',scaleFactor)
    # print('ESM factor2',BBplanet(7.5*u.micron) / BBstar(7.5*u.micron))
    # print('ESM factor3',Rp_cgs**2 / Rstar_cgs**2)
    # print('ESM factor4',1/ 10**(Kmag/5))
    # print('ESM',ESM)

    return round(ESM, 4)


# --------------------------- ----------------------------------------
