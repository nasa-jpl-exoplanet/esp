'''system autofill ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import numpy
import copy
import excalibur.system.constants as syscons
from excalibur.util.logg import calculate_logg
import logging

log = logging.getLogger(__name__)


# ----------------- --------------------------------------------------
# -- FILL BLANK UNCERTAINTY FIELD ------------------------------------
def fillUncertainty(param, param_value, param_uncertainty, error_type):
    '''
    Put in some nominal value for the parameter uncertainty
    '''

    # print()
    # print(param,error_type)
    # print(param_value)
    # print(param_uncertainty)

    autofilled = False
    try:
        fillvalue = float(param_uncertainty)
    except ValueError:
        autofilled = True
        huh = False
        if param == 'spTyp':
            # (spectral type isn't mandatory, so probably never comes here)
            fillvalue = ''
        elif param == 'T*':
            # for stellar temperature, 100 K uncertainty as default
            fillvalue = 100
            fillvalue = 150
            fillvalue = 185  # 90-percentile; 95-percentile is 200
        elif param == 'inc':
            # inclination should be pretty close to 90degrees
            #  let's try 5deg uncertainty
            fillvalue = 5
            fillvalue = 1.93  # 90-percentile; 95-percentile is 2.7
            # make sure inclination range stays within 0-90 degrees
            # actually no, allow inclination to go 0-180, as in the Archive
            # if param_value + fillvalue > 90: fillvalue = 90 - param_value
        elif param == 'ecc':
            # for eccentricity, set uncertainty to 20%, with a minimum of 0.1
            fillvalue = numpy.max([float(param_value) * 2.0e-1, 0.1])
            fillvalue = 0.145  # 90-percentile; 95-percentile is 0.185
        elif param == 'omega':
            # for argument of periastron, set a large uncertainty (120 degrees)
            fillvalue = 120
        elif param == 'FEH*':
            # for metallicity, fallback uncertainty is 0.1 dex
            fillvalue = 0.1
            fillvalue = 0.2
            fillvalue = 0.295  # 90-percentile; 95-percentile is 0.3
        else:
            huh = True
        if huh:
            if param in ['LOGG*', 'logg']:
                # for planet or stellar log-g, fallback uncertainty is 0.1 dex
                fillvalue = 0.1
                fillvalue = 0.2
                fillvalue = 0.3  # 90-percentile; 95-percentile is 0.3
            elif param == 'Hmag' or param == 'Jmag' or param == 'Kmag':
                # 0.1 dex for J,H,K band magnitudes
                # offenders: Kepler-1651, TOI-2583 A, TOI-3976 A, WASP-193
                fillvalue = 0.1
            elif param == 'period':
                # orbital period should be known pretty well. how about 1%?
                fillvalue = float(param_value) * 1.0e-2
                fillvalue = float(param_value) * 1.0e-4
                fillvalue = (
                    float(param_value) * 1.5e-5
                )  # 90-percentile; 95-percentile is 2.7e-5
            elif param == 't0':
                # T_0 is known to much much better than 10%
                #   set it to 10% of the period?
                #   problem though - period is not passed in here
                # for now, set it to 1 hour uncertainty
                fillvalue = 1.0 / 24.0
                fillvalue = 0.01
                fillvalue = 0.0076  # 90-percentile; 95-percentile is 0.01
            elif param in ['rp', 'sma', 'mass', 'R*', 'M*', 'RHO*', 'L*']:
                # set uncertainty to 10% for planet radius, semi-major axis, mass
                #   same for stellar radius, mass, density, luminosity
                fillvalue = float(param_value) * 0.1
                if param in ['rp']:
                    fillvalue = float(param_value) / 3.0
                    fillvalue = (
                        float(param_value) * 0.392
                    )  # 90-percentile; 95-percentile is 0.50
                elif param in ['sma']:
                    fillvalue = float(param_value) / 20.0
                    fillvalue = (
                        float(param_value) * 0.053
                    )  # 90-percentile; 95-percentile is 0.093
                elif param in ['mass']:
                    fillvalue = float(param_value) / 5
                    fillvalue = (
                        float(param_value) * 0.44
                    )  # 90-percentile; 95-percentile is 0.66
                elif param in ['R*']:
                    fillvalue = float(param_value) / 5
                    fillvalue = (
                        float(param_value) * 0.36
                    )  # 90-percentile; 95-percentile is 0.45
            elif param == 'AGE*':
                # age is generally not well known.  have at least 50% uncertainty
                fillvalue = float(param_value) * 0.5
            elif param == 'teq':
                # planet equilibrium temperature to maybe 10%?
                #  (error should really be derived from errors on L*,a_p)
                # a_p error is 5%, so 10% here should be very conservative
                fillvalue = float(param_value) / 10.0
            elif param == 'impact':
                # impact parameter is a number from 0 to 1 (normalized to stellar radius)
                # inclination uncertainty is ~2deg, and rp/sma is ~10, so maybe 0.2 here?
                fillvalue = 0.2
            elif param == 'trandepth':
                # transit depth to 20%? S/N=5 seems like a minimum for detection
                fillvalue = float(param_value) * 0.2
            elif param == 'trandur':
                # transit duration to 20%? S/N=5 seems reasonable I guess
                fillvalue = float(param_value) * 0.2
            elif param == 'dist':
                # 20% for distance, I guess
                # offenders: WASP-103, WASP-105, WASP-20. HATS-58, HATS-12, K2-238
                fillvalue = float(param_value) * 0.2
            else:
                # fallback option is to set uncertainty to 10%
                fillvalue = float(param_value) * 0.1
                print('another PARAM:', param)
                pass
            pass
        # make sure that the upper error is positive and the lower is negative
        if error_type == 'uperr':
            fillvalue = abs(fillvalue)
        elif error_type == 'lowerr':
            fillvalue = -abs(fillvalue)
        else:
            # exit('ERROR: unknown error_type')
            pass
        pass
    return fillvalue, autofilled


# ----------------- --------------------------------------------------
# -- SELECT THE BEST PARAMETER VALUE FROM VARIOUS ARCHIVE VALUES -----
def bestValue(
    valuesorig,
    uperrsorig,
    lowerrsorig,
    refsorig,
    lbl,
    bestref,
    bestpubIndex,
    maximizeSelfConsistency=True,
    selectMostRecent=False,
    verbose=False,
):
    '''
    From a list of parameter values, determine the most trustworthy value
    There are three options here:
    1) use the Exoplanet Archive default
    2) use the most recent publication (for a selected list of parameters, e.g. t0)
    3) use the publication that can fill in needed params self-consistently
    '''

    values = copy.deepcopy(valuesorig)
    uperrs = copy.deepcopy(uperrsorig)
    lowerrs = copy.deepcopy(lowerrsorig)
    refs = copy.deepcopy(refsorig)

    # option to select the most recent publication for period and T0
    if selectMostRecent:
        selectMostRecent = lbl in ('period', 't0')

    if maximizeSelfConsistency:
        # push the best values into the default index ([0])
        # 1) old method. doesn't work for derived values (they have a different ref name)
        for iref, ref in enumerate(refs):
            if ref == bestref:
                values[0] = values[iref]
                uperrs[0] = uperrs[iref]
                lowerrs[0] = lowerrs[iref]
                refs[0] = refs[iref]
        # 2) new method.  use the index pulled from calculate_selfConsistency_metric
        if bestpubIndex != -1:
            if verbose and values[0] != values[bestpubIndex]:
                print('filling in an associated derived value!', lbl)
                print(' old value,ref', values[0], refs[0])
                print(
                    ' new value,ref', values[bestpubIndex], refs[bestpubIndex]
                )
            values[0] = values[bestpubIndex]
            uperrs[0] = uperrs[bestpubIndex]
            lowerrs[0] = lowerrs[bestpubIndex]
            refs[0] = refs[bestpubIndex]
        else:
            if verbose:
                print('Oof! the best publication index is undefined!')
            pass

    if values[0] != '' and not selectMostRecent:
        # step 1: if there is a default value at the start of the list, use that
        # -- exception: for the emphemeris (period and t_0) always use the most recent value --
        bestvalue = values[0]
        bestuperr = uperrs[0]
        bestlowerr = lowerrs[0]
        bestref = refs[0]
    else:
        # step 2: find the most recently published non-blank value
        bestvalue = ''
        bestref = ''
        bestyear = 0
        bestuperr = ''
        bestlowerr = ''
        for value, uperr, lowerr, ref in zip(values, uperrs, lowerrs, refs):
            try:
                year = int(ref[-4:])
            except ValueError:
                # there are some refs without a year; make them lower priority
                year = 1
            # select the most recently published non-blank value
            if value != '' and year >= bestyear:
                # use author name as tie-breaker for 2 with same year
                if year == bestyear and bestref < ref:
                    # print('skipping this same-year ref',year,ref,bestref)
                    pass
                else:
                    bestyear = year
                    bestvalue = value
                    bestuperr = uperr
                    bestlowerr = lowerr
                    bestref = ref

        # if selectMostRecent:
        #     print(lbl,'old default:',values[0],' new most recent:',bestvalue)

    return bestvalue, bestuperr, bestlowerr, bestref


# -------------------------------------------------------------------
def calculate_selfConsistency_metric(data, setRefByHand=False, verbose=False):
    '''
    Determine how many key parameters (currently 12) are provided by each publication
    '''

    # list which parameters should be considered for self-consistency
    starParams = ['R*', 'FEH*', 'T*', 'LOGG*']
    planetParams = ['rp', 'mass', 'sma', 'period', 't0', 'inc', 'ecc', 'omega']

    counts = []
    Npublications = len(data['R*'])

    iplanet = 0
    nplanet = 0
    nplanetsum = 0
    for ipub in range(Npublications):
        counts.append({})
        starRef = data['R*_ref'][ipub]
        counts[ipub]['starref'] = starRef
        counts[ipub]['starrefcount'] = 0
        if starRef == setRefByHand:
            counts[ipub]['starrefcount'] = 100
        for starParam in starParams:
            if starParam + '_ref' not in data.keys():
                print('PROBLEM: missing ref param', starParam + '_ref')
            if data[starParam + '_ref'][ipub] != starRef and not data[
                starParam + '_ref'
            ][ipub].startswith('derived'):
                print(
                    'TROUBLE: different reference for star data',
                    starParam,
                    data[starParam + '_ref'][ipub],
                )
            if data[starParam][ipub] != '':
                counts[ipub]['starrefcount'] += 1

        # it is difficult to find the corresponding planet data
        # the first few indices correspond to the default for each planet
        #  then the rest is groups in batches by planet.  (batch length varies)
        # getting the number of refs for each planet is a good first step
        Nplanets = len(data['planets'])
        Nref = {}
        Nrefsum = 0
        for planet in data['planets']:
            Nref[planet] = len(data[planet]['rp'])
            Nrefsum += Nref[planet]
        if Nrefsum != len(data['R*']):
            log.warning('ERROR: number of star and planet pubs mismatched!')

        if ipub < Nplanets:
            planet = data['planets'][ipub]
            ipubplanet = 0
        else:
            planet = data['planets'][iplanet]
            ipubplanet = ipub - Nplanets + 1 - nplanetsum
            nplanet += 1
            if nplanet >= Nref[planet] - 1:  # jump to the next planet
                nplanetsum += nplanet
                nplanet = 0  # reset the planet publication counter back to zero
                iplanet += 1
        # if verbose: print('ipub cross match',ipub,planet,ipubplanet)
        counts[ipub]['planet'] = planet

        planetRef = data[planet]['rp_ref'][ipubplanet]
        counts[ipub]['planetref'] = planetRef
        counts[ipub]['planetrefcount'] = 0
        if planetRef == setRefByHand:
            counts[ipub]['planetrefcount'] = 100
        if verbose:
            if starRef == planetRef:
                print(
                    'starRef vs planetRef SAME',
                    planet,
                    ipub,
                    ipubplanet,
                    starRef,
                    planetRef,
                )
            else:
                print(
                    'starRef vs planetRef DIFF',
                    planet,
                    ipub,
                    ipubplanet,
                    starRef,
                    planetRef,
                )
        for planetParam in planetParams:
            if planetParam + '_ref' not in data[planet].keys():
                print('PROBLEM: missing ref param', planetParam + '_ref')
            if data[planet][planetParam + '_ref'][
                ipubplanet
            ] != planetRef and not data[planet][planetParam + '_ref'][
                ipubplanet
            ].startswith(
                'derived'
            ):
                print(
                    'TROUBLE: different reference for planet data',
                    planetParam,
                    data[planet][planetParam + '_ref'][ipub],
                )
            if data[planet][planetParam][ipubplanet] != '':
                counts[ipub]['planetrefcount'] += 1

        # sum up the star and planet counts (if it's the same publication)
        if counts[ipub]['planetref'] == counts[ipub]['starref']:
            counts[ipub]['totalcount'] = (
                counts[ipub]['starrefcount'] + counts[ipub]['planetrefcount']
            )
            counts[ipub]['totalref'] = counts[ipub]['starref']
        elif counts[ipub]['planetrefcount'] > counts[ipub]['starrefcount']:
            counts[ipub]['totalcount'] = counts[ipub]['planetrefcount']
            counts[ipub]['totalref'] = counts[ipub]['planetref']
        else:
            counts[ipub]['totalcount'] = counts[ipub]['starrefcount']
            counts[ipub]['totalref'] = counts[ipub]['starref']

    # if verbose:
    #    for count in counts: print('counts',count)

    totalcounts = numpy.array([count['totalcount'] for count in counts])
    bestscore = numpy.max(totalcounts)
    refs = numpy.array([count['totalref'] for count in counts])
    refs_with_max_score = refs[numpy.where(totalcounts == bestscore)]
    if verbose:
        print('refs_with_max_score', refs_with_max_score)
    bestyear = 0
    bestref = 'ZZ'
    for ref in refs_with_max_score:
        try:
            year = int(ref[-4:])
        except ValueError:
            # there are some refs without a year; make them lower priority
            year = 1
        # select the most recently published as tiebreaker
        if year >= bestyear:
            # use author name as tie-breaker for 2 with same year
            if year == bestyear and bestref < ref:
                pass
            else:
                bestyear = year
                bestref = ref
                # ref_with_max_score = ref

    # alternate metric: combine total count for all planets, not just individual star+planet lines
    bestref_oldmethod = bestref
    for ipub in range(Npublications):
        counts[ipub]['totalcount_singleplanet'] = counts[ipub]['totalcount']
    if verbose:
        print('total counts before', [count['totalcount'] for count in counts])
    for ipub in range(Npublications):
        thisref = counts[ipub]['totalref']
        thisplanet = counts[ipub]['planet']
        for otherplanet in data['planets']:
            if otherplanet != thisplanet:
                # loop through all the other pubs to find this same ref but for other planets
                #  there should usually be a match, but not always
                additionalCounts = 0
                for iloop in range(Npublications):
                    if (
                        counts[iloop]['totalref'] == thisref
                        and counts[iloop]['planet'] == otherplanet
                    ):
                        # if verbose: print('Found other planet info for this ref:',
                        #                  ipub,iloop,thisplanet,otherplanet,thisref)
                        additionalCounts = counts[iloop][
                            'totalcount_singleplanet'
                        ]
                # print('additionalCounts',thisplanet,otherplanet,additionalCounts)
                counts[ipub]['totalcount'] += additionalCounts
    if verbose:
        print('total counts after', [count['totalcount'] for count in counts])
    totalcounts = numpy.array([count['totalcount'] for count in counts])
    bestscore = numpy.max(totalcounts)
    refs = numpy.array([count['totalref'] for count in counts])
    refs_with_max_score = refs[numpy.where(totalcounts == bestscore)]
    if verbose:
        print('refs_with_max_score', refs_with_max_score)
    bestyear = 0
    bestref = 'ZZ'
    for ref in refs_with_max_score:
        try:
            year = int(ref[-4:])
        except ValueError:
            # there are some refs without a year; make them lower priority
            year = 1
        # select the most recently published as tiebreaker
        if year >= bestyear:
            # use author name as tie-breaker for 2 with same year
            if year == bestyear and bestref < ref:
                pass
            else:
                bestyear = year
                bestref = ref
                # ref_with_max_score = ref

    if verbose:
        if bestref_oldmethod == bestref:
            print(' (old,new method give the same result)', bestref)
        else:
            print(
                'MULTI-PLANET SELECTION GIVES NEW BEST PUB:',
                bestref_oldmethod,
                bestref,
            )

    # let's keep track of what parameters are missing
    #  should be a lot with planet mass missing
    #  FEH* sometimes missing is not a problem
    # oof actually this is hard. maybe do it during the bestValue() selection

    bestpubIndices = {'star': -1}
    for planet in data['planets']:
        bestpubIndices[planet] = -1
    for ipub in range(Npublications):
        if data['R*_ref'][ipub] == bestref:
            bestpubIndices['star'] = ipub
    for planet in data['planets']:
        for ipub in range(len(data[planet]['rp_ref'])):
            if data[planet]['rp_ref'][ipub] == bestref:
                bestpubIndices[planet] = ipub
    if verbose:
        if bestpubIndices['star'] == -1:
            print('ERROR!: no matching star pub?!?')
        for planet in data['planets']:
            if bestpubIndices[planet] == -1:
                print('OH NO!: no matching planet pub!', planet)

    # return bestscore,ref_with_max_score,bestpubIndices
    return bestscore, bestref, bestpubIndices


# -------------------------------------------------------------------
def estimate_mass_from_radius(radius_Jup):
    '''
    Use an assumed mass-radius relationship to guess the planet mass
    '''

    sscmks = syscons.ssconstants(cgs=True)

    radius_Earth = radius_Jup * sscmks['Rjup'] / sscmks['Rearth']

    ChenKipping_is_acceptable = False

    if ChenKipping_is_acceptable:
        rocky_slope = 0.279
        gas_slope = 0.589
        R_0 = 1.008  # radius for 1 Earth mass
        transition_point_mass = 2.04  # Earth masses
        transition_point_radius = R_0 * transition_point_mass**rocky_slope
        # transition point is 1.23 Earth radii
        # print('power law transition point (R_Earth)',transition_point_radius)

        if radius_Earth < transition_point_radius:
            mass_Earth = (radius_Earth / R_0) ** (1.0 / rocky_slope)
            # at trans: mass_Earth = (transition_point_radius/R_0)**(1./rocky_slope)
            # which simplifies to: mass_Earth = transition_point_mass  check
        else:
            mass_Earth = transition_point_mass * (
                radius_Earth / transition_point_radius
            ) ** (1.0 / gas_slope)
            # at trans: mass_Earth = transition_point_mass

    else:
        # Traub2011 formula, modified to enforce Earth/Jupiter fit
        r0 = 0.375
        r1 = 11.0
        m1 = 1000.0
        w = 1.5
        if radius_Earth > r1:
            # for radius larger than Jupiter, assume Jupiter mass
            mass_Earth = sscmks['Mjup'] / sscmks['Mearth']

        elif radius_Earth < 1:
            # for planets smaller than Earth, assume Earth density
            mass_Earth = radius_Earth**3

        else:
            # regular formula
            mass_Earth = m1 * 10.0 ** (
                -w * ((r1 - radius_Earth) / (radius_Earth - r0)) ** 0.25
            )

    mass_Jup = mass_Earth * sscmks['Mearth'] / sscmks['Mjup']

    return mass_Jup


# -------------------------------------------------------------------
def derive_RHOstar_from_M_and_R(starInfo):
    '''
    If stellar density is blank, fill it in based on R* and M*
    '''

    # get Msun and Rsun definitions, for calculating stellar density from M*,R*
    sscmks = syscons.ssconstants(cgs=True)

    RHO_derived = []
    RHO_lowerr_derived = []
    RHO_uperr_derived = []
    RHO_ref_derived = []

    for R, Rerr1, Rerr2, M, Merr1, Merr2, RHO, RHOerr1, RHOerr2, RHOref in zip(
        starInfo['R*'],
        starInfo['R*_lowerr'],
        starInfo['R*_uperr'],
        starInfo['M*'],
        starInfo['M*_lowerr'],
        starInfo['M*_uperr'],
        starInfo['RHO*'],
        starInfo['RHO*_lowerr'],
        starInfo['RHO*_uperr'],
        starInfo['RHO*_ref'],
    ):

        # check for blank stellar density
        #  (but only update it if M* and R* are both defined)
        if RHO == '' and R != '' and M != '':
            newRHO = (
                float(M)
                * sscmks['Msun']
                / (4.0 * numpy.pi / 3.0 * (float(R) * sscmks['Rsun']) ** 3)
            )
            # RHO_derived.append(str('%6.4f' %newRHO))
            RHO_derived.append(f'{newRHO:6.4f}')
            RHO_ref_derived.append('derived from M*,R*')

            # also fill in the uncertainty on RHO, based on R,M uncertainties
            if Rerr1 == '' or Merr1 == '':
                RHO_lowerr_derived.append('')
            else:
                newRHOfractionalError1 = -numpy.sqrt(
                    (3.0 * float(Rerr1) / float(R)) ** 2
                    + (float(Merr1) / float(M)) ** 2
                )
                # RHO_lowerr_derived.append(str('%6.4f' %(newRHO * newRHOfractionalError1)))
                RHO_lowerr_derived.append(
                    f'{(newRHO * newRHOfractionalError1):6.4f}'
                )
            if Rerr2 == '' or Merr2 == '':
                RHO_uperr_derived.append('')
            else:
                newRHOfractionalError2 = numpy.sqrt(
                    (3.0 * float(Rerr2) / float(R)) ** 2
                    + (float(Merr2) / float(M)) ** 2
                )
                # RHO_uperr_derived.append(str('%6.4f' %(newRHO * newRHOfractionalError2)))
                RHO_uperr_derived.append(
                    f'{(newRHO * newRHOfractionalError2):6.4f}'
                )
        else:
            RHO_derived.append(RHO)
            RHO_lowerr_derived.append(RHOerr1)
            RHO_uperr_derived.append(RHOerr2)
            RHO_ref_derived.append(RHOref)

    return RHO_derived, RHO_lowerr_derived, RHO_uperr_derived, RHO_ref_derived


# -------------------------------------------------------------------
def derive_SMA_from_P_and_Mstar(starInfo, planet_letter):
    '''
    If semi-major axis is blank, fill it in based on the orbital period and star mass
    '''

    # get G, Msun, AU definitions
    sscmks = syscons.ssconstants(cgs=True)

    sma_derived = []
    sma_lowerr_derived = []
    sma_uperr_derived = []
    sma_ref_derived = []

    for (
        Mold,
        Molderr1,
        Molderr2,
        P,
        Perr1,
        Perr2,
        sma,
        sma_err1,
        sma_err2,
        sma_ref,
    ) in zip(
        starInfo['M*'],
        starInfo['M*_lowerr'],
        starInfo['M*_uperr'],
        starInfo[planet_letter]['period'],
        starInfo[planet_letter]['period_lowerr'],
        starInfo[planet_letter]['period_uperr'],
        starInfo[planet_letter]['sma'],
        starInfo[planet_letter]['sma_lowerr'],
        starInfo[planet_letter]['sma_uperr'],
        starInfo[planet_letter]['sma_ref'],
    ):

        # careful: for multi-planet systems,
        #  these stellar and planetary parameters are misaligned.
        # we have to match stellar and planetary params by their reference
        if sma_ref in starInfo['M*_ref']:
            imass = starInfo['M*_ref'].index(sma_ref)
            M = starInfo['M*'][imass]
            Merr1 = starInfo['M*_lowerr'][imass]
            Merr2 = starInfo['M*_uperr'][imass]
            # NOTE: didn't finish testing this completely before python3.12/pymc/etc upgrade
            # if M != Mold: print('new stellar mass',planet_letter, M,'old stellar mass',Mold)
        else:
            # print('STRANGE: this planet ref doesnt exist in star refs:',sma_ref)
            # log.warning('STRANGE: this planet ref doesnt exist in star refs: %s',sma_ref)
            M = Mold
            Merr1 = Molderr1
            Merr2 = Molderr2

        # check for blank stellar density
        #  (but only update it if M* and P are both defined)
        if sma == '' and M != '' and P != '':
            # M = 1
            # P = 365.25
            # giving these Earth parameters results in 0.9999874 AU
            GM = sscmks['G'] * float(M) * sscmks['Msun']
            newsma = (
                GM * (float(P) * sscmks['day'] / 2.0 / numpy.pi) ** 2
            ) ** (1.0 / 3.0)
            newsma /= sscmks['AU']
            # sma_derived.append(str('%6.4f' %newsma))
            sma_derived.append(f'{newsma:6.4f}')
            sma_ref_derived.append('derived from period,M*')

            # also fill in the uncertainty on sma, based on R,M uncertainties
            if Perr1 == '' or Merr1 == '':
                sma_lowerr_derived.append('')
            else:
                newsma_fractionalError1 = -numpy.sqrt(
                    (2.0 / 3.0 * float(Perr1) / float(P)) ** 2
                    + (1.0 / 3.0 * float(Merr1) / float(M)) ** 2
                )
                # sma_lowerr_derived.append(str('%6.4f' %(newsma * newsma_fractionalError1)))
                sma_lowerr_derived.append(
                    f'{(newsma * newsma_fractionalError1):6.4f}'
                )
            if Perr2 == '' or Merr2 == '':
                sma_uperr_derived.append('')
            else:
                newsma_fractionalError2 = numpy.sqrt(
                    (2.0 / 3.0 * float(Perr2) / float(P)) ** 2
                    + (1.0 / 3.0 * float(Merr2) / float(M)) ** 2
                )
                # sma_uperr_derived.append(str('%6.4f' %(newsma * newsma_fractionalError2)))
                sma_uperr_derived.append(
                    f'{(newsma * newsma_fractionalError2):6.4f}'
                )
        else:
            sma_derived.append(sma)
            sma_lowerr_derived.append(sma_err1)
            sma_uperr_derived.append(sma_err2)
            sma_ref_derived.append(sma_ref)

    return sma_derived, sma_lowerr_derived, sma_uperr_derived, sma_ref_derived


# -------------------------------------------------------------------
def derive_LOGGstar_from_R_and_M(starInfo):
    '''
    If stellar log-g is blank, fill it in based on the star's radius and mass
    '''

    # get Msun and Rsun definitions
    sscmks = syscons.ssconstants(cgs=True)

    LOGG_derived = []
    LOGG_lowerr_derived = []
    LOGG_uperr_derived = []
    LOGG_ref_derived = []

    for (
        R,
        Rerr1,
        Rerr2,
        M,
        Merr1,
        Merr2,
        LOGG,
        LOGGerr1,
        LOGGerr2,
        LOGGref,
    ) in zip(
        starInfo['R*'],
        starInfo['R*_lowerr'],
        starInfo['R*_uperr'],
        starInfo['M*'],
        starInfo['M*_lowerr'],
        starInfo['M*_uperr'],
        starInfo['LOGG*'],
        starInfo['LOGG*_lowerr'],
        starInfo['LOGG*_uperr'],
        starInfo['LOGG*_ref'],
    ):

        # check for blank stellar log-g
        #  (but only update it if M* and R* are both defined)
        if LOGG == '' and R != '' and M != '':
            newLOGG = calculate_logg(M, R, sscmks, units='solar')
            # LOGG_derived.append(str('%6.4f' %newLOGG))
            LOGG_derived.append(f'{newLOGG:6.4f}')
            LOGG_ref_derived.append('derived from M*,R*')

            # also fill in the uncertainty on LOGG, based on R,M uncertainties
            if Rerr1 == '' or Merr1 == '':
                LOGG_lowerr_derived.append('')
            else:
                newLOGGfractionalError1 = -numpy.sqrt(
                    (2.0 * float(Rerr1) / float(R)) ** 2
                    + (float(Merr1) / float(M)) ** 2
                )
                # LOGG_lowerr_derived.append(str('%6.4f' %numpy.log10
                LOGG_lowerr_derived.append(
                    f'{numpy.log10(1 + newLOGGfractionalError1):6.4f}'
                )
            if Rerr2 == '' or Merr2 == '':
                LOGG_uperr_derived.append('')
            else:
                newLOGGfractionalError2 = numpy.sqrt(
                    (2.0 * float(Rerr2) / float(R)) ** 2
                    + (float(Merr2) / float(M)) ** 2
                )
                # LOGG_uperr_derived.append(str('%6.4f' %numpy.log10
                LOGG_uperr_derived.append(
                    f'{numpy.log10(1 + newLOGGfractionalError2):6.4f}'
                )
        else:
            LOGG_derived.append(LOGG)
            LOGG_lowerr_derived.append(LOGGerr1)
            LOGG_uperr_derived.append(LOGGerr2)
            LOGG_ref_derived.append(LOGGref)

    return (
        LOGG_derived,
        LOGG_lowerr_derived,
        LOGG_uperr_derived,
        LOGG_ref_derived,
    )


# -------------------------------------------------------------------
def derive_LOGGplanet_from_R_and_M(starInfo, planet_letter, verbose=False):
    '''
    If planetary log-g is blank, fill it in based on the planet's radius and mass
    '''

    # get MJup and RJup definitions
    sscmks = syscons.ssconstants(cgs=True)

    logg_derived = []
    logg_lowerr_derived = []
    logg_uperr_derived = []
    logg_ref_derived = []

    # print()
    # print(starInfo[planet_letter])
    # print(starInfo[planet_letter].keys())
    # print()

    # this allows for calls from overwrite, where the dictionary is filled with floats, not lists
    # print(starInfo[planet_letter]['rp'], type(starInfo[planet_letter]['rp']))
    if not isinstance(starInfo[planet_letter]['rp'], list):
        starInfo[planet_letter]['rp'] = [starInfo[planet_letter]['rp']]
        starInfo[planet_letter]['rp_lowerr'] = [
            starInfo[planet_letter]['rp_lowerr']
        ]
        starInfo[planet_letter]['rp_uperr'] = [
            starInfo[planet_letter]['rp_uperr']
        ]
        starInfo[planet_letter]['mass'] = [starInfo[planet_letter]['mass']]
        starInfo[planet_letter]['mass_lowerr'] = [
            starInfo[planet_letter]['mass_lowerr']
        ]
        starInfo[planet_letter]['mass_uperr'] = [
            starInfo[planet_letter]['mass_uperr']
        ]
    else:
        # this one is different from other subroutines here - the 'logg' field doesn't exist yet
        if 'logg' in starInfo[planet_letter].keys() and verbose:
            print('ERROR: logg field shouldnt exist yet')

    # for R,Rerr1,Rerr2, M,Merr1,Merr2, logg,loggerr1,loggerr2,loggref in zip(
    for R, Rerr1, Rerr2, M, Merr1, Merr2 in zip(
        starInfo[planet_letter]['rp'],
        starInfo[planet_letter]['rp_lowerr'],
        starInfo[planet_letter]['rp_uperr'],
        starInfo[planet_letter]['mass'],
        starInfo[planet_letter]['mass_lowerr'],
        starInfo[planet_letter]['mass_uperr'],
    ):
        # starInfo[planet_letter]['logg'],
        # starInfo[planet_letter]['logg_lowerr'],
        # starInfo[planet_letter]['logg_uperr'],
        # starInfo[planet_letter][planet_letter]['logg_ref']):

        # if logg=='' and R!='' and M!='':
        # no need to check for blank planetary log-g; it doesn't exist in Archive table
        if R != '' and M != '':
            logg = calculate_logg(M, R, sscmks, units='Jupiter')
            # print('M,R,logg',M,R,logg)
            # logg_derived.append(str('%6.4f' %logg))
            logg_derived.append(f'{logg:6.4f}')

            logg_ref_derived.append('derived from Mp,Rp')

            # also fill in the uncertainty on logg, based on R,M uncertainties
            if Rerr1 == '' or Merr1 == '':
                logg_lowerr_derived.append('')
            else:
                loggfractionalError1 = -numpy.sqrt(
                    (2.0 * float(Rerr1) / float(R)) ** 2
                    + (float(Merr1) / float(M)) ** 2
                )
                # print()
                # print('Rp fractional error',float(Rerr1),float(R))
                # print('Rp fractional error',float(Rerr1)/float(R))
                # print('Mp fractional error',float(Merr1),float(M))
                # print('Mp fractional error',float(Merr1)/float(M))
                # print('logg fractional error',-loggfractionalError1)
                # if loggfractionalError1
                # logg_lowerr_derived.append(str('%6.4f' %numpy.log10

                # this conditional avoids log(negative) error
                if loggfractionalError1 > -1:
                    logg_lowerr_derived.append(
                        f'{numpy.log10(1 + loggfractionalError1):6.4f}'
                    )
                else:
                    # HD 23472 Trifonov2019 has large R,M error bars, probably a typo
                    # doesn't matter, since it's not the default publication
                    logg_lowerr_derived.append('-1')
            if Rerr2 == '' or Merr2 == '':
                logg_uperr_derived.append('')
            else:
                loggfractionalError2 = numpy.sqrt(
                    (2.0 * float(Rerr2) / float(R)) ** 2
                    + (float(Merr2) / float(M)) ** 2
                )
                # logg_uperr_derived.append(str('%6.4f' %numpy.log10
                #                              (1 + loggfractionalError2)))
                logg_uperr_derived.append(
                    f'{numpy.log10(1 + loggfractionalError2):6.4f}'
                )
        else:
            logg_derived.append('')
            logg_lowerr_derived.append('')
            logg_uperr_derived.append('')
            logg_ref_derived.append('')

    return (
        logg_derived,
        logg_lowerr_derived,
        logg_uperr_derived,
        logg_ref_derived,
    )


# -------------------------------------------------------------------
def derive_Lstar_from_R_and_T(starInfo):
    '''
    If stellar luminosity is blank, calculate it from the star's radius and temperature
    '''

    # get Tsun definition
    sscmks = syscons.ssconstants(cgs=True)

    Lstar_derived = []
    Lstar_lowerr_derived = []
    Lstar_uperr_derived = []
    Lstar_ref_derived = []

    # print(starInfo['L*'])
    # print(starInfo['L*_lowerr'])
    # print(starInfo['L*_uperr'])

    for (
        R,
        Rerr1,
        Rerr2,
        T,
        Terr1,
        Terr2,
        Lstar,
        Lstarerr1,
        Lstarerr2,
        Lstarref,
    ) in zip(
        starInfo['R*'],
        starInfo['R*_lowerr'],
        starInfo['R*_uperr'],
        starInfo['T*'],
        starInfo['T*_lowerr'],
        starInfo['T*_uperr'],
        starInfo['L*'],
        starInfo['L*_lowerr'],
        starInfo['L*_uperr'],
        starInfo['L*_ref'],
    ):

        # check for blank stellar luminosity
        #  (but only update it if R* and T* are both defined)
        if Lstar == '' and R != '' and T != '':
            newLstar = float(R) ** 2 * (float(T) / sscmks['Tsun']) ** 4
            # print('Lstar derived',newLstar)

            if newLstar > 1.0e-3:
                Lstar_derived.append(f'{newLstar:6.4f}')
            else:
                # WD 1856 has a very low luminosity (it's a white dwarf); needs more sig-figs
                Lstar_derived.append(f'{newLstar:10.8f}')
            Lstar_ref_derived.append('derived from R*,T*')

            # also fill in the uncertainty on Lstar, based on R,M uncertainties
            if Rerr1 == '' or Terr1 == '':
                Lstar_lowerr_derived.append('')
            else:
                newLstarfractionalError1 = -numpy.sqrt(
                    (2.0 * float(Rerr1) / float(R)) ** 2
                    + (4.0 * float(Terr1) / float(T)) ** 2
                )
                Lstar_lowerr_derived.append(
                    f'{(newLstar * newLstarfractionalError1):6.4f}'
                )
            if Rerr2 == '' or Terr2 == '':
                Lstar_uperr_derived.append('')
            else:
                newLstarfractionalError2 = numpy.sqrt(
                    (2.0 * float(Rerr2) / float(R)) ** 2
                    + (4.0 * float(Terr2) / float(T)) ** 2
                )
                Lstar_uperr_derived.append(
                    f'{(newLstar * newLstarfractionalError2):6.4f}'
                )
        else:
            # SPECIAL ATTENTION FOR L*: the Archive gives it logged; correct here
            if Lstar != '':
                Lstar = 10.0 ** float(Lstar)
                # print('Lstar fixed',Lstar)
                # print('Lstar errors before',Lstarerr1,Lstarerr2)
                # there's an Archive bug where some upper error bars are missing, eg HAT-P-40
                if Lstarerr1 == '' and Lstarerr2 != '':
                    Lstarerr1 = -float(Lstarerr2)
                if Lstarerr1 != '':
                    Lstarerr1 = Lstar * (10 ** float(Lstarerr1) - 1)
                    Lstarerr1 = f'{Lstarerr1:6.4f}'
                if Lstarerr2 != '':
                    Lstarerr2 = Lstar * (10 ** float(Lstarerr2) - 1)
                    Lstarerr2 = f'{Lstarerr2:6.4f}'
                Lstar = f'{Lstar:6.4f}'
                # print('Lstar errors after ',Lstarerr1,Lstarerr2)
            Lstar_derived.append(Lstar)
            Lstar_lowerr_derived.append(Lstarerr1)
            Lstar_uperr_derived.append(Lstarerr2)
            Lstar_ref_derived.append(Lstarref)

    return (
        Lstar_derived,
        Lstar_lowerr_derived,
        Lstar_uperr_derived,
        Lstar_ref_derived,
    )


# -------------------------------------------------------------------
def derive_Teqplanet_from_Lstar_and_sma(starInfo, planet_letter, verbose=False):
    '''
    If planet T_equilibrium is blank, calculate it from star luminosity, planet sma
    '''

    # get Lsun definition.  (not needed if we scale to solar)
    # sscmks = syscons.ssconstants(cgs=True)
    # F_1AU = sscmks['Lsun'] / 4./numpy.pi / sscmks['AU']**2
    # sigrad = 5.6704e-5  # cgs
    # T_1AU = (F_1AU / 4. / sigrad)**0.25
    # print('T_eq for Earth',T_1AU)
    T_1AU = 278.33  # K
    # NOTE: this equilibrium temperature assumes albedo=0

    # this allows for calls from overwrite, where the dictionary is filled with floats, not lists
    # print('sma',starInfo[planet_letter]['sma'], type(starInfo[planet_letter]['sma']))
    if not isinstance(starInfo[planet_letter]['sma'], list):
        starInfo[planet_letter]['sma'] = [starInfo[planet_letter]['sma']]
        starInfo[planet_letter]['sma_lowerr'] = [
            starInfo[planet_letter]['sma_lowerr']
        ]
        starInfo[planet_letter]['sma_uperr'] = [
            starInfo[planet_letter]['sma_uperr']
        ]

    Teq_derived = []
    Teq_lowerr_derived = []
    Teq_uperr_derived = []
    Teq_ref_derived = []
    for (
        Lstar,
        Lstarerr1,
        Lstarerr2,
        sma,
        smaerr1,
        smaerr2,
        Teq,
        Teqerr1,
        Teqerr2,
        Teqref,
    ) in zip(
        starInfo['L*'],
        starInfo['L*_lowerr'],
        starInfo['L*_uperr'],
        starInfo[planet_letter]['sma'],
        starInfo[planet_letter]['sma_lowerr'],
        starInfo[planet_letter]['sma_uperr'],
        starInfo[planet_letter]['teq'],
        starInfo[planet_letter]['teq_lowerr'],
        starInfo[planet_letter]['teq_uperr'],
        starInfo[planet_letter]['teq_ref'],
    ):

        # check for blank equilibrium temperature
        #  (but only update it if L* and sma are both defined)
        # 6/12/24 actually let's always update it, if possible;
        #   published values are inconsistently derived (e.g. albedo)
        # if Teq=='' and Lstar!='' and sma!='':
        if Lstar != '' and sma != '':
            newTeq = T_1AU * float(Lstar) ** 0.25 / float(sma) ** 0.5
            # print('Lstar derived',newTeq)
            if verbose and Teq != '':
                print(
                    'updating a published Teq (new,new/old):',
                    planet_letter,
                    newTeq,
                    newTeq / Teq,
                )

            Teq_derived.append(f'{newTeq:6.4f}')
            Teq_ref_derived.append('derived from L*,sma')

            # also fill in the uncertainty on Teq, based on R,M uncertainties
            if Lstarerr1 == '' or smaerr1 == '':
                Teq_lowerr_derived.append('')
            else:
                newTeqfractionalError1 = -numpy.sqrt(
                    (0.25 * float(Lstarerr1) / float(Lstar)) ** 2
                    + (0.5 * float(smaerr1) / float(sma)) ** 2
                )
                Teq_lowerr_derived.append(
                    f'{(newTeq * newTeqfractionalError1):6.4f}'
                )
            if Lstarerr2 == '' or smaerr2 == '':
                Teq_uperr_derived.append('')
            else:
                newTeqfractionalError2 = numpy.sqrt(
                    (0.25 * float(Lstarerr2) / float(Lstar)) ** 2
                    + (0.5 * float(smaerr2) / float(sma)) ** 2
                )
                Teq_uperr_derived.append(
                    f'{(newTeq * newTeqfractionalError2):6.4f}'
                )
        else:
            Teq_derived.append(Teq)
            Teq_lowerr_derived.append(Teqerr1)
            Teq_uperr_derived.append(Teqerr2)
            Teq_ref_derived.append(Teqref)

    return Teq_derived, Teq_lowerr_derived, Teq_uperr_derived, Teq_ref_derived


# -------------------------------------------------------------------
def derive_inclination_from_impactParam(starInfo, planet_letter):
    '''
    If planet inclination is blank, calculate it from impact param, star radius, semi-major axis
    '''

    # get Rsun definition
    sscmks = syscons.ssconstants(cgs=True)

    inc_derived = []
    inc_lowerr_derived = []
    inc_uperr_derived = []
    inc_ref_derived = []

    for (
        Rstar,
        Rstarerr1,
        Rstarerr2,
        sma,
        smaerr1,
        smaerr2,
        impact,
        impacterr1,
        impacterr2,
        inc,
        incerr1,
        incerr2,
        incref,
    ) in zip(
        starInfo['R*'],
        starInfo['R*_lowerr'],
        starInfo['R*_uperr'],
        starInfo[planet_letter]['sma'],
        starInfo[planet_letter]['sma_lowerr'],
        starInfo[planet_letter]['sma_uperr'],
        starInfo[planet_letter]['impact'],
        starInfo[planet_letter]['impact_lowerr'],
        starInfo[planet_letter]['impact_uperr'],
        starInfo[planet_letter]['inc'],
        starInfo[planet_letter]['inc_lowerr'],
        starInfo[planet_letter]['inc_uperr'],
        starInfo[planet_letter]['inc_ref'],
    ):

        # check for blank inclination
        #  (but only update it if impact, R*, and sma are all defined)
        if inc == '' and Rstar != '' and sma != '' and impact != '':

            cosinc = (
                float(impact) * float(Rstar) * sscmks['Rsun/AU'] / float(sma)
            )

            # ok TOI-2669 is weird (inc=76degrees) because the star is puffy (4.1 RSun)
            # print('R* (RSun)',out['priors']['R*'])
            # print('R* (AU)',out['priors']['R*']*sscmks['Rsun/AU'])
            # print('ap (AU)',out['priors'][p]['sma'])
            newinc = numpy.arccos(cosinc) * 180 / numpy.pi
            # print('inclination derived from impact parameter:',newinc)

            inc_derived.append(f'{newinc:6.4f}')
            inc_ref_derived.append('derived from impact parameter')

            # also fill in the uncertainty on inc, based on impact,R*,sma uncertainties
            if Rstarerr1 == '' or smaerr1 == '' or impacterr1 == '':
                inc_lowerr_derived.append('')
            else:
                cosincfractionalError1 = -numpy.sqrt(
                    (float(Rstarerr1) / float(Rstar)) ** 2
                    + (float(impacterr1) / float(impact)) ** 2
                    + (float(smaerr1) / float(sma)) ** 2
                )
                inc_lowerr_derived.append(
                    f'{(cosincfractionalError1 * 180 / numpy.pi):6.4f}'
                )
            if Rstarerr2 == '' or smaerr2 == '' or impacterr2 == '':
                inc_uperr_derived.append('')
            else:
                cosincfractionalError2 = numpy.sqrt(
                    (float(Rstarerr2) / float(Rstar)) ** 2
                    + (float(impacterr2) / float(impact)) ** 2
                    + (float(smaerr2) / float(sma)) ** 2
                )
                inc_uperr_derived.append(
                    f'{(cosincfractionalError2 * 180 / numpy.pi):6.4f}'
                )
        else:
            inc_derived.append(inc)
            inc_lowerr_derived.append(incerr1)
            inc_uperr_derived.append(incerr2)
            inc_ref_derived.append(incref)

    return inc_derived, inc_lowerr_derived, inc_uperr_derived, inc_ref_derived


# -------------------------------------------------------------------
def derive_impactParam_from_inclination(starInfo, planet_letter):
    '''
    If planet impact parameter is blank, calculate it from inclination, star radius, semi-major axis
    '''

    # get Rsun definition
    sscmks = syscons.ssconstants(cgs=True)

    imp_derived = []
    imp_lowerr_derived = []
    imp_uperr_derived = []
    imp_ref_derived = []

    for (
        Rstar,
        Rstarerr1,
        Rstarerr2,
        sma,
        smaerr1,
        smaerr2,
        inc,
        incerr1,
        incerr2,
        imp,
        imperr1,
        imperr2,
        impref,
    ) in zip(
        starInfo['R*'],
        starInfo['R*_lowerr'],
        starInfo['R*_uperr'],
        starInfo[planet_letter]['sma'],
        starInfo[planet_letter]['sma_lowerr'],
        starInfo[planet_letter]['sma_uperr'],
        starInfo[planet_letter]['inc'],
        starInfo[planet_letter]['inc_lowerr'],
        starInfo[planet_letter]['inc_uperr'],
        starInfo[planet_letter]['impact'],
        starInfo[planet_letter]['impact_lowerr'],
        starInfo[planet_letter]['impact_uperr'],
        starInfo[planet_letter]['impact_ref'],
    ):

        # check for blank impact parameter
        #  (but only update it if inclination, R*, and sma are all defined)
        if imp == '' and Rstar != '' and sma != '' and inc != '':

            cosinc = numpy.cos(float(inc) / 180.0 * numpy.pi)
            newimp = cosinc * float(sma) / float(Rstar) / sscmks['Rsun/AU']
            # print()
            # print(' cosinc',cosinc)
            # print(' sma/R*',float(sma) / float(Rstar)/sscmks['Rsun/AU'])
            # print('impact derived from inclination:',newimp)
            newimp = numpy.abs(newimp)
            # if newimp > 1: print('STRANGE impact parameter',newimp)
            # if newimp > 1: newimp = 1
            newimp = min(newimp, 1)

            imp_derived.append(f'{newimp:6.4f}')
            imp_ref_derived.append('derived from inclination')

            # also fill in the uncertainty on impact param, based on inc,R*,sma uncertainties
            if Rstarerr1 == '' or smaerr1 == '' or incerr1 == '':
                imp_lowerr_derived.append('')
            else:
                impfractionalError1 = -numpy.sqrt(
                    (float(Rstarerr1) / float(Rstar)) ** 2
                    + (float(incerr1) * numpy.pi / 180.0 / cosinc) ** 2
                    + (float(smaerr1) / float(sma)) ** 2
                )
                imp_lowerr_derived.append(
                    f'{(newimp * impfractionalError1):6.4f}'
                )
            if Rstarerr2 == '' or smaerr2 == '' or incerr2 == '':
                imp_uperr_derived.append('')
            else:
                impfractionalError2 = numpy.sqrt(
                    (float(Rstarerr2) / float(Rstar)) ** 2
                    + (float(incerr2) * numpy.pi / 180.0 / cosinc) ** 2
                    + (float(smaerr2) / float(sma)) ** 2
                )
                imp_uperr_derived.append(
                    f'{(newimp * impfractionalError2):6.4f}'
                )
        else:
            imp_derived.append(imp)
            imp_lowerr_derived.append(imperr1)
            imp_uperr_derived.append(imperr2)
            imp_ref_derived.append(impref)

    return imp_derived, imp_lowerr_derived, imp_uperr_derived, imp_ref_derived


# -------------------------------------------------------------------
def derive_sma_from_ars(starInfo, planet_letter):
    '''
    If planet semi-major axis is blank, calculate it a/Rp (e.g. Stassun)
    '''

    # get Rsun definition
    sscmks = syscons.ssconstants(cgs=True)

    sma_derived = []
    sma_lowerr_derived = []
    sma_uperr_derived = []
    sma_ref_derived = []

    for (
        Rstar,
        Rstarerr1,
        Rstarerr2,
        ars,
        arserr1,
        arserr2,
        sma,
        smaerr1,
        smaerr2,
        smaref,
    ) in zip(
        starInfo['R*'],
        starInfo['R*_lowerr'],
        starInfo['R*_uperr'],
        starInfo[planet_letter]['ars'],
        starInfo[planet_letter]['ars_lowerr'],
        starInfo[planet_letter]['ars_uperr'],
        starInfo[planet_letter]['sma'],
        starInfo[planet_letter]['sma_lowerr'],
        starInfo[planet_letter]['sma_uperr'],
        starInfo[planet_letter]['sma_ref'],
    ):

        # check for blank semi-major axis
        #  (but only update it if R* and a/R* are defined)
        if sma == '' and Rstar != '' and ars != '':
            newsma = float(ars) * float(Rstar) * sscmks['Rsun/AU']

            sma_derived.append(f'{newsma:6.4f}')
            sma_ref_derived.append('derived from a/R*')

            # also fill in the uncertainty on sma, based on R*,a/R* uncertainties
            #  hmm, these are going to be larger than they should
            #   since the two uncertainties are correlated
            #  it's probably better to just use the a/R* uncertainty straight up
            # no wait, forget that.  they are actually two independent measurements
            if Rstarerr1 == '' or arserr1 == '':
                sma_lowerr_derived.append('')
            else:
                smafractionalError1 = -numpy.sqrt(
                    (float(Rstarerr1) / float(Rstar)) ** 2
                    + (float(arserr1) / float(ars)) ** 2
                )
                # smafractionalError1 = -numpy.sqrt((float(arserr1)/float(ars))**2)
                sma_lowerr_derived.append(
                    f'{(newsma * smafractionalError1):6.4f}'
                )
            if Rstarerr2 == '' or arserr2 == '':
                sma_uperr_derived.append('')
            else:
                smafractionalError2 = numpy.sqrt(
                    (float(Rstarerr2) / float(Rstar)) ** 2
                    + (float(arserr2) / float(ars)) ** 2
                )
                # smafractionalError2 = numpy.sqrt((float(arserr2)/float(ars))**2)
                sma_uperr_derived.append(
                    f'{(newsma * smafractionalError2):6.4f}'
                )
        else:
            sma_derived.append(sma)
            sma_lowerr_derived.append(smaerr1)
            sma_uperr_derived.append(smaerr2)
            sma_ref_derived.append(smaref)

    return sma_derived, sma_lowerr_derived, sma_uperr_derived, sma_ref_derived


# -------------------------------------------------------------------


def fill_in_some_blank_omegas(starInfo, planet_letter):
    '''
    If omega is blank and eccentricity is zero, omega is undefined (fill with zero)
    '''

    omega_filled = []
    omega_lowerr_filled = []
    omega_uperr_filled = []
    omega_ref_filled = []

    for ecc, omega, omegalowerr, omegauperr, omegaref in zip(
        starInfo[planet_letter]['ecc'],
        starInfo[planet_letter]['omega'],
        starInfo[planet_letter]['omega_lowerr'],
        starInfo[planet_letter]['omega_uperr'],
        starInfo[planet_letter]['omega_ref'],
    ):
        # print('omega,ecc',omega,ecc)
        circular = False
        if ecc != '' and float(ecc) == 0:
            circular = True
        if omega == '' and circular:
            omega_filled.append('0')
            omega_ref_filled.append('derived from ecc=0 (undefined)')
            omega_lowerr_filled.append('100')
            omega_uperr_filled.append('100')
        else:
            omega_filled.append(omega)
            omega_lowerr_filled.append(omegalowerr)
            omega_uperr_filled.append(omegauperr)
            omega_ref_filled.append(omegaref)

    return (
        omega_filled,
        omega_lowerr_filled,
        omega_uperr_filled,
        omega_ref_filled,
    )


# -------------------------------------------------------------------
def fixZeroUncertainties(starInfo, starParam, planetParam):
    '''
    If any uncertainty is zero, remove it.  We don't want chi-squared=infinity
    '''

    for param in starParam:
        # print('param check',param)
        if param != 'spTyp':
            for i in range(len(starInfo[param + '_uperr'])):
                if str(starInfo[param + '_uperr'][i]) and str(
                    starInfo[param + '_lowerr'][i]
                ):
                    if (
                        float(starInfo[param + '_uperr'][i]) == 0
                        or float(starInfo[param + '_lowerr'][i]) == 0
                    ):
                        # print('zero uncertainty fixed',i,param,
                        #       starInfo[param+'_uperr'][i],
                        #       starInfo[param+'_lowerr'][i])
                        starInfo[param + '_uperr'][i] = ''
                        starInfo[param + '_lowerr'][i] = ''

    for param in planetParam:
        if param != 'logg':
            # print()
            # print('param check',param)
            for planet in starInfo['planets']:
                for i in range(len(starInfo[planet][param + '_uperr'])):
                    if str(starInfo[planet][param + '_uperr'][i]) and str(
                        starInfo[planet][param + '_lowerr'][i]
                    ):
                        # print(float(starInfo[planet][param+'_uperr'][i]),
                        #       float(starInfo[planet][param+'_lowerr'][i]))
                        if (
                            float(starInfo[planet][param + '_uperr'][i]) == 0
                            or float(starInfo[planet][param + '_lowerr'][i])
                            == 0
                        ):
                            # print('zero uncertainty fixed',i,planet+':'+param,
                            #       starInfo[planet][param+'_uperr'][i],
                            #       starInfo[planet][param+'_lowerr'][i])
                            starInfo[planet][param + '_uperr'][i] = ''
                            starInfo[planet][param + '_lowerr'][i] = ''

    return starInfo


# -------------------------------------------------------------------
def checkValidData(starInfo, starParam, planetParam):
    '''
    Verify that the needed parameters exist in the incoming target state vector
    '''

    missingParams = []

    for param in starParam:
        if param not in starInfo.keys():
            missingParams.append(param)
        # if param!='spTyp':    # spectral type doesnt have error bars, but actually they are here ok
        if param + '_uperr' not in starInfo.keys():
            missingParams.append(param + '_uperr')
        if param + '_lowerr' not in starInfo.keys():
            missingParams.append(param + '_lowerr')

    for param in planetParam:
        if (
            param != 'logg'
        ):  # planet logg is derived later; not loaded from Exoplanet Archive
            for planet in starInfo['planets']:
                if param not in starInfo[planet].keys():
                    missingParams.append(planet + ':' + param)
                if param + '_uperr' not in starInfo[planet].keys():
                    missingParams.append(planet + ':' + param + '_uperr')
                if param + '_lowerr' not in starInfo[planet].keys():
                    missingParams.append(planet + ':' + param + '_lowerr')

    if missingParams:
        # print()
        # for planet in starInfo['planets']:
        #     print(' target param check',starInfo[planet].keys())
        # print()
        # print('missing parameters in the target state vector:',missingParams)
        return False
    return True


# -------------------------------------------------------------------


def fixPublishedLimits(
    systemInfo, starLimitReplacements, planetLimitReplacements, verbose=False
):
    '''
    Check for and fix any upper/lower limits
    Fix them as prescribed by the input replacement values
    '''

    for param in starLimitReplacements:
        # print()
        # print('star param',param)
        # print('values before',param,systemInfo[param])
        # print('lim before',param,systemInfo[param+'_lim'])
        # print('refs before',param,systemInfo[param+'_ref'])

        for ipl in range(len(systemInfo[param])):
            if (
                systemInfo[param + '_lim'][ipl] == '1'
                or systemInfo[param + '_lim'][ipl] == '-1'
            ):
                if verbose:
                    print(
                        'Fixing limit',
                        param,
                        systemInfo[param][ipl],
                        systemInfo[param + '_lim'][ipl],
                        systemInfo[param + '_ref'][ipl],
                    )
                if starLimitReplacements[param] != 'keep value':
                    systemInfo[param][ipl] = starLimitReplacements[param]
                systemInfo[param + '_lim'][ipl] = '0'

        # print('values after',param,systemInfo[param])
        # print('lim after',param,systemInfo[param+'_lim'])
        # print('refs after',param,systemInfo[param+'_ref'])

    for param in planetLimitReplacements:
        # print()
        # print('param',param)

        for p in systemInfo['planets']:
            # print('values before',p,param,systemInfo[p][param])
            # print('lim before',p,param,systemInfo[p][param+'_lim'])
            # print('refs before',p,param,systemInfo[p][param+'_ref'])

            for ipl in range(len(systemInfo[p][param])):
                if (
                    systemInfo[p][param + '_lim'][ipl] == '1'
                    or systemInfo[p][param + '_lim'][ipl] == '-1'
                ):
                    if verbose:
                        print(
                            'Fixing limit',
                            p,
                            param,
                            systemInfo[p][param][ipl],
                            systemInfo[p][param + '_lim'][ipl],
                            systemInfo[p][param + '_ref'][ipl],
                        )
                    if planetLimitReplacements[param] != 'keep value':
                        systemInfo[p][param][ipl] = planetLimitReplacements[
                            param
                        ]
                    systemInfo[p][param + '_lim'][ipl] = '0'
                    if param == 'ecc':
                        systemInfo[p][param + '_ref'][
                            ipl
                        ] = 'derived from eccentricity upper limit'

            # print('values after',p,param,systemInfo[p][param])
            # print('lim after',p,param,systemInfo[p][param+'_lim'])
            # print('refs after',p,param,systemInfo[p][param+'_ref'])

    # check for eccentricity values that are really upper limits
    #  use less than 2-sigma or 3-sigma as criteria?
    for p in systemInfo['planets']:
        # print('eccentricity check',p,systemInfo[p]['ecc'][ipl],
        #      'uperr',systemInfo[p]['ecc_uperr'][ipl],
        #      'loerr',systemInfo[p]['ecc_lowerr'][ipl])
        for ipl in range(len(systemInfo[p][param])):
            # note that this fix-limits call comes before the fill-uncertainty call
            #  so a blank ecc might be set to zero, but the lowerr is still blank
            if (
                systemInfo[p]['ecc'][ipl] != ''
                and float(systemInfo[p]['ecc'][ipl]) > 0
                and systemInfo[p]['ecc_lowerr'][ipl] != ''
            ):
                # print('eccentricity check',p,systemInfo[p]['ecc'],
                #       systemInfo[p]['ecc_uperr'],systemInfo[p]['ecc_lowerr'])
                ecc_SNR = -float(systemInfo[p]['ecc'][ipl]) / float(
                    systemInfo[p]['ecc_lowerr'][ipl]
                )
                if ecc_SNR < 2:
                    if verbose:
                        print(
                            'SYSTEM: should be an ecc upper limit!? SNR=%s',
                            ecc_SNR,
                        )
                    systemInfo[p]['ecc'][ipl] = '0'
                    systemInfo[p]['ecc_ref'][
                        ipl
                    ] = 'derived from eccentricity upper limit'

    # print()
    return systemInfo


# -------------------------------------------------------------------
