'''ariel ariel ds'''

# Heritage code shame:
# pylint: disable=invalid-name,no-member
# pylint: disable=too-many-branches,too-many-locals,too-many-statements
#
# no-member is for astropy.units

# -- IMPORTS -- ------------------------------------------------------
import os
import numpy as np
import logging

import excalibur

import h5py
from astropy.io.misc.hdf5 import read_table_hdf5

# asdf
# from arielrad.api.run_target import run_target
import astropy.units as u

log = logging.getLogger(__name__)


# ---------------------------- ---------------------------------------
def calculate_ariel_instrument(
    target,
    system_params,
    ancil_params,
    runtime_params,
    verbose=False,
):
    '''
    Use ArielRad to calculate instrument uncertainty as a function of wavelength

    Uncertainty is for a single visit;
    number of observed transits is taken into account later
    '''
    tier = runtime_params.tier
    thorngren = runtime_params.thorngrenMassMetals
    chachan = runtime_params.chachanMassMetals

    planet_letter = target[-1]

    # check for missing data; these should be filled in via system/overwriter
    if system_params[planet_letter]['trandur'] == '':
        log.error(
            'ArielRad Input Error: MISSING TRANSIT DURATION!!! %s', target
        )
        system_params[planet_letter]['trandur'] = 2.0
    if system_params['dist'] == '':
        log.error('ArielRad Input Error: MISSING DISTANCE!!! %s', target)
        system_params['dist'] = '666'
    if system_params[planet_letter]['impact'] == '':
        log.error('ArielRad Input Error: MISSING IMPACT PARAM!!! %s', target)
        system_params[planet_letter]['impact'] = '0'
    if system_params['M*'] == '':
        log.error('ArielRad Input Error: MISSING STAR MASS!!! %s', target)
        system_params['M*'] = '1'

    # SYSTEM PARAMS NEEDED FOR ARIELRAD:
    arielrad_params = {
        'star name': target[:-2],
        'star m': system_params['M*'] * u.M_sun,
        'star r': system_params['R*'] * u.R_sun,
        'star t': system_params['T*'] * u.K,
        'star z': 0.02,
        'star d': system_params['dist'] * u.pc,
        'planet name': target,
        'planet r': system_params[planet_letter]['rp'] * u.R_jup,
        'planet m': system_params[planet_letter]['mass'] * u.M_jup,
        'planet teff': float(system_params[planet_letter]['teq']) * u.K,
        'planet P': system_params[planet_letter]['period'] * u.d,
        'planet a': system_params[planet_letter]['sma'] * u.AU,
        'planet albedo': 0.1,
        'planet T14': system_params[planet_letter]['trandur'] * u.h,
    }
    if chachan:
        arielrad_params['planet mmw'] = (
            ancil_params[planet_letter]['mmw_chachan'] * u.g / u.mol
        )
    elif thorngren:
        arielrad_params['planet mmw'] = (
            ancil_params[planet_letter]['mmw_thorngren'] * u.g / u.mol
        )
    else:
        arielrad_params['planet mmw'] = (
            ancil_params[planet_letter]['mmw_min'] * u.g / u.mol
        )
    # print('SYSTEM PARAMS FOR ARIELRAD:', arielrad_params)

    # Run the simulation
    noise_table, output_dict, nobs, _ = [{}, {}, {}, {}]
    noise_table, output_dict, nobs, info = run_target(
        target_dict=arielrad_params,
        run_config='/proj/sdp/data/arielrad/runConfig.xml',
        obs_mode='transit',
        verbose=verbose,
    )
    if verbose and noise_table == {}:
        print('looks like arielrad calculation was not run')

    # print('output_dict',output_dict.keys())
    # ['tier1', 'tier2', 'tier3', 'SNRTab']
    # print('noise_table',noise_table.keys())
    # ['ch_name', 'wavelength', 'left_bin_edge', 'right_bin_edge', 'total_noise', 'noise_on_transit_floor']
    # print('nobs',nobs)
    # print('version info:', info)

    if tier == 1:
        nVisits = nobs['tier1']
    elif tier == 3:
        nVisits = nobs['tier3']  # divide this by 5
    else:
        if tier != 2:
            log.warning('--< Unknown Ariel Tier!: %s >--', tier)
        nVisits = nobs['tier2']

    ariel_instrument = {
        'nVisits': nVisits,
        'transitDuration': output_dict['SNRTab']['T14'].value,
        'wavelength': noise_table['wavelength'].value,
        'wavelow': noise_table['left_bin_edge'].value,
        'wavehigh': noise_table['right_bin_edge'].value,
        'noise': noise_table['noise_on_transit_floor'].value,
    }
    # print('NUMBER OF VISITS', nVisits)
    # print(
    #    'noisespectrum median,stdev',
    #    np.median(ariel_instrument['noise']),
    #    np.std(ariel_instrument['noise']),
    # )

    for iwave in range(len(ariel_instrument['wavelow']) - 1):
        # multiply by number slightly above 1 to deal with numerical precision error
        if (
            ariel_instrument['wavehigh'][iwave]
            > ariel_instrument['wavelow'][iwave + 1] * 1.00001
        ):
            # print(
            #    'spectral channels overlap!!',
            #    iwave,
            #    ariel_instrument['wavehigh'][iwave],
            #    ariel_instrument['wavelow'][iwave + 1],
            # )
            log.info(
                '--< ARIELSIM adjusting wavelength grid: %s wave=%s >--',
                target,
                ariel_instrument['wavelength'][iwave],
            )
            ariel_instrument['wavehigh'][iwave] = (
                ariel_instrument['wavelow'][iwave + 1] * 0.99999
            )

    if not ariel_instrument:
        log.warning(
            '--< ARIELSIM: target not in SNR files; failing to simulate spectrum  for %s >--',
            target,
        )

    return ariel_instrument


# ---------------------------- ---------------------------------------
def load_ariel_instrument(target, runtime_params):
    '''
    Load in the output from ArielRad - uncertainty as a function of wavelength

    Uncertainty is for a single visit;
    number of observed transits is taken into account later
    '''
    tier = runtime_params.tier
    arielRad_version = runtime_params.arielRad
    thorngren = runtime_params.thorngrenMassMetals
    chachan = runtime_params.chachanMassMetals

    noise_model_dir = excalibur.context['data_dir'] + '/ariel/'

    # noise_model_filenames = ['arielRad_02aug2023.h5']
    # noise_model_filenames = ['arielRad_07aug2023_mmwFixed.h5']
    # noise_model_filenames = ['arielRad_07aug2023_mmwExcal.h5']
    # noise_model_filenames = ['arielRad_3mar2024_MCStargetsAdded.h5']
    # noise_model_filenames = ['arielRad_14aug2024_mmwFixed.h5']
    # 'mmwFixed' means mmw=2.3
    # 'mmwThorngren' means mmw comes from our mass-metallicity relation
    # noise_model_filenames = ['arielRad_14aug2024_mmwThorngren.h5']

    if chachan:
        mmwModel = '_mmwChachan'
    elif thorngren:
        mmwModel = '_mmwThorngren'
    else:
        mmwModel = '_mmwFixed'

    # uh oh, arielrad is crashing, with too many files open
    #  have to split the file into two parts
    noise_model_filenames = [
        'arielRad_' + arielRad_version + mmwModel + '-part1.h5',
        'arielRad_' + arielRad_version + mmwModel + '-part2.h5',
    ]
    # print('FILENAMES:',noise_model_filenames)

    ariel_instrument = None
    for noise_model_filename in noise_model_filenames:
        # if the file is missing, note the error and move on
        if not os.path.isfile(noise_model_dir + noise_model_filename):
            log.warning(
                '--< ARIELSIM ERROR: SNR file missing for %s  %s >--',
                target,
                noise_model_filename,
            )
        else:
            with h5py.File(
                noise_model_dir + noise_model_filename, 'r'
            ) as arielRad_results:
                # print('ArielRad keys',arielRad_results.keys())

                targets = arielRad_results['errorbars_evaluated'].keys()

                if target not in targets:
                    # print('NOTE: target not in Ariel SNR file; failing to simulate spectrum',target)
                    # log.warning('--< ARIELSIM: target not in SNR file; failing to simulate spectrum  %s >--',target)
                    #  ariel_instrument = None
                    # print('  not in this part:',noise_model_filename)
                    pass
                else:
                    # print('  found it in this part:',noise_model_filename)

                    # use 'SNR' table to determine the required number of transits
                    # SNR is not an hdf5 table.  just access it like a normal dict
                    SNRtable = arielRad_results['SNR']['SNRTab_to_group']
                    # print('SNR options',SNRtable.keys())

                    if tier == 1:
                        nTransitsCH0 = SNRtable['AIRS-CH0-T1-nTransits'][
                            'value'
                        ][()]
                        nTransitsCH1 = SNRtable['AIRS-CH1-T1-nTransits'][
                            'value'
                        ][()]
                    elif tier == 3:
                        nTransitsCH0 = SNRtable['AIRS-CH0-T3-nTransits'][
                            'value'
                        ][()]
                        nTransitsCH1 = SNRtable['AIRS-CH1-T3-nTransits'][
                            'value'
                        ][()]
                    else:
                        if tier != 2:
                            log.warning('--< Unknown Ariel Tier!: %s >--', tier)
                        nTransitsCH0 = SNRtable['AIRS-CH0-T2-nTransits'][
                            'value'
                        ][()]
                        nTransitsCH1 = SNRtable['AIRS-CH1-T2-nTransits'][
                            'value'
                        ][()]
                    T14 = SNRtable['T14']['value'][()]
                    planetNames = SNRtable['planetName']['value'][()]
                    planetNames = np.array(
                        [name.decode('UTF-8') for name in planetNames]
                    )

                    thisplanetIndex = np.where(target == planetNames)
                    nVisits = 666
                    transitDuration = 6666
                    if len(thisplanetIndex) == 0:
                        log.warning(
                            '--< ArielRad #-of-visits missing: %s >--', target
                        )
                    elif len(thisplanetIndex) > 1:
                        log.warning(
                            '--< ArielRad has multiple target matches?!: %s >--',
                            target,
                        )
                    elif np.isfinite(
                        nTransitsCH0[thisplanetIndex]
                    ) and np.isfinite(nTransitsCH1[thisplanetIndex]):
                        nVisits = np.min(
                            [
                                nTransitsCH0[thisplanetIndex],
                                nTransitsCH1[thisplanetIndex],
                            ]
                        )
                        nVisits = int(np.ceil(nVisits))
                        transitDuration = T14[thisplanetIndex][0]
                    else:
                        log.warning(
                            '--< ArielRad has non-finite # of visits: %s >--',
                            target,
                        )
                    log.info(
                        '--< ArielRad/Tier-%s requires %s visits for %s >--',
                        str(tier),
                        str(nVisits),
                        target,
                    )

                    noiseSpectrum = read_table_hdf5(
                        arielRad_results['errorbars_evaluated'][target],
                        path='table',
                    )
                    # print('noiseSpectrum options',noiseSpectrum.keys())
                    ariel_instrument = {
                        'nVisits': nVisits,
                        'transitDuration': transitDuration,
                        'wavelength': noiseSpectrum['Wavelength'].value,
                        'wavelow': noiseSpectrum['LeftBinEdge'].value,
                        'wavehigh': noiseSpectrum['RightBinEdge'].value,
                        # 'wavebinsize':(noiseSpectrum['RightBinEdge'].value -
                        #                noiseSpectrum['LeftBinEdge'].value),
                        'noise': noiseSpectrum['NoiseOnTransitFloor'].value,
                    }
                    # print('NUMBER OF VISITS', nVisits)
                    # print('noisespectrum median,stdev',
                    #      np.median(noiseSpectrum['NoiseOnTransitFloor']),
                    #      np.std(noiseSpectrum['NoiseOnTransitFloor']))

                    for iwave in range(len(ariel_instrument['wavelow']) - 1):
                        # multiply by number slightly above 1 to deal with numerical precision error
                        if (
                            ariel_instrument['wavehigh'][iwave]
                            > ariel_instrument['wavelow'][iwave + 1] * 1.00001
                        ):
                            # print('spectral channels overlap!!',iwave,
                            #       ariel_instrument['wavehigh'][iwave],
                            #       ariel_instrument['wavelow'][iwave+1])
                            # log.info('--< ARIELSIM adjusting wavelength grid: %s wave=%s >--',
                            #             target,ariel_instrument['wavelength'][iwave])
                            ariel_instrument['wavehigh'][iwave] = (
                                ariel_instrument['wavelow'][iwave + 1] * 0.99999
                            )

    if not ariel_instrument:
        log.warning(
            '--< ARIELSIM: target not in SNR files; failing to simulate spectrum  for %s >--',
            target,
        )

    return ariel_instrument
