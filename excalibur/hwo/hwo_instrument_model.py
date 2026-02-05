'''hwo hwo ds'''

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

import astropy.units as u

log = logging.getLogger(__name__)


# ---------------------------- ---------------------------------------


def load_hwo_instrument(target, runtime_params):
    '''
    Load in the HWO instrument response - uncertainty as a function of wavelength

    Uncertainty is for a single visit;
    number of observed transits is taken into account later
    '''
    thorngren = runtime_params.thorngrenMassMetals
    chachan = runtime_params.chachanMassMetals

    noise_model_dir = excalibur.context['data_dir'] + '/hwo/'

    noise_model_filenames = [
        'arielRad_29oct2024_mmwThorngren-part1.h5',
        'arielRad_29oct2024_mmwThorngren-part2.h5',
    ]
    # print('FILENAMES:',noise_model_filenames)

    hwo_instrument = None
    for noise_model_filename in noise_model_filenames:
        # if the file is missing, note the error and move on
        if not os.path.isfile(noise_model_dir + noise_model_filename):
            log.warning(
                '--< HWOSIM ERROR: SNR file missing for %s  %s >--',
                target,
                noise_model_filename,
            )
        else:
            with h5py.File(
                noise_model_dir + noise_model_filename, 'r'
            ) as hwomodel:

                targets = hwomodel['errorbars_evaluated'].keys()

                if target not in targets:
                    print(
                        'NOTE: target not in HWO SNR file; failing to simulate spectrum',
                        target,
                    )
                    log.warning(
                        '--< HWOSIM: target not in SNR file; failing to simulate spectrum  %s >--',
                        target,
                    )
                    hwo_instrument = None
                    print('  not in this part:', noise_model_filename)
                    pass
                else:
                    # print('  found it in this part:',noise_model_filename)

                    # use 'SNR' table to determine the required number of transits
                    # SNR is not an hdf5 table.  just access it like a normal dict
                    SNRtable = hwomodel['SNR']['SNRTab_to_group']
                    # print('SNR options',SNRtable.keys())

                    # assume tier-2
                    nTransitsCH0 = SNRtable['AIRS-CH0-T2-nTransits']['value'][
                        ()
                    ]
                    nTransitsCH1 = SNRtable['AIRS-CH1-T2-nTransits']['value'][
                        ()
                    ]
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
                            '--< HWO #-of-visits missing: %s >--', target
                        )
                    elif len(thisplanetIndex) > 1:
                        log.warning(
                            '--< HWO has multiple target matches?!: %s >--',
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
                            '--< HWO has non-finite # of visits: %s >--',
                            target,
                        )
                    log.info(
                        '--< HWO requires %s visits for %s >--',
                        str(nVisits),
                        target,
                    )

                    noiseSpectrum = read_table_hdf5(
                        hwomodel['errorbars_evaluated'][target],
                        path='table',
                    )
                    # print('noiseSpectrum options',noiseSpectrum.keys())
                    hwo_instrument = {
                        'nVisits': nVisits,
                        'transitDuration': transitDuration,
                        'wavelength': noiseSpectrum['Wavelength'].value,
                        'wavelow': noiseSpectrum['LeftBinEdge'].value,
                        'wavehigh': noiseSpectrum['RightBinEdge'].value,
                        'noise': noiseSpectrum['NoiseOnTransitFloor'].value,
                    }

                    for iwave in range(len(hwo_instrument['wavelow']) - 1):
                        # multiply by number slightly above 1 to deal with numerical precision error
                        if (
                            hwo_instrument['wavehigh'][iwave]
                            > hwo_instrument['wavelow'][iwave + 1] * 1.00001
                        ):
                            log.info(
                                '--< HWOSIM adjusting wavelength grid: %s wave=%s >--',
                                target,
                                hwo_instrument['wavelength'][iwave],
                            )
                            hwo_instrument['wavehigh'][iwave] = (
                                hwo_instrument['wavelow'][iwave + 1] * 0.99999
                            )

    if not hwo_instrument:
        log.warning(
            '--< HWOSIM: target not in SNR files; failing to simulate spectrum  for %s >--',
            target,
        )

    return hwo_instrument
