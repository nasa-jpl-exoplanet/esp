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
import scipy.constants as const

log = logging.getLogger(__name__)


# ---------------------------- ---------------------------------------


# def load_hwo_instrument(target):
#     '''
#     Load in the HWO instrument response - uncertainty as a function of wavelength

#     Uncertainty is for a single visit;
#     number of observed transits is taken into account later
#     '''

#     noise_model_dir = excalibur.context['data_dir'] + '/hwo/'

#     noise_model_filenames = [
#         'arielRad_29oct2024_mmwThorngren-part1.h5',
#         'arielRad_29oct2024_mmwThorngren-part2.h5',
#     ]
#     # print('FILENAMES:',noise_model_filenames)

#     hwo_instrument = None
#     for noise_model_filename in noise_model_filenames:
#         # if the file is missing, note the error and move on
#         if not os.path.isfile(noise_model_dir + noise_model_filename):
#             log.warning(
#                 '--< HWOSIM ERROR: SNR file missing for %s  %s >--',
#                 target,
#                 noise_model_filename,
#             )
#         else:
#             with h5py.File(
#                 noise_model_dir + noise_model_filename, 'r'
#             ) as hwomodel:

#                 targets = hwomodel['errorbars_evaluated'].keys()

#                 if target not in targets:
#                     pass
#                 else:
#                     # print('  found it in this part:',noise_model_filename)

#                     # use 'SNR' table to determine the required number of transits
#                     # SNR is not an hdf5 table.  just access it like a normal dict
#                     SNRtable = hwomodel['SNR']['SNRTab_to_group']
#                     # print('SNR options',SNRtable.keys())

#                     # assume tier-2
#                     nTransitsCH0 = SNRtable['AIRS-CH0-T2-nTransits']['value'][
#                         ()
#                     ]
#                     nTransitsCH1 = SNRtable['AIRS-CH1-T2-nTransits']['value'][
#                         ()
#                     ]
#                     T14 = SNRtable['T14']['value'][()]
#                     planetNames = SNRtable['planetName']['value'][()]
#                     planetNames = np.array(
#                         [name.decode('UTF-8') for name in planetNames]
#                     )

#                     thisplanetIndex = np.where(target == planetNames)
#                     nVisits = 666
#                     transitDuration = 6666
#                     if len(thisplanetIndex) == 0:
#                         log.warning(
#                             '--< HWO #-of-visits missing: %s >--', target
#                         )
#                     elif len(thisplanetIndex) > 1:
#                         log.warning(
#                             '--< HWO has multiple target matches?!: %s >--',
#                             target,
#                         )
#                     elif np.isfinite(
#                         nTransitsCH0[thisplanetIndex]
#                     ) and np.isfinite(nTransitsCH1[thisplanetIndex]):
#                         nVisits = np.min(
#                             [
#                                 nTransitsCH0[thisplanetIndex],
#                                 nTransitsCH1[thisplanetIndex],
#                             ]
#                         )
#                         nVisits = int(np.ceil(nVisits))
#                         transitDuration = T14[thisplanetIndex][0]
#                     else:
#                         log.warning(
#                             '--< HWO has non-finite # of visits: %s >--',
#                             target,
#                         )
#                     log.info(
#                         '--< HWO requires %s visits for %s >--',
#                         str(nVisits),
#                         target,
#                     )

#                     noiseSpectrum = read_table_hdf5(
#                         hwomodel['errorbars_evaluated'][target],
#                         path='table',
#                     )
#                     # print('noiseSpectrum options',noiseSpectrum.keys())
#                     hwo_instrument = {
#                         'nVisits': nVisits,
#                         'transitDuration': transitDuration,
#                         'wavelength': noiseSpectrum['Wavelength'].value,
#                         'wavelow': noiseSpectrum['LeftBinEdge'].value,
#                         'wavehigh': noiseSpectrum['RightBinEdge'].value,
#                         'noise': noiseSpectrum['NoiseOnTransitFloor'].value,
#                     }

#                     for iwave in range(len(hwo_instrument['wavelow']) - 1):
#                         # multiply by number slightly above 1 to deal with numerical precision error
#                         if (
#                             hwo_instrument['wavehigh'][iwave]
#                             > hwo_instrument['wavelow'][iwave + 1] * 1.00001
#                         ):
#                             log.info(
#                                 '--< HWOSIM adjusting wavelength grid: %s wave=%s >--',
#                                 target,
#                                 hwo_instrument['wavelength'][iwave],
#                             )
#                             hwo_instrument['wavehigh'][iwave] = (
#                                 hwo_instrument['wavelow'][iwave + 1] * 0.99999
#                             )

#     if not hwo_instrument:
#         log.warning(
#             '--< HWOSIM: target not in SNR files; failing to simulate spectrum  for %s >--',
#             target,
#         )

#     # wavelength grid parameters
#     Nwave = 300
#     # wavemin = 0.4
#     # wavemax = 1.0
#     wavemin = 0.2
#     wavemax = 1.8

#     wavegridedges = np.linspace(wavemin, wavemax, Nwave + 1)
#     # dwave = wavegridedges[1] - wavegridedges[0]
#     # print('dwave', dwave)

#     wavelow = wavegridedges[:-1]
#     wavehigh = wavegridedges[1:]
#     # waves = wavegridedges[:-1] + dwave/2  # same thing as below
#     waves = (wavelow + wavehigh) / 2
#     # print('wave range', waves[0], waves[10], waves[-1])

#     # instrument noise parameters
#     noiselevel_ppm = 100.0
#     noise = np.array([noiselevel_ppm / 1.0e6] * Nwave)

#     # transitDuration isnt used; this is just for future checking
#     #  (the noise is calculated based on this transit duration)
#     #
#     # nVisits is used; the noise will be multiplied by sqrt(N)

#     hwo_instrument = {
#         'nVisits': 1,
#         'transitDuration': 3600,
#         'wavelength': waves,
#         'wavelow': wavelow,
#         'wavehigh': wavehigh,
#         'noise': noise,
#     }

#     return hwo_instrument


def load_hwo_instrument(target, system_params):
    '''
    Load in the HWO instrument response - uncertainty as a function of wavelength

    Uncertainty is for a single visit;
    number of observed transits is taken into account later
    '''

    def reso_range(start, finish, res):
        '''
        Builds a wavelength grid with specified endpoints and resolution (returns the midpoints of each bin)
        '''
        wl_low = [start]
        res = 1.0 / res
        wl_high = [start + (start * res)]
        while wl_high[-1] < finish:
            wl_low.append(wl_high[-1])
            wl_high.append(wl_low[-1] + (wl_low[-1] * res))

        bins = np.array([wl_low, wl_high]).T

        return bins

    print('available system parameters:', system_params.keys())
    Tstar = system_params['T*']
    Rstar = system_params['R*']
    d = system_params['dist']

    planet_letter = target[-1]
    print('available planet parameters:', system_params[planet_letter].keys())
    transit_duration = system_params[planet_letter]['trandur']
    print('transit duration (hours)', transit_duration)

    # assumed instrument parameters
    tel_diameter = 7.2  # circumscribed diameter of EAC1
    obsc_fact = 0.9  # correction factor for tel area, due to obscuration (if on-axis) and segemented mirrors
    tel_area = (
        np.pi * (tel_diameter / 2) ** 2 * obsc_fact
    )  # telescope collecting area

    wavelow = 0.2 * 1e-6  # minimum wavelength
    wavehigh = 1.8 * 1e-6  # maximum wavelength
    waveres = 500  # wavelength resolution
    wavelength_bins = reso_range(wavelow, wavehigh, waveres)  # wavelength grid
    wavelength = 0.5 * (wavelength_bins[:, 0] + wavelength_bins[:, 1])
    bin_widths = wavelength_bins[:, 1] - wavelength_bins[:, 0]

    optical_thruput = 0.626  # from 2 XeLiF and 17 Ag surfaces measured at 1 micron; based on EAC1 (no 0.5 from polarizer/dichroic, no IFS loss)

    # star properties
    BB = (
        2
        * const.h
        * const.c**2
        / wavelength**5
        / (np.exp(const.h * const.c / (wavelength * const.k * Tstar)) - 1.0)
    )  # flux per wl per solid angle (blackbody spectral radiance)
    F_atEarth = BB * np.pi * (Rstar / d) ** 2  # flux received at Earth

    # transit
    transit_in = transit_duration * 3600
    transit_out = transit_in * 1.0

    # collect photons
    tel_power = F_atEarth * bin_widths * tel_area * optical_thruput
    photon_energy = const.h * const.c / wavelength
    photon_rate = tel_power / photon_energy
    photons_in_transit = (
        photon_rate * transit_in
    )  # total photons in each bin during the transit
    photons_out_transit = (
        photon_rate * transit_out
    )  # total photons in each bin out of transit

    # noise
    shotnoise = np.sqrt(
        1 / photons_in_transit + 1 / photons_out_transit
    )  # shot noise per bin

    nVisits = 1  # number of visits

    hwo_instrument = {
        'nVisits': nVisits,
        'transitDuration': transit_in,
        'wavelength': wavelength * 1e6,
        'wavelow': wavelength_bins[:, 0] * 1e6,
        'wavehigh': wavelength_bins[:, 1] * 1e6,
        'noise': shotnoise,
    }

    return hwo_instrument
