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
import astropy.constants as const
import astropy.units as u

log = logging.getLogger(__name__)


# ---------------------------- ---------------------------------------

def load_hwo_instrument(target, system_params):
    '''
    Load in the HWO instrument response - uncertainty as a function of wavelength

    Uncertainty is for a single visit;
    number of observed transits is taken into account later
    '''

    def reso_range(start, finish, res):
        '''
        Builds a wavelength grid with specified endpoints and resolution
        '''
        wl_low = [start]
        res = 1.0 / res
        wl_high = [start + (start * res)]
        while wl_high[-1] < finish:
            wl_low.append(wl_high[-1])
            wl_high.append(wl_low[-1] + (wl_low[-1] * res))
        bins = np.array([wl_low, wl_high]).T

        return bins

    Tstar = system_params['T*'] * u.K
    Rstar = system_params['R*'] * const.R_sun
    d = system_params['dist'] * u.pc
    print('Rstar', Rstar)
    print('d', d)
    
    planet_letter = target[-1]
    # print('available planet parameters:', system_params[planet_letter].keys())
    transit_duration = system_params[planet_letter]['trandur'] * u.hr
    print('transit duration (hours)', transit_duration)

    # assumed instrument parameters
    tel_diameter = 7.2 * u.m  # circumscribed diameter of EAC1
    obsc_fact = 0.9  # correction factor for tel area, due to obscuration (if on-axis) and segemented mirrors
    tel_area = (
        np.pi * (tel_diameter / 2) ** 2 * obsc_fact
    )  # telescope collecting area

    wavelow = 0.2  # minimum wavelength
    wavehigh = 1.8  # maximum wavelength
    waveres = 500  # wavelength resolution
    wavelength_bins = reso_range(wavelow, wavehigh, waveres)  # wavelength grid
    wavelength = 0.5 * (wavelength_bins[:, 0] + wavelength_bins[:, 1])
    bin_widths = wavelength_bins[:, 1] - wavelength_bins[:, 0]

    wavelength = wavelength * u.um
    bin_widths = bin_widths * u.um

    optical_thruput = 0.626  # from 2 XeLiF and 17 Ag surfaces measured at 1 micron; based on EAC1 (no 0.5 from polarizer/dichroic, no IFS loss)
    
    # star properties
    BB = (
        2
        * const.h
        * const.c**2
        / wavelength**5
        / (np.exp(const.h * const.c / (wavelength * const.k_B * Tstar)) - 1.0)
    )  # flux per wl per solid angle (blackbody spectral radiance)
    print('BB', BB.decompose())
    F_atEarth = BB * np.pi * (Rstar / d) ** 2  # flux received at Earth
    print('flux', F_atEarth.decompose())

    # transit
    transit_in = transit_duration
    transit_out = transit_in * 1.0

    # collect photons
    tel_power = F_atEarth * bin_widths * tel_area * optical_thruput
    photon_energy = const.h * const.c / wavelength
    photon_rate = tel_power / photon_energy
    # print ('photon_rate', photon_rate)
    photon_rate = photon_rate.decompose()
    print ('photon_rate', photon_rate.decompose())
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
    shotnoise = shotnoise.decompose()
    print('fractional uncertainty on the stellar flux',
          np.median(shotnoise), np.std(shotnoise))
    
    nVisits = 1  # number of visits

    hwo_instrument = {
        'nVisits': nVisits,
        'transitDuration': transit_in,
        'wavelength': wavelength / u.um,
        'wavelow': wavelength_bins[:, 0],
        'wavehigh': wavelength_bins[:, 1],
        'noise': shotnoise,
    }

    return hwo_instrument
