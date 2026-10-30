'''data core ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,many-positional-arguments,too-many-locals

# -- IMPORTS --------------------------------------------------------
import numpy as np
from astropy.io import fits


def read_jwst_noise_parameters(calibrated_paths):
    """
    Read per-detector times and scalar NIRSpec gain estimates. This is from STSci
    """
    average_gain = {"NRS1": 0.996, "NRS2": 1.137}
    integration_time = {}
    gains = {}

    for path in calibrated_paths:
        header = fits.getheader(path, 0)
        detector_name = header["DETECTOR"]
        detector_time = header["EFFINTTM"]
        detector_gain = average_gain[detector_name] * header["GAINFACT"]

        if detector_name in integration_time:
            if (
                integration_time[detector_name] != detector_time
                or gains[detector_name] != detector_gain
            ):
                raise ValueError(
                    f'Integration time or gain varies between files for {detector_name}.'
                )

        integration_time[detector_name] = detector_time
        gains[detector_name] = detector_gain

    return integration_time, gains


def extract_jwst_photon_noise(source_rate, integration_time, gain):
    """
    Estimates the photon noise from JWST observations
    """
    electron_count = source_rate * integration_time * gain
    electron_count = np.where(
        electron_count > 0, electron_count, np.nan
    )  # change negative values to nan
    photon_noise_abs = np.sqrt(electron_count)
    photon_noise_rel = 1.0 / np.sqrt(electron_count)

    output_dict = {
        "electron_count": electron_count,
        "photon_noise_abs": photon_noise_abs,
        "photon_noise_rel": photon_noise_rel,
    }

    return output_dict


def calculate_jwst_native_photon_noise(
    source_rate,
    detector,
    valid,
    integration_time,
    gains,
):
    """
    Calculates JWST native photon noise
    """
    photon_noise_by_detector = {}

    for detector_name in np.unique(detector):
        selection = detector == detector_name

        detector_rate = np.where(
            valid[selection],
            source_rate[selection],
            np.nan,
        )

        if isinstance(integration_time, dict):
            detector_time = integration_time[detector_name]
        else:
            detector_time = integration_time

        photon_noise_by_detector[detector_name] = extract_jwst_photon_noise(
            detector_rate,
            detector_time,
            gains[detector_name],
        )

    return photon_noise_by_detector


def bin_jwst_photon_noise(
    electron_count,
    wavelength,
    in_transit,
    out_of_transit,
    wavelength_low,
    wavelength_high,
):
    """
    Bin electron counts and calculate photon limited transit depth errors
    """
    binned_noise = []
    for lower, upper in zip(wavelength_low, wavelength_high):
        wavelength_selection = (wavelength >= lower) & (wavelength < upper)
        electrons_in = np.nansum(
            electron_count[wavelength_selection & in_transit[:, None]]
        )
        electrons_out = np.nansum(
            electron_count[wavelength_selection & out_of_transit[:, None]]
        )
        # make sure the counts exist
        if electrons_in > 0 and electrons_out > 0:
            noise = np.sqrt(1.0 / electrons_in + 1.0 / electrons_out)
        else:
            noise = np.nan

        binned_noise.append(noise)
    return np.asarray(binned_noise)


def calculate_jwst_photon_noise(
    source_rate,
    wavelength,
    detector,
    valid,
    in_transit,
    out_of_transit,
    planet_spectrum,
    integration_time,
    gains,
    visits,
):
    """
    Does the actual photon noise calculation by gathering the noise and binning to jwst spectrum
    """
    if visits.shape != detector.shape:
        raise ValueError('Visit labels must have one entry per integration.')

    native_noise = calculate_jwst_native_photon_noise(
        source_rate, detector, valid, integration_time, gains
    )
    binned_photon_noise = {}

    for detector_key, detector_data in planet_spectrum.items():
        detector_name = str(detector_key)
        detector_selection = detector == detector_name

        binned_photon_noise[detector_name] = {}

        for visit, visit_data in detector_data.items():
            # Each fitted spectrum uses only its own visit's integrations.
            visit_selection = visits[detector_selection] == int(visit)
            if not np.any(visit_selection):
                raise ValueError(
                    f'No calibration integrations for {detector_name} visit {visit}.'
                )

            binned_photon_noise[detector_name][visit] = bin_jwst_photon_noise(
                native_noise[detector_name]["electron_count"],
                wavelength[detector_selection],
                in_transit[detector_selection] & visit_selection,
                out_of_transit[detector_selection] & visit_selection,
                visit_data["WBlow"],
                visit_data["WBup"],
            )

    return binned_photon_noise
