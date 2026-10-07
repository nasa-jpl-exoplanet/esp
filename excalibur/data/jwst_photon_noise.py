
import numpy as np

def extract_jwst_photon_noise(source_rate, integration_time, gain):
    """
    Estimates the photon noise from JWST observations
    """
    electron_count = source_rate * integration_time * gain
    electron_count = np.where(electron_count>0,electron_count,np.nan) # change negative values to nan
    photon_noise_abs = np.sqrt(electron_count)
    photon_noise_rel = 1.0 / np.sqrt(electron_count)

    output_dict = {
        "electron_count": electron_count,
        "photon_noise_abs": photon_noise_abs,
        "photon_noise_rel": photon_noise_rel
    }

    return output_dict

def bin_jwst_photon_noise(electron_count, wavelength, in_transit, out_of_transit, wavelength_low, wavelength_high):
    """
    Bin electron counts and calculate photon limited transit depth errors
    """
    binned_noise = []
    for lower, upper in zip(wavelength_low, wavelength_high):
        wavelength_selection = (wavelength >= lower) & (wavelength < upper)
        electrons_in = np.nansum(electron_count[wavelength_selection & in_transit[:,None]])
        electrons_out = np.nansum(electron_count[wavelength_selection & out_of_transit[:,None]])
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
):
    """
    Does the actual photon noise calculation by gathering the noise and binning to jwst spectrum
    """
    binned_photon_noise = {}

    for detector_key, detector_data in planet_spectrum.items():
        detector_name = str(detector_key)
        detector_selection = detector == detector_name

        detector_rate = np.where(
            valid[detector_selection],
            source_rate[detector_selection],
            np.nan,
        )

        if isinstance(integration_time, dict):
            detector_integration_time = integration_time[detector_name]
        else:
            detector_integration_time = integration_time

        photon_noise = extract_jwst_photon_noise(
            detector_rate,
            detector_integration_time,
            gains[detector_name],
        )

        binned_photon_noise[detector_name] = {}

        for visit, visit_data in detector_data.items():
            binned_photon_noise[detector_name][visit] = (
                bin_jwst_photon_noise(
                    photon_noise["electron_count"],
                    wavelength[detector_selection],
                    in_transit[detector_selection],
                    out_of_transit[detector_selection],
                    visit_data["WBlow"],
                    visit_data["WBup"],
                )
            )

    return binned_photon_noise
