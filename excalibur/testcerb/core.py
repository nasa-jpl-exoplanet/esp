'''testcerb core ds'''

# Heritage code shame:
# asdfpylint: disable=invalid-name
# asdfpylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
# import logging

from excalibur.ariel.core import simulate_spectra as ariel_simulate_spectra

# ----------------- --------------------------------------------------
# -- SIMULATE ARIEL SPECTRA ------------------------------------------


def simulate_spectra(target, system_dict, runtime_params, out, verbose=False):
    '''
    Simulate Ariel spectra, adding noise based on the Ariel instrument model
    '''

    status = ariel_simulate_spectra(
        target, system_dict, runtime_params, out, verbose=verbose
    )

    return status
