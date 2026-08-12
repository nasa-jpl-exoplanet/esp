'''cerberus teagrid ds'''

# Heritage code shame:
# asfpylint: disable=invalid-name
# asdfpylint: disable=no-member

import os
import numpy as np
import excalibur
from scipy.interpolate import RegularGridInterpolator

INTERP_TEA_DIR = os.path.join(
    excalibur.context['data_dir'], 'CERBERUS/INTERP_TEA/'
)


def get_TEA_grid():

    species_name = [
        'CH4',
        'CO2',
        'CO',
        'H2O',
        'H2',
        'H2S',
        'He',
        'N2',
        'NH3',
        'O3',
        'SO2',
        'HCN',
        'TIO',
        'C2H2',
        'N2O',
        'NO',
        'O2',
        'OH',
    ]

    temp = np.load(INTERP_TEA_DIR + 'grid_parameters/temperature.npy')
    pressure = np.load(INTERP_TEA_DIR + 'grid_parameters/pressure.npy')
    XtoH = np.load(INTERP_TEA_DIR + 'grid_parameters/XtoH.npy')
    CtoO = np.load(INTERP_TEA_DIR + 'grid_parameters/CtoO.npy')

    interp_tea = {}
    for molecule in species_name:
        grid_4d = np.load(INTERP_TEA_DIR + molecule + '.npy')
        interp_mol = RegularGridInterpolator(
            (temp, pressure, XtoH, CtoO), grid_4d
        )
        # temporary drop the cubic spline.  takes some cpu during debugging
        # interp_mol = RegularGridInterpolator(
        #    (temp, pressure, XtoH, CtoO), grid_4d, method='cubic'
        # )
        interp_tea[molecule] = interp_mol

    return interp_tea
