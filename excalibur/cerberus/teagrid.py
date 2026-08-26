'''cerberus teagrid ds'''

# Heritage code shame:
# pylint: disable=invalid-name

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
        # 'H2CO',
        'H2O',
        'H2',
        'H2S',
        'He',
        'N2',
        'NH3',
        'O3',
        # 'PH3',
        'SO2',
        'HCN',
        'TIO',
        'C2H2',
        'N2O',
        'NO',
        'O2',
        'OH',
    ]

    interp_tea = {}

    temperature = np.load(INTERP_TEA_DIR + 'grid_parameters/temperature.npy')
    pressure = np.load(INTERP_TEA_DIR + 'grid_parameters/pressure.npy')
    XtoH = np.load(INTERP_TEA_DIR + 'grid_parameters/XtoH.npy')
    CtoO = np.load(INTERP_TEA_DIR + 'grid_parameters/CtoO.npy')

    # print('T range', temperature[0], temperature[-1])  # 300-3000
    # print('XtoH range', XtoH[0], XtoH[-1])  # 0.1-100
    # print('CtoO range', CtoO[0], CtoO[-1])  # 0.1-10
    # save these ranges; check whether interp is going outside of range
    # no wait, don't bother.  use the standard edge flags in the interpolator
    # interp_tea['Trange'] = (temperature[0], temperature[-1])

    for molecule in species_name:
        grid_4d = np.load(INTERP_TEA_DIR + molecule + '.npy')
        interp_mol = RegularGridInterpolator(
            (temperature, pressure, XtoH, CtoO),
            grid_4d,
            bounds_error=False,
            fill_value=None,
            method='cubic',
        )
        # temporary drop the cubic spline?  (takes some cpu during debugging)
        interp_tea[molecule] = interp_mol

    return interp_tea
