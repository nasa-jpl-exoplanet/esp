'''cerberus teagrid ds'''

# Heritage code shame:
# pylint: disable=invalid-name

import os
import numpy as np
import excalibur
from excalibur.util.cerberus import calcTEA
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
        )
        # temporary drop the cubic spline.  takes some cpu during debugging
        #    method='cubic',
        interp_tea[molecule] = interp_mol

    return interp_tea


def grid_generation(parameters, species, verbose=False):
    # Parameter 1 : temperature in K
    # Parameter 2 : pressure in bar
    # Parameter 3 : metallicities, log10 values
    # Parameter 4 : CtoOs, log10 values
    # Species : list of species with the names from TEA
    # ex : ['CH4_g, N2_ref, ...]

    temp = parameters[:, 0]
    pressure = parameters[:, 1]
    metallicities = parameters[:, 2]
    CtoOs = parameters[:, 3]

    # parameters saving
    np.save(INTERP_TEA_DIR + 'grid_parameters/temperature.npy', temp)
    np.save(INTERP_TEA_DIR + 'grid_parameters/pressure.npy', pressure)
    np.save(INTERP_TEA_DIR + 'grid_parameters/XtoH.npy', metallicities)
    np.save(INTERP_TEA_DIR + 'grid_parameters/CtoO.npy', CtoOs)

    temp_grid, pressure_grid = np.meshgrid(temp, pressure, indexing='ij')
    temp_comb = temp_grid.flatten()
    pressure_comb = pressure_grid.flatten()

    # species names for calling the cross sections
    # TIO for the cross section library
    species_name = [
        'TIO' if el == 'TiO_g' else el.split('_')[0] for el in species
    ]

    # dictionary creation
    dic = {
        molecule: np.zeros(
            (len(temp), len(pressure), len(metallicities), len(CtoOs))
        )
        for molecule in species_name
    }

    # call TEA for each XtoH and CtoO
    for XtoH in metallicities:
        for CtoO in CtoOs:
            if verbose:
                # to allow to know at what steps we are for the computation
                print("metallicity :", 10.0**XtoH)
                print("C/O :", 10.0**CtoO)
            # TEA computation for all the temperatures and pressures
            mixratioprofiles = calcTEA(
                temp_comb,
                pressure_comb,
                species,
                metallicity=10.0**XtoH,
                C_O=0.55 * 10.0**CtoO,
            )
            for molecule in mixratioprofiles:
                dic[molecule][:, :, i, j] = np.reshape(
                    mixratioprofiles[molecule], ((len(temp), len(pressure)))
                )

    # grid saving for each molecule
    for molecule in mon_dico:
        np.save(INTERP_TEA_DIR + molecule + '.npy', dic[molecule])
