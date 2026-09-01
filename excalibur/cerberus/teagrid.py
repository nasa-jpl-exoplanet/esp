'''cerberus teagrid ds'''

# Heritage code shame:
# pylint: disable=invalid-name,too-many-locals

import os
import numpy as np
import logging

import excalibur
from excalibur.util.cerberus import calcTEA
from scipy.interpolate import RegularGridInterpolator

log = logging.getLogger(__name__)

INTERP_TEA_DIR = os.path.join(
    excalibur.context['data_dir'], 'CERBERUS/INTERP_TEA/'
)


def get_TEA_grid(modelName=None):

    if modelName is None:
        # log.info('no model name for TEA grid')
        modelDir = INTERP_TEA_DIR
    else:
        # log.info('model name for TEA grid:', modelName)
        modelDir = INTERP_TEA_DIR + modelName + '/'
        if not os.path.isdir(modelDir):
            log.error('TEA grid directory is missing!')
            # print('TEA grid directory is missing!')
            modelDir = INTERP_TEA_DIR

    # (species should be saved and read in; this is a backup default list)
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

    temperature = np.load(modelDir + 'grid_parameters/temperature.npy')
    pressure = np.load(modelDir + 'grid_parameters/pressure.npy')
    XtoH = np.load(modelDir + 'grid_parameters/XtoH.npy')
    CtoO = np.load(modelDir + 'grid_parameters/CtoO.npy')
    if os.path.isfile(modelDir + 'grid_parameters/species.npy'):
        species_name = np.load(modelDir + 'grid_parameters/species.npy')
    else:
        log.warning('TEA species list is missing!')

    # print('T', temperature)
    # print('P', pressure)
    # print('XtoH range', XtoH)
    # print('CtoO range', CtoO)
    # print('T range', temperature[0], temperature[-1])  # 300-3000
    # print('P range', pressure[0], pressure[-1])  #
    # print('XtoH range', XtoH[0], XtoH[-1])  # 0.1-100
    # print('CtoO range', CtoO[0], CtoO[-1])  # 0.1-10
    # print('# of T', len(temperature))
    # print('# of P', len(pressure))
    # print('# of XtoH', len(XtoH))
    # print('# of CtoO', len(CtoO))
    # save these ranges; check whether interp is going outside of range
    # no wait, don't bother.  use the standard edge flags in the interpolator
    # interp_tea['Trange'] = (temperature[0], temperature[-1])

    for molecule in species_name:
        grid_4d = np.load(modelDir + molecule + '.npy')
        # print('grid shape', grid_4d.shape, molecule)
        interp_tea[molecule] = RegularGridInterpolator(
            (temperature, pressure, XtoH, CtoO),
            grid_4d,
            bounds_error=False,
            fill_value=None,
            method='cubic',  # comment out during debugging (linear is faster)
        )

    return interp_tea


def grid_generation(parameters, species, modelName=None, verbose=False):
    # Parameter 1 : temperature in K
    # Parameter 2 : pressure in bar
    # Parameter 3 : metallicities, log10 values
    # Parameter 4 : CtoOs, log10 values
    # Species : list of species names from TEA, e.g. CH4_g, N2_ref

    if modelName is None:
        modelDir = INTERP_TEA_DIR
    else:
        modelDir = INTERP_TEA_DIR + modelName + '/'
        if not os.path.isdir(modelDir):
            os.mkdir(modelDir)
            os.chmod(modelDir, 0o777)  # set permissions to read-write-x
            os.mkdir(modelDir + 'grid_parameters/')
            os.chmod(modelDir + 'grid_parameters/', 0o777)
            # print('creating', modelDir)
        # else:
        #    print(modelDir, 'already exists')

    # save the grid ranges for each parameter
    parameter_names = ['temperature', 'pressure', 'XtoH', 'CtoO']
    for iparam, param in enumerate(parameter_names):
        filename = modelDir + 'grid_parameters/' + param + '.npy'
        np.save(filename, parameters[iparam])
        os.chmod(filename, 0o666)  # set permissions to read-write

    temp_grid, pressure_grid = np.meshgrid(
        parameters[0], parameters[1], indexing='ij'
    )
    temp_comb = temp_grid.flatten()
    pressure_comb = pressure_grid.flatten()

    # species names for calling the cross sections
    #  TIO is a special case; lower-case 'i' is capitalized in the XOMOL dir
    species_name = [
        'TIO' if el == 'TiO_g' else el.split('_')[0] for el in species
    ]
    # save the list of molecules
    filename = modelDir + 'grid_parameters/species.npy'
    np.save(filename, species_name)
    os.chmod(filename, 0o666)  # set permissions to read-write

    # empty initialization
    TEA_mixratio_grids = {
        molecule: np.zeros(
            (
                len(parameters[0]),
                len(parameters[1]),
                len(parameters[2]),
                len(parameters[3]),
            )
        )
        for molecule in species_name
    }

    # call TEA for each XtoH and CtoO
    for i, XtoH in enumerate(parameters[2]):
        for j, CtoO in enumerate(parameters[3]):
            if verbose:
                print('metallicity, [C/O] :', XtoH, CtoO)

            # TEA computation for all the temperatures and pressures
            mixratioprofiles = calcTEA(
                temp_comb,
                pressure_comb,
                species,
                metallicity=10.0**XtoH,
                C_O=0.55 * 10.0**CtoO,
            )
            for molecule in mixratioprofiles:
                TEA_mixratio_grids[molecule][:, :, i, j] = np.reshape(
                    mixratioprofiles[molecule],
                    ((len(parameters[0]), len(parameters[1]))),
                )

    # save the grid for each molecule
    for molecule in TEA_mixratio_grids:
        filename = modelDir + molecule + '.npy'
        np.save(filename, TEA_mixratio_grids[molecule])
        os.chmod(filename, 0o666)  # set permissions to read-write

    return
