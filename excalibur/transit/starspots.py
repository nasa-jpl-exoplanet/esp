'''transit starspots ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=duplicate-code
# pylint: disable=too-many-arguments,too-many-branches,too-many-instance-attributes,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
# import dawgie

# import excalibur.data.core as datcore
# import excalibur.system.core as syscore
# import excalibur.util.cerberus as crbutil
from excalibur.cerberus.plotting import rebin_data
from excalibur.util.plotters import (
    save_plot_tosv,
    #    save_plot_myfit,
    #    plot_residual_fft,
    add_scale_height_labels,
)
from excalibur.transit.core import vecistar
from excalibur.transit.spotmodel.spotmodel import SpotModel

# import pandas as pd
# from matplotlib import pyplot

# import os
# import main
# import importlib
# importlib.reload(main)
# from scipy.interpolate import griddata

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import logging

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
def starspots(fin, wht, spc, out):
    '''
    Viktor Sumida's starspot model
    '''

    # print('whitelight planetletters',wht['data'].keys())
    # print('spectrum planetletters',spc['data'].keys())
    # print('whitelight info',wht['data']['b'].keys())
    # print('spectrum info',spc['data']['b'].keys())

    # 1) make sure we have all the input parameters

    # print('fin keys',fin['priors'].keys())

    Rstar = fin['priors']['R*']
    Tstar = fin['priors']['T*']
    Mstar = fin['priors']['M*']
    Lstar = fin['priors']['L*']
    print('star R,T,L,M:  ', Rstar, Tstar, Lstar, Mstar)

    exospec = False

    # planetletters = fin['priors']['planets']
    # some of the planets available in fin (system parameters) don't have transits
    # use either transit.whitelight or transit.spectrum for the list of planets
    planetletters = spc['data'].keys()
    for planetletter in planetletters:
        print('STARSPOTS: big loop over each planet letter', planetletter)

        Rplanet = fin['priors'][planetletter]['rp']
        inc = fin['priors'][planetletter]['inc']
        period = fin['priors'][planetletter]['period']
        sma = fin['priors'][planetletter]['sma']
        print('planet' + planetletter, 'R,inc,P,a:', Rplanet, inc, period, sma)

        # limb dark
        limb_coeffs = wht['data'][planetletter]['whiteld']
        print('limb darkening parameters from whitelight       ', limb_coeffs)
        print(
            'limb darkening parameters from spectrum (median)',
            np.median(spc['data'][planetletter]['LD'], axis=0),
        )
        # print('limb darkening parameters from spectrum (mean)  ',
        #      np.mean(spc['data'][planetletter]['LD'],axis=0))

        # this is the spectrum without any starspot correction
        transitdata = {}
        transitdata['wavelength'] = spc['data'][planetletter]['WB']
        transitdata['depth'] = spc['data'][planetletter]['ES'] ** 2
        transitdata['error'] = (
            2
            * spc['data'][planetletter]['ES']
            * spc['data'][planetletter]['ESerr']
        )

        # for key in spc['data'][planetletter].keys():
        #    if key!='Teq':
        #        print('len check on spc',key)
        #        print('len check on spc',key,len(spc['data'][planetletter][key]))

        # bins the data (for plotting only, not for science analysis)
        transitdata = rebin_data(transitdata)

        # 2) save whatever inputs might be helpful later
        out['data'][planetletter] = {}
        out['data'][planetletter]['WB'] = transitdata['wavelength']
        out['data'][planetletter]['ES'] = transitdata['depth']
        out['data'][planetletter]['ESerr'] = transitdata['error']
        out['data'][planetletter]['RSTAR'] = Rstar
        out['data'][planetletter]['TSTAR'] = Tstar
        out['data'][planetletter]['LD'] = limb_coeffs

        # 3) for each planet, calculate starspot model based on the input parameters

        #  ************************************************************
        #               Viktor - insert your code here
        #  ************************************************************

        # I think that your result will be an array of spectra on a 2-D grid
        #  one axis is filling factor and
        #  one axis is the spot temperature
        # ?????
        spotmodel = np.zeros((100, 100, 100))

        # 4) save the results
        out['data'][planetletter]['spotmodel'] = spotmodel

        # 5) for each planet, make a plot of the spectrum

        myfig, ax = plt.subplots(figsize=(8, 6))

        ax.errorbar(
            transitdata['wavelength'],
            1e2 * transitdata['depth'],
            fmt='.',
            yerr=1e2 * transitdata['error'],
            color='lightgray',
        )
        ax.errorbar(
            transitdata['binned_wavelength'],
            1e2 * transitdata['binned_depth'],
            fmt='o',
            yerr=1e2 * transitdata['binned_error'],
            color='k',
        )

        plt.title('planet ' + planetletter)
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('$(R_p/R_*)^2$ [%]'))
        add_scale_height_labels(
            spc['data'][planetletter], transitdata['depth'], ax, myfig
        )

        # print('saving plot as testsave.png')
        # plt.savefig('/proj/data/bryden/testsave.png')

        out['data'][planetletter]['plot_starspot_spectrum'] = save_plot_tosv(
            myfig
        )
        plt.close(myfig)

        # 6) make a plot of the limb darkening as a function of wavelength

        myfig, ax = plt.subplots(figsize=(6, 4))

        LD = np.array(spc['data'][planetletter]['LD'])
        radii = np.linspace(0, 1, 111)
        for iwave in range(len(transitdata['wavelength'])):
            # vecistar is normalized such that an integral over the star is 1
            # multiply by pi.R^2 to get the star-averaged limb-darkening
            limbdarkening = (
                vecistar(
                    radii,
                    LD[iwave, 0],
                    LD[iwave, 1],
                    LD[iwave, 2],
                    LD[iwave, 3],
                )
                * np.pi
            )
            wavelength_colors = mpl.cm.rainbow(
                float(iwave) / (len(transitdata['wavelength']) - 1)
            )
            ax.plot(radii, limbdarkening, color=wavelength_colors)

        plt.xlabel('Radius [$R_{\\star}$]', fontsize=14)
        # plt.ylabel(str('$I(R)/I(0)$'),fontsize=14)
        plt.ylabel('limb darkening (normalized)', fontsize=14)
        plt.title('planet ' + planetletter, fontsize=14)
        plt.xlim([0, 1])
        ylims = ax.get_ylim()
        # plt.ylim([0,ylims[1]])
        plt.ylim([0.5, ylims[1]])

        myfig.tight_layout()

        #        print('saving plot as testsave.png')
        #        plt.savefig('/proj/data/bryden/testsave.png')

        out['data'][planetletter]['plot_starspot_limbdarkening'] = (
            save_plot_tosv(myfig)
        )
        plt.close(myfig)

        # 7) also plot the limb darkening coefficients

        myfig = plt.figure(figsize=(8, 6))

        LD = np.array(spc['data'][planetletter]['LD'])
        # print('len check',len(spc['data'][planetletter]['LD']))
        # print('len check',LD.shape)
        # print('do these match now?!?',len(transitdata['wavelength']), len(LD[:-1,0]))
        ax = myfig.add_subplot(2, 2, 1)
        ax.plot(transitdata['wavelength'], LD[:-1, 0], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #1'))
        # plt.title('planet '+planetletter)

        ax = myfig.add_subplot(2, 2, 2)
        ax.plot(transitdata['wavelength'], LD[:-1, 1], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #2'))
        ax = myfig.add_subplot(2, 2, 3)

        ax.plot(transitdata['wavelength'], LD[:-1, 2], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #3'))

        ax = myfig.add_subplot(2, 2, 4)
        ax.plot(transitdata['wavelength'], LD[:-1, 3], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #4'))
        # plt.title('planet '+planetletter)

        myfig.tight_layout()

        # print('saving plot as testsave.png')
        # plt.savefig('/proj/data/bryden/testsave.png')

        out['data'][planetletter]['plot_limbCoeffs'] = save_plot_tosv(myfig)
        plt.close(myfig)

        exospec = True
        out['STATUS'].append(True)

        #    VIKTOR CODE BELOW

        target = planetletter  # Exoplanet name

        # Code Parameters

        max_intensity = 1000
        matrix_radius = 700
        matrix_size = 2 * matrix_radius + 100
        # Geoff: NOT USED
        # timeInterval = 3
        # en to pt
        intensidadeMaxima = max_intensity
        raio = matrix_radius
        tamanhoMatriz = matrix_size

        # Main Parameters
        # Star
        # raioStar = 0.2159  # [R_sun]
        # massStar = 0.1844  # [M_sun]
        # tempStar = 3100  # [K]
        raioStar = Rstar  # [R_sun]
        massStar = Mstar  # [M_sun]
        tempStar = Tstar  # [K]

        # Planet
        # raioPlanetaRj = 0.1543  # [Rj] (in in jupiter's radius)
        # periodo = 24.73723  # [days]
        # anguloInclinacao = 89.86  # [deg]
        # semiEixoUA = 0.0946  # [AU]
        raioPlanetaRj = Rplanet  # [Rj] (in in jupiter's radius)
        periodo = period  # [days]
        anguloInclinacao = inc  # [deg]
        semiEixoUA = sma  # [AU]
        ecc = 0
        anom = 0

        # Plot Options
        plot_anim = False
        plot_star = False
        plot_graph = False

        # Limb-darkening coefficients and wavelengths
        #  Geoff: we're only considering a single set of LD coeff
        #         so maybe clean this up later and remove num_elements
        #          and remove the loop inside of spotmodel/main
        c1 = [limb_coeffs[0]]
        c2 = [limb_coeffs[1]]
        c3 = [limb_coeffs[2]]
        c4 = [limb_coeffs[3]]
        num_elements = len(c1)
        # set the limb-darkening profile
        profile = '4-parameter'  # Geoff:  *** ??!?!? ***

        lambdaEff = [0.5, 1.0, 1.5]  # <-- fix this later after testing!

        # Starspots/Faculae
        include_starspots = True  # Caution! Do not change to False
        lat = np.array([20, 20, 20])  # [deg]
        longt = np.array([-20, 0, 20])  # [deg]
        quantidade = len(lat)

        # Spot simulation parameters
        # define a grid of ff_spot and T_spot values
        ff_spot_min = 0.01
        ff_spot_max = 0.1
        T_spot_min = (0.418 * tempStar + 1620) - 1200
        T_spot_max = tempStar - 50
        num_ff_spot_simulations = 5
        num_T_spot_simulations = 5

        # Facula simulation parameters
        # define a grid of ff_fac and T_fac values
        ff_fac_min = 0.01
        ff_fac_max = 0.1
        T_fac_min = tempStar + 50
        T_fac_max = tempStar + 1000
        num_ff_fac_simulations = 5
        num_T_fac_simulations = 5

        other_params = {
            'num_elements': num_elements,
            'profile': profile,
            'c1': c1,
            'c2': c2,
            'c3': c3,
            'c4': c4,
            'lambdaEff': lambdaEff,
            'target': target,
            'raio': raio,
            'intensidadeMaxima': intensidadeMaxima,
            'tamanhoMatriz': tamanhoMatriz,
            'raioStar': raioStar,
            'ecc': ecc,
            'anom': anom,
            'tempStar': tempStar,
            'starspots': include_starspots,
            'quantidade': quantidade,
            'lat': lat,
            'longt': longt,
            'semiEixoUA': semiEixoUA,
            'massStar': massStar,
            'plot_anim': plot_anim,
            'periodo': periodo,
            'anguloInclinacao': anguloInclinacao,
            'raioPlanetaRj': raioPlanetaRj,
            'plot_graph': plot_graph,
            'plot_star': plot_star,
        }

        # 1) RUN THE UNSPOTTED SCENARIO FIRST (ff=0), tempSpot = NaN
        count = 0
        while count == 0:
            include_starspots = False
            print(
                "\nRunning the unspotted scenario (ff=0) first, so it appears first in the results file."
            )

            unspotted_params = other_params.copy()
            unspotted_params['r'] = 0.0  # => ff=0 => unspotted
            unspotted_params['tempSpot'] = float(
                'nan'
            )  # no valid temperature for a spot

            # Here we call run_simulations with ff_min=0, ff_max=0, T_spot_min=NaN, T_spot_max=NaN
            # and only 1 step for each, so it simulates a single unspotted point.
            run_simulations(
                ff_min=0.0,
                ff_max=0.0,
                T_spot_min=float('nan'),
                T_spot_max=float('nan'),
                num_ff_simulations=1,
                num_T_spot_simulations=1,
                other_params=unspotted_params,
                result_type="unspotted",
            )
            count = count + 1
        include_starspots = True

        # 2) Run simulations for spots
        spotresults = run_simulations(
            ff_spot_min,
            ff_spot_max,
            T_spot_min,
            T_spot_max,
            num_ff_spot_simulations,
            num_T_spot_simulations,
            other_params,
            "spot",
        )
        print('RESULTS FOR SPOTS:', spotresults)

        # 3) Run simulations for faculae
        facresults = run_simulations(
            ff_fac_min,
            ff_fac_max,
            T_fac_min,
            T_fac_max,
            num_ff_fac_simulations,
            num_T_fac_simulations,
            other_params,
            "faculae",
        )
        print('RESULTS FOR Faculae:', facresults)

    return exospec


def run_simulations(
    ff_min,
    ff_max,
    T_spot_min,
    T_spot_max,
    num_ff_simulations,
    num_T_spot_simulations,
    other_params,
    result_type,
):
    """
    Executes simulations for a grid of ff and T_spot values.

    Parameters:
        ff_min (float): Minimum value for ff.
        ff_max (float): Maximum value for ff.
        T_spot_min (float): Minimum value for T_spot.
        T_spot_max (float): Maximum value for T_spot.
        num_ff_simulations (int): Number of simulations in ff.
        num_T_spot_simulations (int): Number of simulations in T_spot.
        other_params (dict): Fixed parameters for the main program execution.
        result_type (str): Type of the simulation ("spot", "faculae", or "unspotted").

    Returns:
        (simulated_ff, simulated_T_spot): Arrays of ff and T_spot used in the simulations.
    """
    grid_ff = np.linspace(ff_min, ff_max, num_ff_simulations)
    grid_T_spot = np.linspace(T_spot_min, T_spot_max, num_T_spot_simulations)

    for ff in grid_ff:
        for T_spot in grid_T_spot:
            if not 0 <= ff <= 1:
                raise ValueError(f"ff={ff} is out of range [0, 1].")

            # If there are multiple starspots, each spot receives a proportional filling factor
            quantidade = max(
                other_params.get('quantidade', 1), 1
            )  # Ensure it's at least 1
            ff_per_spot = (
                ff / quantidade
            )  # Distribute the filling factor equally among the spots
            spot_radius = np.sqrt(
                ff_per_spot
            )  # Convert to individual spot radius

            iteration_params = other_params.copy()
            iteration_params['r'] = spot_radius
            iteration_params['tempSpot'] = T_spot

            print(
                f"Running simulation for {result_type} with ff={ff} and T={T_spot}"
            )

            # Executes the SpotModel with the updated parameters
            SpotModel(
                target=other_params['target'],
                num_elements=iteration_params['num_elements'],
                profile=iteration_params['profile'],
                c1=iteration_params['c1'],
                c2=iteration_params['c2'],
                c3=iteration_params['c3'],
                c4=iteration_params['c4'],
                lambdaEff=iteration_params['lambdaEff'],
                raio=iteration_params['raio'],
                intensidadeMaxima=iteration_params['intensidadeMaxima'],
                tamanhoMatriz=iteration_params['tamanhoMatriz'],
                raioStar=iteration_params['raioStar'],
                ecc=iteration_params['ecc'],
                anom=iteration_params['anom'],
                tempStar=iteration_params['tempStar'],
                starspots=iteration_params['starspots'],
                quantidade=iteration_params['quantidade'],
                lat=iteration_params['lat'],  # Mantido como array
                longt=iteration_params['longt'],  # Mantido como array
                r=iteration_params['r'],
                semiEixoUA=iteration_params['semiEixoUA'],
                massStar=iteration_params['massStar'],
                plot_anim=iteration_params['plot_anim'],
                periodo=iteration_params['periodo'],
                anguloInclinacao=iteration_params['anguloInclinacao'],
                raioPlanetaRj=iteration_params['raioPlanetaRj'],
                plot_graph=iteration_params['plot_graph'],
                plot_star=iteration_params['plot_star'],
                tempSpot=iteration_params['tempSpot'],
            )

    return grid_ff, grid_T_spot
