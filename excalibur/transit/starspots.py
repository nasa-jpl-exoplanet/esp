'''transit starspots ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=duplicate-code
# pylint: disable=too-many-arguments,too-many-branches,too-many-instance-attributes,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------

from excalibur.transit.spotmodel.Spotmodel import SpotModel
from excalibur.transit.spotmodel.plotters import plot_transit_depths
from excalibur.transit.core import vecistar
from excalibur.cerberus.plotting import rebin_data
from excalibur.util.plotters import add_scale_height_labels
# save_plot_tosv

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import logging

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# def starspots(fin, wht, spc, out):
def starspots(fin, spc, out):
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
    # Lstar = fin['priors']['L*']
    # print('star R,T,M:  ', Rstar, Tstar, Mstar)

    spotssolved = False

    # planetletters = fin['priors']['planets']
    # some of the planets available in fin (system parameters) don't have transits
    # use either transit.whitelight or transit.spectrum for the list of planets
    planetletters = spc['data'].keys()
    for planetletter in planetletters:
        # print('STARSPOTS: big loop over each planet letter', planetletter)

        Rplanet = fin['priors'][planetletter]['rp']
        inc = fin['priors'][planetletter]['inc']
        period = fin['priors'][planetletter]['period']
        sma = fin['priors'][planetletter]['sma']
        # assume zero eccentricity for planet orbit
        ecc = 0
        anom = 0
        # print('planet' + planetletter, 'R,inc,P,a:', Rplanet, inc, period, sma)

        # limb darkening (ignore whitelight, we want the wavelength dependence)
        # limb_coeffs_whitelight = wht['data'][planetletter]['whiteld']
        # print('limb darkening parameters (whitelight)', limb_coeffs_whitelight)
        limb_coeffs = np.array(spc['data'][planetletter]['LD'])
        # print(
        #    'limb darkening parameters from spectrum (median)',
        #    np.median(spc['data'][planetletter]['LD'], axis=0),
        # )
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
        # out['data'][planetletter]['LD'] = limb_coeffs_whitelight
        out['data'][planetletter]['LD'] = limb_coeffs

        # 3) for each planet, calculate starspot model based on the input parameters
        # 4) save the results

        # 5) for each planet, make some plots (spectrum, spotmodel)

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

        # no longer saving the data spectrum; not needed here
        # out['data'][planetletter]['plot_starspot_spectrum'] = save_plot_tosv(
        #    myfig)
        plt.close(myfig)

        # 6) make a plot of the limb darkening as a function of wavelength

        myfig, ax = plt.subplots(figsize=(6, 4))

        radii = np.linspace(0, 1, 111)
        for iwave in range(len(transitdata['wavelength'])):
            # vecistar is normalized such that an integral over the star is 1
            # multiply by pi.R^2 to get the star-averaged limb-darkening
            limbdarkening = (
                vecistar(
                    radii,
                    limb_coeffs[iwave, 0],
                    limb_coeffs[iwave, 1],
                    limb_coeffs[iwave, 2],
                    limb_coeffs[iwave, 3],
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

        # no longer saving the input limbdarkening;
        #  maybe include below alongside the limb-darkened transit profile
        # out['data'][planetletter]['plot_starspot_limbdarkening'] = (
        #     save_plot_tosv(myfig))
        plt.close(myfig)

        # 7) also plot the limb darkening coefficients

        myfig = plt.figure(figsize=(8, 6))

        # print('len check',len(spc['data'][planetletter]['LD']))
        # print('len check',limb_coeffs.shape)
        # print('do these match now?!?',len(transitdata['wavelength']), len(limb_coeffs[:-1,0]))
        ax = myfig.add_subplot(2, 2, 1)
        ax.plot(transitdata['wavelength'], limb_coeffs[:-1, 0], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #1'))
        # plt.title('planet '+planetletter)

        ax = myfig.add_subplot(2, 2, 2)
        ax.plot(transitdata['wavelength'], limb_coeffs[:-1, 1], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #2'))
        ax = myfig.add_subplot(2, 2, 3)

        ax.plot(transitdata['wavelength'], limb_coeffs[:-1, 2], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #3'))

        ax = myfig.add_subplot(2, 2, 4)
        ax.plot(transitdata['wavelength'], limb_coeffs[:-1, 3], color='k')
        plt.xlabel(str('Wavelength [$\\mu$m]'))
        plt.ylabel(str('Limb darkening coeff #4'))
        # plt.title('planet '+planetletter)

        myfig.tight_layout()

        # print('saving plot as testsave.png')
        # plt.savefig('/proj/data/bryden/testsave.png')

        # no longer saving the input limbdarkening coeffs; not that informative'
        # out['data'][planetletter]['plot_limbCoeffs'] = save_plot_tosv(myfig)
        plt.close(myfig)

        spotssolved = True
        out['STATUS'].append(True)

        #    VIKTOR CODE BELOW

        target = planetletter  # Exoplanet name

        # Code Parameters
        max_intensity = 1000
        matrix_radius = 700
        matrix_size = 2 * matrix_radius + 100
        # Geoff: NOT USED
        # timeInterval = 3

        # Plot Options
        plot_anim = False
        plot_star = False
        plot_graph = False

        #  just a few wavelengths during debugging
        # wavelengths = [0.5, 1.0, 1.5]
        wavelengths = transitdata['wavelength']
        num_wavelengths = len(wavelengths)

        # Limb-darkening coefficients
        #  simple case during debugging:
        # c1 = [limb_coeffs_whitelight[0]] * num_wavelengths
        # c2 = [limb_coeffs_whitelight[1]] * num_wavelengths
        # c3 = [limb_coeffs_whitelight[2]] * num_wavelengths
        # c4 = [limb_coeffs_whitelight[3]] * num_wavelengths
        c1 = limb_coeffs[:, 0]
        c2 = limb_coeffs[:, 1]
        c3 = limb_coeffs[:, 2]
        c4 = limb_coeffs[:, 3]

        # Starspots/Faculae
        include_starspots = True  # Caution! Do not change to False
        lat = np.array([20, 20, 20])  # [deg]
        longt = np.array([-20, 0, 20])  # [deg]
        quantidade = len(lat)

        # Spot simulation parameters
        # define a grid of ff_spot and T_spot values
        ff_spot_min = 0.01
        ff_spot_max = 0.1
        T_spot_min = (0.418 * Tstar + 1620) - 1200
        T_spot_max = Tstar - 50
        num_ff_spot_simulations = 4
        num_T_spot_simulations = 5

        # Facula simulation parameters
        # define a grid of ff_fac and T_fac values
        ff_fac_min = 0.01
        ff_fac_max = 0.1
        T_fac_min = Tstar + 50
        T_fac_max = Tstar + 1000
        num_ff_fac_simulations = 4
        num_T_fac_simulations = 5

        other_params = {
            'num_wavelengths': num_wavelengths,
            'c1': c1,
            'c2': c2,
            'c3': c3,
            'c4': c4,
            'lambdaEff': wavelengths,
            'target': target,
            'raio': matrix_radius,
            'intensidadeMaxima': max_intensity,
            'tamanhoMatriz': matrix_size,
            'raioStar': Rstar,
            'ecc': ecc,
            'anom': anom,
            'tempStar': Tstar,
            'starspots': include_starspots,
            'quantidade': quantidade,
            'lat': lat,
            'longt': longt,
            'semiEixoUA': sma,
            'massStar': Mstar,
            'plot_anim': plot_anim,
            'periodo': period,
            'anguloInclinacao': inc,
            'raioPlanetaRj': Rplanet,
            'plot_graph': plot_graph,
            'plot_star': plot_star,
        }

        # 1) RUN THE UNSPOTTED SCENARIO FIRST (ff=0), tempSpot = NaN
        count = 0
        while count == 0:
            include_starspots = False
            # print('Running the unspotted scenario (ff=0) first')

            unspotted_params = other_params.copy()
            unspotted_params['r'] = 0.0  # => ff=0 => unspotted
            # no valid temperature for a spot
            unspotted_params['tempSpot'] = float('nan')

            # run_simulations with ff_min=ff_max=0, T_spot_min=T_spot_max=NaN
            # and only 1 step for each, so it simulates a single unspotted point.
            ff_grid, T_grid, wave_grid, transit_depths_juststar, oneplot = (
                run_simulations(
                    ff_min=0.0,
                    ff_max=0.0,
                    T_spot_min=float('nan'),
                    T_spot_max=float('nan'),
                    num_ff_simulations=1,
                    num_T_spot_simulations=1,
                    other_params=unspotted_params,
                    # result_type="unspotted",
                )
            )
            # print('GRID FOR UNSPOTTED:', count, ff_grid, T_grid)
            count = count + 1
        out['data'][planetletter]['ff_juststar'] = ff_grid
        out['data'][planetletter]['T_juststar'] = T_grid
        out['data'][planetletter]['waves_juststar'] = wave_grid
        out['data'][planetletter]['depths_juststar'] = transit_depths_juststar

        include_starspots = True

        # 2) Run simulations for spots
        ff_grid, T_grid, wave_grid, transit_depths_spots, oneplot = (
            run_simulations(
                ff_spot_min,
                ff_spot_max,
                T_spot_min,
                T_spot_max,
                num_ff_spot_simulations,
                num_T_spot_simulations,
                other_params,
                # "spot",
            )
        )
        # print('GRID FOR SPOTS:', ff_grid, T_grid)
        out['data'][planetletter]['ff_spots'] = ff_grid
        out['data'][planetletter]['T_spots'] = T_grid
        out['data'][planetletter]['waves_spots'] = wave_grid
        out['data'][planetletter]['depths_spots'] = transit_depths_spots

        # 3) Run simulations for faculae
        ff_grid, T_grid, wave_grid, transit_depths_fac, _ = run_simulations(
            ff_fac_min,
            ff_fac_max,
            T_fac_min,
            T_fac_max,
            num_ff_fac_simulations,
            num_T_fac_simulations,
            other_params,
            # "faculae",
        )
        # print('GRID FOR FACULAE:', ff_grid, T_grid)
        out['data'][planetletter]['ff_fac'] = ff_grid
        out['data'][planetletter]['T_fac'] = T_grid
        out['data'][planetletter]['waves_fac'] = wave_grid
        out['data'][planetletter]['depths_fac'] = transit_depths_fac

        # print('ff   ',ff_grid.shape, ff_grid)
        # print('Tspot',T_grid.shape, T_grid)
        # print('waves',wave_grid.shape, wave_grid)
        # print('modelResult depth shape', transit_depths.shape)
        # print('modelResult depth', transit_depths)

        out['data'][planetletter]['plot_starspot_transitdepths'] = (
            plot_transit_depths(
                ff_grid,
                T_grid,
                wave_grid,
                transit_depths_spots,
                transit_depths_juststar,
            )
        )
        out['data'][planetletter]['plot_starspot_deltadepths'] = (
            plot_transit_depths(
                ff_grid,
                T_grid,
                wave_grid,
                transit_depths_spots,
                transit_depths_juststar,
                subtractStar=True,
            )
        )
        # show one of the lightcurve plots (put it after the transitdepth result)
        #  it covers a range of wavelengths for a single ff+Tspot model
        out['data'][planetletter]['plot_starspot_lightcurves'] = oneplot

    return spotssolved


def run_simulations(
    ff_min,
    ff_max,
    T_spot_min,
    T_spot_max,
    num_ff_simulations,
    num_T_spot_simulations,
    other_params,
    # result_type,
):
    """
    Executes simulations for a grid of ff and T_spot
     ff = filling factor
     T_spot = spot temperature

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

    transit_depths = []
    for ff in grid_ff:
        transit_depths.append([])

        if not 0 <= ff <= 1:
            raise ValueError(f"ff={ff} is out of range [0, 1].")

        # If there are multiple starspots, each receives a proportional filling factor
        # Ensure it's at least 1
        quantidade = max(other_params.get('quantidade', 1), 1)
        # Distribute the filling factor equally among the spots
        ff_per_spot = ff / quantidade
        # Convert to individual spot radius
        spot_radius = np.sqrt(ff_per_spot)

        for T_spot in grid_T_spot:
            transit_depths[-1].append([])

            iteration_params = other_params.copy()
            iteration_params['r'] = spot_radius
            iteration_params['tempSpot'] = T_spot

            # print(
            #     f"Running simulation for {result_type} with ff={ff} and T={T_spot}"
            # )

            # executes the SpotModel with the chosen parameters (ff and T_spot)
            oneModel = SpotModel(iteration_params)

            # print('modelResult ff', oneModel.ff)
            # print('modelResult T', oneModel.T)
            # print('modelResult wave', oneModel.wavearray)
            # print('modelResult depth', oneModel.depth)
            # depths as a func of wavelength
            transit_depths[-1][-1] = oneModel.depth
            # print('transitdepths',transit_depths)

    # print('transitdepths before arraying',transit_depths)
    # print('transitdepths before arraying',len(transit_depths))
    # print('transitdepths before arraying',len(transit_depths[0]))
    # print('transitdepths before arraying',len(transit_depths[0][0]))
    return (
        grid_ff,
        grid_T_spot,
        oneModel.wavearray,
        np.array(transit_depths),
        oneModel.plot_lightcurves,
    )
