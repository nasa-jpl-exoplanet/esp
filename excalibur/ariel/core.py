'''ariel core ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements,too-many-arguments,too-many-positional-arguments,too-many-function-args

# -- IMPORTS -- ------------------------------------------------------
import logging

from collections import namedtuple

# import excalibur
import excalibur.system.core as syscore
import excalibur.util.cerberus as crbutil

from excalibur.ariel.metallicity import (
    massMetalRelation,
    massMetalRelationDisp,
    randomCtoO_linear,
)
from excalibur.ariel.clouds import fixedCloudParameters, randomCloudParameters
from excalibur.ariel.ariel_instrument_model import load_ariel_instrument

# from excalibur.ariel.ariel_instrument_model import (
#    load_ariel_instrument,
#    calculate_ariel_instrument,
# )
from excalibur.ariel.forward_models import make_cerberus_atmos
from excalibur.cerberus.core import myxsecs
from excalibur.ariel.plotters import (
    plot_spectrum,
    plot_spectrum_topmolecules,
    plot_depthprobed,
)

# import os
# import sys

# import pickle
import numpy as np
import scipy.constants as cst

log = logging.getLogger(__name__)

ArielParams = namedtuple(
    'ariel_params_from_runtime',
    [
        'tier',
        'arielRad',
        'SNRfactor',
        'randomSeed',
        'randomCloudProperties',
        'thorngrenMassMetals',
        'chachanMassMetals',
        'includeMetallicityDispersion',
        'metallicityDispersion',
        'CtoOdaSilva',
        'CtoOaverage',
        'CtoOdispersion',
        'knownspecies',
        'cialist',
        'xmollist',
        'nlevels',
        'solrad',
        'Hsmax',
        'lbroadening',
        'lshifting',
        'isothermal',
    ],
)


# ----------------- --------------------------------------------------
# -- SIMULATE ARIEL SPECTRA ------------------------------------------
def calc_mmw_Hs(pressureArray, temperature, logg, X2Hr=0, useTEA=False):
    '''
    calculate the mean molecular weight and scale height
    '''

    # INCLUDE C/O RATIO HERE????  ASDF

    if useTEA:
        # log.error('TEA removed for now')
        tempCoeffs = [0, temperature, 0, 1, 0, -1, 1, 0, -1, 1]  # isothermal
        mixratioprofiles = crbutil.calcTEA(
            tempCoeffs, pressureArray, metallicity=10.0**X2Hr
        )
        # have to take the average! (same as done in crbce)
        mixratio = {}
        for molecule in mixratioprofiles:
            mixratio[molecule] = np.log10(
                np.mean(10.0 ** mixratioprofiles[molecule])
            )

        mmw, fH2, fHe = crbutil.getmmw(mixratio)

        # (uncomment to compare TEC vs TEA)
        # print('TEA:', mixratio, fH2, fHe)
        # print('TEA: water profile', mixratioprofiles['H2O'])
        # mixratio, mixratioprofiles, fH2, fHe = crbutil.crbce(
        #    pressureArray, temperature, X2Hr=X2Hr
        # )
        # print('TEC:', mixratio, fH2, fHe)
        # print('TEC: water profile', mixratioprofiles['H2O'])

    else:
        # mixratio, mixratioprofiles, fH2, fHe = crbutil.crbce(
        mixratio, _, fH2, fHe = crbutil.crbce(
            pressureArray, temperature, X2Hr=X2Hr
        )
        mmw, fH2, fHe = crbutil.getmmw(
            mixratio,
            protosolar=False,
            fH2=fH2,
            fHe=fHe,
        )
    # print('mmw      (inside)', mmw, fH2, fHe)

    # print('mixratio (inside)', mixratio, fH2, fHe)
    # X2Hr=cheq['XtoH'])
    # assume solar C/O and N/O for now
    # C2Or=cheq['CtoO'], N2Or=cheq['NtoO'])

    mmw_kg = mmw * cst.m_p  # [kg]
    Hs = (
        cst.Boltzmann * temperature / mmw_kg / 1e-2 / (10.0 ** float(logg))
    )  # [m]

    return mmw, Hs


def simulate_spectra(
    target,
    system_dict,
    ancil_dict,
    runtime_params,
    out,
    verbose=False,
):
    '''
    Simulate Ariel spectra, adding noise based on the Ariel instrument model
    Mulitple spectra are now calculated, allowing a choice within cerberus.atmos fitting
    1) only Cerberus atmosphere models now; Taurex option removed Dec.2024
    2) with or without clouds
    3) two models for metallicity/mmw (mmw=2.3 or FINESSE mass-metallicity relation)
    4) TEC vs TEA (DISEQ is not implemented yet)
    '''
    # print(runtime_params)
    # print('metallicity dispersion?',runtime_params.includeMetallicityDispersion)

    testTarget = bool(target.startswith('test'))

    # select Tier-1 or Tier-2 for spectra SNR
    tier = runtime_params.tier

    sscmks = syscore.ssconstants(mks=True)

    system_params = system_dict['priors']
    ancil_params = ancil_dict['data']

    # maybe save/read the xslib file from disk.
    # for debugging at least, since it's a lot faster that way
    # xslibSaveDir = os.path.join(excalibur.context['data_dir'], 'bryden/')

    # specify which models should be calculated (use these as keys within data)
    atmosModels = [
        'cerberus',
        'cerberusTEA',
        'cerberuslowmmw',
        'cerberusNoclouds',
        'cerberusTEANoclouds',
        'cerberuslowmmwNoclouds',
    ]

    if testTarget:
        atmosModels = [
            'cerberus',
            'cerberusTEA',
            'cerberusNoclouds',
            'cerberusTEANoclouds',
        ]

    solarCtoO = 0.54951

    out['data']['models'] = atmosModels
    # save target,planet names, for plotting (in states.py)
    out['data']['target'] = target
    out['data']['planets'] = system_params['planets']

    # save the cross-section library from one model to the next (no need to re-calculate)
    # but be careful, it has to be defined for each planet letter
    xslib = {'data': {}, 'STATUS': [False]}

    completed_at_least_one_planet = False
    for planet_letter in system_params['planets']:
        out['data'][planet_letter] = {}

        # set the random seed as a function of target name
        #  (otherwise all spectra have identical noise realizations!)
        #  (also they would have the same C/O and metallicity offset!)
        intFromTarget = 1  # arbitrary initialization for the random seed
        # for the test targets, use a new set of noise for each target
        if testTarget:
            intFromTarget += int(target[-3:])
        for char in target + ' ' + planet_letter:
            intFromTarget = (
                runtime_params.randomSeed * intFromTarget + ord(char)
            ) % 1000000
        np.random.seed(intFromTarget)

        # check for Ariel targets that are not in excalibur
        # from excalibur.target.targetlists import targetlist_ArielMCSknown_transitCategory
        # targs = targetlist_ArielMCSknown_transitCategory()
        # import excalibur.target.edit as trgedit
        # excaliburTargets = trgedit.targetlist.__doc__
        # for targ in targs:
        #     if targ[:-2] not in excaliburTargets:
        #         print('ADD NEW TARGET:',targ)

        # load in the wavelength bins and the noise model
        # there is a separate SNR file for each planet
        targetplanet = target + ' ' + planet_letter
        if testTarget:
            # use HD 209458 SNR as a default for test cases
            targetplanet = 'HD 209458 b'

        # select old ArielRad (results read from a table) or new (results calculated internally)
        oldArielRad = True
        # oldArielRad = False
        if oldArielRad:
            ariel_instrument = load_ariel_instrument(
                targetplanet,
                system_params,
                ancil_params,
                runtime_params,
            )
        else:
            # ariel_instrument = calculate_ariel_instrument(
            ariel_instrument = load_ariel_instrument(
                targetplanet,
                system_params,
                ancil_params,
                runtime_params,
                # verbose=verbose,  # put this back in for calculate_ariel_instrument
            )

        if ariel_instrument:
            # asdf : LATER : add in uncertainty scatter to these model parameters
            model_params = {
                'T*': system_params['T*'],
                'R*': system_params['R*'],
                'Rp': system_params[planet_letter]['rp'],
                'Mp': system_params[planet_letter]['mass'],
                'logg': system_params[planet_letter]['logg'],
                'sma': system_params[planet_letter]['sma'],
                'Teq': system_params[planet_letter]['teq'],
            }

            # Calculate the atmosphere scale height
            #  cerberus wants it, to normalize the spectrum
            #  (some of this code is copied from transit/core)
            #  NOTE : sma is system, but should be model value.  EDIT LATER?
            eqtemp = model_params['T*'] * np.sqrt(
                model_params['R*']
                * sscmks['Rsun/AU']
                / (2.0 * model_params['sma'])
            )
            pgrid = np.exp(
                np.arange(
                    np.log(runtime_params.solrad) - runtime_params.Hsmax,
                    np.log(runtime_params.solrad)
                    + runtime_params.Hsmax / runtime_params.nlevels,
                    runtime_params.Hsmax / (runtime_params.nlevels - 1),
                )
            )
            pressure = pgrid[::-1]

            # planet metallicity is from an assumed mass-metallicity relation
            #  with scatter included
            # planet metallicity should be defined relative to the stellar metallicity
            metallicity_star_dex = system_params['FEH*']
            M_p = model_params['Mp']
            if runtime_params.includeMetallicityDispersion:
                # make sure that the random mass is fixed for each target planet
                np.random.seed(intFromTarget + 1234)
                metallicity_planet_dex = massMetalRelationDisp(
                    metallicity_star_dex,
                    M_p,
                    thorngren=runtime_params.thorngrenMassMetals,
                    chachan=runtime_params.chachanMassMetals,
                    dispersion=runtime_params.metallicityDispersion,
                )
            else:
                metallicity_planet_dex = massMetalRelation(
                    metallicity_star_dex,
                    M_p,
                    thorngren=runtime_params.thorngrenMassMetals,
                    chachan=runtime_params.chachanMassMetals,
                )
            # print('metallicity_star_dex',metallicity_star_dex)
            # print('metallicity_planet_dex',metallicity_planet_dex)
            # metallicity_planet_dex_nonrandom = massMetalRelation(metallicity_star_dex, M_p,
            #              thorngren=runtime_params.thorngrenMassMetals)
            #              chachan=runtime_params.chachanMassMetals)
            # print('metallicity_planet_dex (non random)',metallicity_planet_dex_nonrandom)
            # print('planet mass',M_p)

            # make sure that the random C/O is fixed for each target planet
            np.random.seed(intFromTarget + 12345)
            # this is linear C/O, not dex

            # option to use da Silva 2024 C/O trend as the baseline,
            #  (before adding on some dispersion)
            if runtime_params.CtoOdaSilva:
                CtoOstar = ancil_params['CO*']
                CtoO_planet_linear = randomCtoO_linear(
                    logCtoOaverage=CtoOstar + np.log10(solarCtoO),
                    logCtoOdispersion=runtime_params.CtoOdispersion,
                )
                # oldCtoO_planet_linear = randomCtoO_linear(
                #    logCtoOaverage=runtime_params.CtoOaverage,
                #    logCtoOdispersion=runtime_params.CtoOdispersion,
                # )
                # print('CtoO_planet_linear new,old',
                #      CtoO_planet_linear,oldCtoO_planet_linear)
            else:
                CtoO_planet_linear = randomCtoO_linear(
                    logCtoOaverage=runtime_params.CtoOaverage,
                    logCtoOdispersion=runtime_params.CtoOdispersion,
                )
                # print('CtoO_planet_linear',CtoO_planet_linear)

            # Load the instrument model and rescale based on #-of-transits
            uncertainties = ariel_instrument['noise']
            # use #-of-visits our ArielRad calculation, not from Edwards table
            visits = ariel_instrument['nVisits']
            # print('# of visits:',visits,'  tier',tier,'  ',target+' '+planet_letter)

            uncertainties /= np.sqrt(float(visits))

            # allow for arbitrary scaling of the spectrum SNR during testing
            if verbose:
                print('SNR adjustment factor:', runtime_params.SNRfactor)
            if runtime_params.SNRfactor:
                uncertainties /= runtime_params.SNRfactor

            # ________LOOP OVER ALL SELECTED MODELS_______
            opticalDepthProfiles = {}
            for atmosModel in atmosModels:
                if verbose:
                    print()
                    print('starting Atmospheric Model:', atmosModel)
                useTEA = bool('TEA' in atmosModel)
                if useTEA:
                    chemistry = 'TEA'
                else:
                    chemistry = 'TEC'

                # ABUNDANCES
                if 'lowmmw' in atmosModel:
                    # print(' - using a low mmw')
                    model_params['metallicity'] = 0.0  # dex
                    model_params['C/O'] = 0.0  # [C/O] (relative to solar)
                    model_params['N/O'] = 0.0  # [N/O] (relative to solar)
                else:
                    model_params['metallicity*'] = metallicity_star_dex
                    # model_params['metallicity'] = metallicity_star_dex + metallicity_planet_dex
                    # stellar metallicity has already been added (in metallicity.py)
                    model_params['metallicity'] = metallicity_planet_dex

                    # model_params['metallicity'] = 0.0
                    # print('model_params metallicity SET BY HAND!!',
                    #    model_params['metallicity'])

                    # planet C/O ratio is assumed to be solar
                    #  (0.54951 is the default in ACEChemistry, so it actually has no effect)
                    # actually, let's consider a distribution of C/O, as done for FINESSE
                    model_params['C/O'] = np.log10(
                        CtoO_planet_linear / solarCtoO
                    )
                    # print('C/O model param',model_params['C/O'])
                    if runtime_params.CtoOdaSilva:
                        model_params['N/O'] = ancil_params['NO*']
                    else:
                        model_params['N/O'] = 0

                # check whether this planet+metallicity combo is convergent/bound atmosphere
                _, Hs = calc_mmw_Hs(
                    pressure,
                    eqtemp,
                    model_params['logg'],
                    X2Hr=model_params['metallicity'],
                    useTEA=useTEA,
                )
                HoverRp = Hs / (model_params['Rp'] * sscmks['Rjup'])
                # skip non-converging atmospheres!!
                # this limit corresponds to the maximum scale height
                #   allowing for ~1e11 pressure range
                # (our actual range is typically 10 bar to 1 microbar,
                #   which requires 0.06 max H/R)
                if HoverRp > 0.04:
                    log.warning(
                        '--< WARNING UNBOUND ATMOS: %s %s ; scale height / planet radius = %s %s >--',
                        target,
                        planet_letter,
                        HoverRp,
                        atmosModel,
                    )

                fluxDepth_by_molecule = {}

                if 'cerberus' in atmosModel:
                    # CLOUD PARAMETERS
                    if 'Noclouds' in atmosModel:
                        model_params['CTP'] = (
                            3.0  # cloud deck is very deep - 1000 bars
                        )
                        model_params['HScale'] = (
                            -10.0
                        )  # small number means essentially no haze
                        model_params['HLoc'] = 0.0
                        model_params['HThick'] = 0.0
                    else:
                        # CLOUD TOP PRESSURE
                        # model_params['CTP'] = -2.  # 0.01 bar = 10 mbar
                        # HAZE SCAT CROSS SECTION SCALE FACTOR
                        # model_params['HScale'] = 0.  # some haze (1x times the nominal Jupiter model)
                        # model_params['HLoc'] = 0.
                        # model_params['HThick'] = 0.

                        # use the median from Estrella 2022 or a random selection?
                        if runtime_params.randomCloudProperties:
                            # make sure that random cloud properties are fixed between runs
                            np.random.seed(intFromTarget + 123456)
                            cloudParams = randomCloudParameters()
                            for paramName in [
                                'CTP',
                                'HScale',
                                'HLoc',
                                'HThick',
                            ]:
                                model_params[paramName] = cloudParams[paramName]
                        else:
                            # median values from Estrela et al 2022, Table 2 (TEC column)
                            cloudParams = fixedCloudParameters()
                            for paramName in [
                                'CTP',
                                'HScale',
                                'HLoc',
                                'HThick',
                            ]:
                                model_params[paramName] = cloudParams[paramName]

                    if 'R' in atmosModel:
                        iwave1 = atmosModel.index('R')
                        # print('atmosModel',atmosModel)
                        Nhires = atmosModel[iwave1 + 1 :]
                        # Nhires = int(atmosModel[iwave1+1:])
                        # print('Nhires',Nhires)
                        wavelength_um = np.logspace(-0.5, 1.0, int(Nhires))
                    else:
                        wavelength_um = ariel_instrument['wavelength']

                    if not xslib['STATUS'][-1]:
                        tempspc = {
                            'data': {planet_letter: {'WB': wavelength_um}}
                        }
                        if verbose:
                            print('CALCulating cross-sections START')

                        # Armen - you can use this to save a little time debugging maybe
                        # import pickle
                        # if 0:
                        _ = myxsecs(tempspc, runtime_params, xslib)
                        #    file = open('xslibsave.pkl', 'bw')
                        #    pickle.dump(xslib, file)
                        #    file.close()
                        # else:
                        #    file = open('xslibsave.pkl', 'br')
                        #    xslib = pickle.load(file)

                        if verbose:
                            print('CALCulating cross-sections DONE')
                    else:
                        # make sure that it exists for this planet letter
                        if planet_letter in xslib['data']:
                            # print('XSLIB ALREADY CALCULATED',planet_letter)
                            pass
                        else:
                            # print('XSLIB TRANSFERRED TO NEW PLANETLETTER')
                            existingPlanetLetter = list(xslib['data'].keys())[
                                -1
                            ]
                            # print('existingPlanetLetter',existingPlanetLetter)
                            xslib['data'][planet_letter] = xslib['data'][
                                existingPlanetLetter
                            ]
                    (
                        fluxDepth,
                        fluxDepth_by_molecule,
                        pressures,
                        opticalDepthProfiles,
                        _,
                        # moleculeProfiles,
                    ) = make_cerberus_atmos(
                        runtime_params,
                        wavelength_um,
                        model_params,
                        xslib,
                        planet_letter,
                        chemistry=chemistry,
                    )
                    # pressures should be the same thing as pressure
                    if np.any(pressures != pressure):
                        print('ERROR: pressure grids dont match!')

                elif 'taurex' in atmosModel:
                    log.error('ERROR: taurex no longer an option')
                else:
                    log.error('ERROR: unknown model')

                # bads = np.where(np.isnan(fluxDepth))
                # if len(bads[0]) > 0:
                #    print(
                #        '   bad fluxDepths:',
                #        len(bads[0]),
                #        'out of',
                #        len(fluxDepth),
                #    )

                # REBIN to the Ariel spectral resolution
                wavelength_um_rebin = []
                fluxDepth_rebin = []
                fluxDepth_by_molecule_rebin = {}
                molecules = fluxDepth_by_molecule.keys()
                for molecule in molecules:
                    fluxDepth_by_molecule_rebin[molecule] = []
                for wavelow, wavehigh in zip(
                    ariel_instrument['wavelow'],
                    ariel_instrument['wavehigh'],
                ):
                    thisWaveBin = np.where(
                        (wavelength_um >= wavelow) & (wavelength_um <= wavehigh)
                    )
                    fluxDepth_rebin.append(np.average(fluxDepth[thisWaveBin]))
                    wavelength_um_rebin.append(
                        np.average(wavelength_um[thisWaveBin])
                    )
                    for molecule in molecules:
                        # print('molecule:',molecule)
                        fluxDepth_by_molecule_rebin[molecule].append(
                            np.average(
                                fluxDepth_by_molecule[molecule][thisWaveBin]
                            )
                        )

                wavelength_um_rebin = np.array(wavelength_um_rebin)
                fluxDepth_rebin = np.array(fluxDepth_rebin)

                for molecule in molecules:
                    fluxDepth_by_molecule_rebin[molecule] = np.array(
                        fluxDepth_by_molecule_rebin[molecule]
                    )
                # print ('Nan check on raw',np.isnan(np.max(fluxDepth)))
                # print ('Nan check on rebin',np.isnan(np.max(fluxDepth_rebin)))

                # ADD OBSERVATIONAL NOISE TO THE TRUE SPECTRUM
                # set the random seed so instrument noise is fixed for each planet
                np.random.seed(intFromTarget + 1234567)
                fluxDepth_observed = fluxDepth_rebin + np.random.normal(
                    scale=uncertainties
                )

                # SAVE THE RESULTS
                # if planet_letter not in out['data'].keys():
                #    out['data'][planet_letter] = {}

                # careful - ES and ESerr are supposed to be radii, not transit dept
                #  need to take a sqrt of them
                # careful2 - watch out for sqrt of negative numbers
                signedSqrt = np.sign(fluxDepth_observed) * np.sqrt(
                    np.abs(fluxDepth_observed)
                )
                # move WB location down a level (it is independent of atmosModel)
                out['data'][planet_letter]['WB'] = wavelength_um_rebin
                out['data'][planet_letter][atmosModel] = {
                    # 'WB':wavelength_um_rebin,
                    'ES': signedSqrt,
                    'ESerr': 0.5 * uncertainties / signedSqrt,
                }

                # cerberus also wants the scale height, to normalize the spectrum
                #  keep Hs as it's own param (not inside of system_ or model_param
                #   so that it matches the formatting for regular spectra itk)
                #
                # rescale the scale height, for scaling plot axes
                Hs_forplot = Hs / (model_params['R*'] * sscmks['Rsun'])
                out['data'][planet_letter][atmosModel]['Hs'] = Hs_forplot

                # save the true spectrum (both raw and binned)
                out['data'][planet_letter][atmosModel]['true_spectrum'] = {
                    'fluxDepth': fluxDepth_rebin,
                    'wavelength_um': wavelength_um_rebin,
                    'fluxDepth_norebin': fluxDepth,
                    'wavelength_norebin': wavelength_um,
                }

                # also save the Tier level and the number of visits; add these to the plot
                out['data'][planet_letter][atmosModel]['tier'] = tier
                out['data'][planet_letter][atmosModel]['visits'] = visits

                # save the parameters used to create the spectrum. some could be useful
                # should save both observed value and truth with scatter added in
                #  'system_params' = the info in system() task
                #  'model_params' = what is actually used to create forward model
                out['data'][planet_letter][atmosModel]['system_params'] = {
                    'R*': system_params['R*'],
                    'T*': system_params['T*'],
                    'Rp': system_params[planet_letter]['rp'],
                    'Teq': system_params[planet_letter]['teq'],
                    'Mp': system_params[planet_letter]['mass'],
                }

                out['data'][planet_letter][atmosModel][
                    'model_params'
                ] = model_params.copy()
                # if testTarget:
                #    print('model_params in ariel:', model_params)

                # convert to percentage depth (just for plotting, not for the saved results)
                fluxDepth = 100 * fluxDepth
                fluxDepth_rebin = 100 * fluxDepth_rebin
                fluxDepth_observed = 100 * fluxDepth_observed
                # careful - uncertainties are reused, so don't change them permanently
                uncertainties_percent = 100 * uncertainties
                molecules = list(fluxDepth_by_molecule_rebin.keys())
                for molecule in molecules:
                    fluxDepth_by_molecule_rebin[molecule] = (
                        fluxDepth_by_molecule_rebin[molecule] * 100
                    )
                out['data'][planet_letter][atmosModel]['plot_simspectrum'] = (
                    plot_spectrum(
                        target,
                        planet_letter,
                        tier,
                        visits,
                        wavelength_um_rebin,
                        fluxDepth_rebin,
                        fluxDepth_observed,
                        uncertainties_percent,
                        molecules,
                        fluxDepth_by_molecule_rebin,
                        out['data'][planet_letter][atmosModel]['Hs'],
                        verbose=verbose,
                    )
                )
                out['data'][planet_letter][atmosModel][
                    'plot_simspectrum_topmolecules'
                ] = plot_spectrum_topmolecules(
                    target,
                    planet_letter,
                    tier,
                    visits,
                    wavelength_um_rebin,
                    ariel_instrument['wavelow'],
                    ariel_instrument['wavehigh'],
                    fluxDepth_rebin,
                    fluxDepth_observed,
                    uncertainties_percent,
                    molecules,
                    fluxDepth_by_molecule_rebin,
                    out['data'][planet_letter][atmosModel]['Hs'],
                    verbose=verbose,
                )
                out['data'][planet_letter][atmosModel]['plot_depthprobed'] = (
                    plot_depthprobed(
                        target,
                        planet_letter,
                        model_params,
                        wavelength_um_rebin,
                        pressure,
                        opticalDepthProfiles,
                        verbose=verbose,
                    )
                )

                #  ***** Armen will make a plot showing this parameter: moleculeProfiles ********

                completed_at_least_one_planet = True

    # print('completed_at_least_one_planet',completed_at_least_one_planet)
    if completed_at_least_one_planet:
        out['STATUS'].append(True)

    return out['STATUS'][-1]


# ------------------------- ------------------------------------------
