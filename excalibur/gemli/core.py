'''gemli core ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements
#  more for customDist pymc method:
# pylint: disable=invalid-name,cell-var-from-loop
# pylint: disable=duplicate-code

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import joblib
from xgboost import XGBRegressor

import excalibur
import excalibur.system.core as syscore
from excalibur.target.targetlists import get_target_lists
from excalibur.cerberus.core import hazelib
from excalibur.cerberus.forward_model import crbFM
from excalibur.cerberus.plotters import (
    rebin_data,
    plot_corner,
    plot_spectrumfit,
    plot_fits_vs_truths,
    plot_fit_uncertainties,
    plot_mass_vs_metals,
)
from excalibur.gemli.plotters import (
    plot_ML_fits_vs_truths,
    plot_ML_spectrumfit,
)

import logging
import os
import numpy as np
from collections import defaultdict

# from collections import namedtuple

log = logging.getLogger(__name__)

# GemliResultsParams = namedtuple(
#    'gemli_results_params_from_runtime',
#    [
#        'nrandomwalkers',
#        'randomseed',
#        'knownspecies',
#        'cialist',
#        'xmollist',
#        'nlevels',
#        'Hsmax',
#        'solrad',
#        'cornerBins',
#        'lbroadening',
#        'lshifting',
#        'isothermal',
#    ],
# )

# GemliAnalysisParams = namedtuple(
#    'gemli_analysis_params_from_runtime',
#    [
#        'tier',
#        'boundTeq',
#        'boundAbundances',
#        'boundCTP',
#        'boundHLoc',
#        'boundHScale',
#        'boundHThick',
#    ],
# )


# ---------------------------------- ---------------------------------
def mlfitversion():
    '''
    V1.0.0:
    '''
    return dawgie.VERSION(1, 0, 0)


# ------------------------------ -------------------------------------
def mlfit(
    trgt,
    filt,
    runtime_params,
    fin,
    anc,
    xsl,
    atm,
    spc,
    out,
    only_these_planets=None,
    verbose=False,
):
    '''
    Plot out the results from gemli.atmos()
    trgt [INPUT]: target name
    filt [INPUT]: filter
    fin [INPUT]: system.finalize.parameters
    anc [INPUT]: ancillary.finalize.parameters
    xls [INPUT]: cerberus.xslib.data
    atm [INPUT]: cerberus.atmos.data
    spectrum [INPUT]: ariel.simspectrum.parameters
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    ssc = syscore.ssconstants(mks=True)

    cerbatmos = atm['data']
    cerbxsl = xsl['data']

    crbhzlib = {'PROFILE': []}
    hazedir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/HAZE')
    hazelib(crbhzlib, hazedir=hazedir, verbose=False)

    if verbose:
        print('starting MLfit for target:', trgt, filt)

    completed_at_least_one_planet = False

    for p in fin['priors']['planets']:

        # TEC,TEA params - X/H, C/O, N/O
        # disEq params - HCN, CH4, C2H2, CO2, H2CO

        # check whether this planet was analyzed
        # (some planets are skipped, because they have an unbound atmosphere)
        if only_these_planets and p not in only_these_planets:
            log.info(
                '--< GEMLI.MLFIT: skipping non-tier2 planet %s %s >--',
                trgt,
                p,
            )
        elif p not in cerbatmos.keys():
            log.warning(
                '>-- GEMLI.MLFIT: this planet is missing cerb fit: %s %s',
                trgt,
                p,
            )

        else:
            out['data'][p] = {}

            # **** NEW CODE START HERE ****
            # print('START MLFit for planet:', p)

            # if verbose:
            #    print('Ariel-sim INPUT PARAMETERS')
            #    print('  ', fin['priors'].keys())
            #    print('  ', fin['priors'][p].keys())
            #    print('  Rp', fin['priors'][p]['rp'])
            #    print('  Mp', fin['priors'][p]['mass'])
            #    print('  Rs', fin['priors']['R*'])
            #    print('  ', spc['data'][p].keys())
            #    print()

            ML_param_names = [
                'Teq',
                'Rp',
                'mlpH2O',
                'mlpCH4',
                'mlpCO',
                'mlpCO2',
                'mlpNH3',
            ]
            ML_param_names_forprint = [
                'T$_{eq}$ (K)',
                'R$_p$ (log R$_\\oplus$)',
                'H$_2$O (log ppm)',
                'CH$_4$ (log ppm)',
                'CO (log ppm)',
                'CO$_2$ (log ppm)',
                'NH$_3$ (log ppm)',
            ]

            ML_inputdir = os.path.join(
                excalibur.context['data_dir'], 'gemli/export_model'
            )

            # load models
            ML_models = []
            for param_name in ML_param_names:
                m = XGBRegressor()
                m.load_model(os.path.join(ML_inputdir, f'{param_name}.json'))
                ML_models.append(m)

            # check how well the ML models perform over a range of test data
            # load test data
            data = np.load(os.path.join(ML_inputdir, 'test_data.npz'))
            test_spectra = data['X_test_2']
            input_params = data['y_test_2']

            # load scaler and use it to normalize the data
            scaler = joblib.load(os.path.join(ML_inputdir, 'scaler.pkl'))
            test_spectra_norm = scaler.transform(test_spectra)

            # predict results for the test sets of parameters
            MLfit_params = np.column_stack(
                [model.predict(test_spectra_norm) for model in ML_models]
            )

            def features_from_one_spectrum(fluxDepth, Rs, Mp):
                """
                fluxDepth: 1D array-like, same length/order as training spectra
                Rs: stellar radius, in Rsun
                Mp: planet mass in Mjup
                returns: (1, n_features) array ready for scaler.transform()
                """
                x = np.asarray(fluxDepth, dtype=float)
                x_mean = x.mean()
                x_std = x.std()

                # per-spectrum normalization (row-wise)
                x_norm = (x - x_mean) / x_std

                Rp_proxy = Rs * np.sqrt(x_mean)

                # stack features
                X_features = np.concatenate(
                    [x_norm, [Rp_proxy, x_std, Rs, Mp]]
                ).reshape(1, -1)

                return X_features

            # def predict_params_from_spectrum(
            #    fluxDepth, Rs, Mp, models, scaler, ML_param_names
            # ):
            #    X_features = features_from_one_spectrum(fluxDepth, Rs, Mp)
            #    X_scaled = scaler.transform(X_features)
            #    preds = [m.predict(X_scaled)[0] for m in models]
            #    return dict(zip(ML_param_names, preds))

            # asdf
            realspectrum = True
            realspectrum = False
            if realspectrum:
                Rs = fin['priors']['R*']
                Mp = fin['priors'][p]['mass']

                # oops actually which spectrum to use?
                # maybe just use the cerberus one?
                # but we want to do ML without cerberus, right? so use ariel-sim

                # (could try using 'cerberusTEA' here)
                # useArielSpectrum = False
                # print(' uppkeys', spc['data'][p]['cerberus'].keys())
                # print(' truekeys', spc['data'][p]['cerberus']['true_spectrum'].keys())
                true_spectrum = spc['data'][p]['cerberus']['true_spectrum']
                ML_spectrum = true_spectrum['fluxDepth']
                print('spec,err', np.mean(ML_spectrum), 0)
                print(' # of waves', len(ML_spectrum))

                ML_spectrum = spc['data'][p]['cerberus']['ES'] ** 2
                ML_spectrum_error = (
                    2
                    * spc['data'][p]['cerberus']['ES']
                    * spc['data'][p]['cerberus']['ESerr']
                )
                # print('dict check', spc['data'][p]['cerberus'].keys())
                print(
                    'spec,err', np.mean(ML_spectrum), np.mean(ML_spectrum_error)
                )
                another_spectrum = cerbatmos[p]['SPECTRUM'] ** 2
                another_spectrum_error = (
                    2 * cerbatmos[p]['SPECTRUM'] * cerbatmos[p]['ERRORS']
                )
                print(
                    'spec,err',
                    np.mean(another_spectrum),
                    np.mean(another_spectrum_error),
                )
                print(' # of waves', len(another_spectrum))
            else:
                # instead of input ariel-sim info, try one of the test spectra
                jtest = 123

                # number of wavelength bins
                numwaves = 52
                # print('raw # of wavelengths', len(exemple))  # 56?!
                ML_spectrum = test_spectra[jtest][:numwaves]
                # print('spectrumexample', test_spectrum)
                Rs = test_spectra[jtest][numwaves + 2]
                Mp = test_spectra[jtest][numwaves + 3]
            if verbose:
                print(' trying Rs,Mp = ', Rs, Mp)

            # pred = predict_params_from_spectrum(
            #     fluxDepth_ex, Rs, Mp, models, scaler, ML_param_names
            # )
            ML_spectrum_features = features_from_one_spectrum(
                ML_spectrum, Rs, Mp
            )
            ML_spectrum_scaled = scaler.transform(ML_spectrum_features)
            ML_param_results = [
                model.predict(ML_spectrum_scaled)[0] for model in ML_models
            ]
            ML_param_results = dict(zip(ML_param_names, ML_param_results))
            if verbose:
                for k in ML_param_names:
                    print(
                        f'MLfit result for {k:7s}: {ML_param_results[k]: .6g}'
                    )

            ML_mixratio = {}
            for param in ML_param_names:
                if param.startswith('mlp'):
                    ML_mixratio[param[3:]] = ML_param_results[param]
            print('ML_mixratio for forwardmodel call', ML_mixratio)

            true_spectrum = spc['data'][p]['cerberus']['true_spectrum']
            # print('true keys', true_spectrum.keys())
            truth_spectrum = {
                'depth': true_spectrum['fluxDepth'],
                'wavelength': true_spectrum['wavelength_um'],
            }
            # print('obs keys', spc['data'][p]['cerberus'].keys())
            transitdata = {}
            transitdata['wavelength'] = spc['data'][p]['WB']
            transitdata['depth'] = spc['data'][p]['cerberus']['ES'] ** 2
            transitdata['error'] = (
                2
                * spc['data'][p]['cerberus']['ES']
                * spc['data'][p]['cerberus']['ESerr']
            )
            transitdata = rebin_data(transitdata)

            # calculate the spectrum based on the best-fit ML model
            #  the best-fit model provided Tp, Rp, & mixing ratios
            #  set clouds/haze to zero
            fmc = crbFM().crbmodel(
                ML_param_results['Teq'],
                10.,
                hazescale=1.e-10,
                hazeloc=1.0,
                hazethick=1.0,
                mixratio=ML_mixratio,
                cheq=None,
                rp0=ML_param_results['Rp'],
                xsecs=cerbxsl[p]['XSECS'],
                qtgrid=cerbxsl[p]['QTGRID'],
                wgrid=transitdata['wavelength'],
                orbp=fin['priors'],
                hzlib=crbhzlib,
                planet=p,
                knownspecies=runtime_params.knownspecies,
                cialist=runtime_params.cialist,
                xmollist=runtime_params.xmollist,
                lbroadening=runtime_params.lbroadening,
                lshifting=runtime_params.lshifting,
                nlevels=runtime_params.nlevels,
                Hsmax=runtime_params.Hsmax,
                solrad=runtime_params.solrad,
            )
            spectrum = fmc.spectrum

            # add offset to match data (i.e. modify Rp)
            okPart = np.where(np.isfinite(transitdata['depth']))
            ML_best_fit = spectrum[okPart] + np.average(
                (transitdata['depth'][okPart] - spectrum[okPart]),
                weights=1 / transitdata['error'][okPart] ** 2,
            )
            # ML_best_fit = truth_spectrum['depth'] * 1.01

            # plot the best-fit-by-ML model vs the data / truth
            out['data'][p]['plot_MLspectrum'], _ = plot_ML_spectrumfit(
                transitdata,
                truth_spectrum,
                ML_best_fit,
                ML_param_names,
                ML_param_names_forprint,
                ML_param_results,
                fin['priors'],
                anc['data'][p],
                filt,
                trgt,
                p,
            )

            # **** NEW CODE END HERE ****

            # limit results to just the TEC model?  No, not for HST/G141
            # but do verify that TEC exists at least
            if (
                'TEC' not in cerbatmos[p]['MODELPARNAMES'].keys()
                and 'TEA' not in cerbatmos[p]['MODELPARNAMES'].keys()
            ):
                log.warning('>-- %s', 'TROUBLE: theres no TEC/TEA fit!?')
                return False

            # there was a bug before where PHOTOCHEM was passed in for Ariel
            # just in case, filter out models that are missing
            models = []
            for model_name in cerbatmos[p]['MODELPARNAMES']:
                if model_name in cerbatmos[p].keys():
                    models.append(model_name)
            for model_name in models:
                all_traces = []
                all_keys = []
                for key in cerbatmos[p][model_name]['MCTRACE']:
                    # print('going through keys in MCTRACE',key)
                    all_traces.append(cerbatmos[p][model_name]['MCTRACE'][key])
                    if model_name == 'TEC':
                        if key in ('TEC[0]', 'TEC'):
                            all_keys.append('[X/H]')
                        elif key == 'TEC[1]':
                            all_keys.append('[C/O]')
                        elif key == 'TEC[2]':
                            all_keys.append('[N/O]')
                        else:
                            all_keys.append(key)
                    elif model_name == 'TEA':
                        if key in ('TEA[0]', 'TEA'):
                            all_keys.append('[X/H]')
                        elif key == 'TEA[1]':
                            all_keys.append('[C/O]')
                        elif key == 'TEA[2]':
                            all_keys.append('[N/O]')
                        else:
                            all_keys.append(key)
                    elif model_name == 'PHOTOCHEM':
                        if key == 'PHOTOCHEM[0]':
                            all_keys.append('HCN')
                        elif key == 'PHOTOCHEM[1]':
                            all_keys.append('CH4')
                        elif key == 'PHOTOCHEM[2]':
                            all_keys.append('C2H2')
                        elif key == 'PHOTOCHEM[3]':
                            all_keys.append('CO2')
                        elif key == 'PHOTOCHEM[4]':
                            all_keys.append('H2CO')
                        else:
                            all_keys.append(key)
                    else:
                        all_keys.append(key)
                # print('all_keys',all_keys)

                # make note of the bounds placed on each parameter
                if 'prior_ranges' in cerbatmos[p][model_name].keys():
                    prior_ranges = cerbatmos[p][model_name]['prior_ranges']
                else:
                    prior_ranges = {}

                fit_cloud_parameters = 'CTP' in all_keys
                fit_n_to_o = '[N/O]' in all_keys
                fit_c_to_o = '[C/O]' in all_keys
                fit_t = 'T' in all_keys

                # save the relevant info
                transitdata = {}
                transitdata['wavelength'] = cerbatmos[p]['WAVELENGTH']
                transitdata['depth'] = cerbatmos[p]['SPECTRUM'] ** 2
                transitdata['error'] = (
                    2 * cerbatmos[p]['SPECTRUM'] * cerbatmos[p]['ERRORS']
                )

                truth_spectrum = None
                truth_params = None
                if 'sim' in filt:
                    if 'TRUTH_SPECTRUM' in cerbatmos[p].keys():
                        # print('NOT ERROR: true spectrum found in atmos output')
                        truth_spectrum = {
                            'depth': cerbatmos[p]['TRUTH_SPECTRUM'],
                            'wavelength': cerbatmos[p]['TRUTH_WAVELENGTH'],
                        }
                        truth_params = cerbatmos[p]['TRUTH_MODELPARAMS']
                    else:
                        log.error(
                            'ERROR: true spectrum is missing from the atmos output'
                        )
                elif 'TRUTH_SPECTRUM' in cerbatmos[p].keys():
                    log.error(
                        'ERROR: true spectrum is present for non-simulated data'
                    )

                if fit_t:
                    tprtrace = cerbatmos[p][model_name]['MCTRACE']['T']
                    tpr = np.median(tprtrace)
                else:
                    if ('TRUTH_MODELPARAMS' in cerbatmos[p]) and (
                        'Teq' in cerbatmos[p]['TRUTH_MODELPARAMS']
                    ):
                        # print('truth params',cerbatmos[p]['TRUTH_MODELPARAMS'])
                        tpr = cerbatmos[p]['TRUTH_MODELPARAMS']['Teq']
                    else:
                        tpr = 666

                mdplist = [
                    key
                    for key in cerbatmos[p][model_name]['MCTRACE']
                    if model_name in key
                ]
                # print('mdplist',mdplist)
                mdptrace = []
                for key in mdplist:
                    mdptrace.append(cerbatmos[p][model_name]['MCTRACE'][key])
                if fit_cloud_parameters:
                    ctptrace = cerbatmos[p][model_name]['MCTRACE']['CTP']
                    hazescaletrace = cerbatmos[p][model_name]['MCTRACE'][
                        'HScale'
                    ]
                    hazeloctrace = cerbatmos[p][model_name]['MCTRACE']['HLoc']
                    hazethicktrace = cerbatmos[p][model_name]['MCTRACE'][
                        'HThick'
                    ]
                    ctp = np.median(ctptrace)
                    hazescale = np.median(hazescaletrace)
                    hazeloc = np.median(hazeloctrace)
                    hazethick = np.median(hazethicktrace)
                    # print('fit results; CTP:', ctp)
                    # print('fit results; HScale:', hazescale)
                    # print('fit results; HLoc:', hazeloc)
                    # print('fit results; HThick:', hazethick)
                else:
                    ctp = cerbatmos[p]['TRUTH_MODELPARAMS']['CTP']
                    hazescale = cerbatmos[p]['TRUTH_MODELPARAMS']['HScale']
                    hazeloc = cerbatmos[p]['TRUTH_MODELPARAMS']['HLoc']
                    hazethick = cerbatmos[p]['TRUTH_MODELPARAMS']['HThick']
                mdp = np.median(np.array(mdptrace), axis=1)
                # print('fit results; T:',tpr)
                # print('fit results; mdplist:',mdp)

                rp0 = fin['priors'][p]['rp'] * ssc['Rjup']

                if model_name in ['TEC', 'TEA']:
                    # if len(mdp)!=3: log.warning('--< Expecting 3 molecules for TEQ model! >--')
                    mixratio = None
                    tceqdict = {}
                    tceqdict['XtoH'] = float(mdp[0])
                    if fit_c_to_o:
                        tceqdict['CtoO'] = float(mdp[1])
                    else:
                        if ('TRUTH_MODELPARAMS' in cerbatmos[p]) and (
                            'CtoO' in cerbatmos[p]['TRUTH_MODELPARAMS']
                        ):
                            # print('truth params',cerbatmos[p]['TRUTH_MODELPARAMS'])
                            tceqdict['CtoO'] = cerbatmos[p][
                                'TRUTH_MODELPARAMS'
                            ]['CtoO']
                        else:
                            # default is Solar
                            tceqdict['CtoO'] = 0.0

                    if fit_n_to_o:
                        tceqdict['NtoO'] = float(mdp[2])
                    else:
                        if ('TRUTH_MODELPARAMS' in cerbatmos[p]) and (
                            'NtoO' in cerbatmos[p]['TRUTH_MODELPARAMS']
                        ):
                            # print('truth params',cerbatmos[p]['TRUTH_MODELPARAMS'])
                            tceqdict['NtoO'] = cerbatmos[p][
                                'TRUTH_MODELPARAMS'
                            ]['NtoO']
                        else:
                            # default is Solar
                            tceqdict['NtoO'] = 0.0
                elif model_name == 'PHOTOCHEM':
                    if len(mdp) != 5:
                        log.warning(
                            '--< Expecting 5 molecules for PHOTOCHEM model! >--'
                        )
                    tceqdict = None
                    mixratio = {}
                    mixratio['HCN'] = float(mdp[0])
                    mixratio['CH4'] = float(mdp[1])
                    mixratio['C2H2'] = float(mdp[2])
                    mixratio['CO2'] = float(mdp[3])
                    mixratio['H2CO'] = float(mdp[4])

                else:
                    log.warning(
                        '--< Expecting TEQ, TEC, or PHOTOCHEM model! >--'
                    )

                param_values_median = (
                    tpr,
                    ctp,
                    hazescale,
                    hazeloc,
                    hazethick,
                    tceqdict,
                    mixratio,
                )

                # print('median fmc',np.nanmedian(fmc))
                fmc = np.zeros(transitdata['depth'].size)
                fmc = crbFM().crbmodel(
                    float(tpr),
                    float(ctp),
                    hazescale=float(hazescale),
                    hazeloc=float(hazeloc),
                    hazethick=float(hazethick),
                    mixratio=mixratio,
                    cheq=tceqdict,
                    rp0=rp0,
                    xsecs=cerbxsl[p]['XSECS'],
                    qtgrid=cerbxsl[p]['QTGRID'],
                    wgrid=transitdata['wavelength'],
                    orbp=fin['priors'],
                    hzlib=crbhzlib,
                    planet=p,
                    knownspecies=runtime_params.knownspecies,
                    cialist=runtime_params.cialist,
                    xmollist=runtime_params.xmollist,
                    lbroadening=runtime_params.lbroadening,
                    lshifting=runtime_params.lshifting,
                    nlevels=runtime_params.nlevels,
                    Hsmax=runtime_params.Hsmax,
                    solrad=runtime_params.solrad,
                )
                spectrum = fmc.spectrum

                # add offset to match data (i.e. modify Rp)
                okPart = np.where(np.isfinite(transitdata['depth']))
                patmos_model = spectrum[okPart] + np.average(
                    (transitdata['depth'][okPart] - spectrum[okPart]),
                    weights=1 / transitdata['error'][okPart] ** 2,
                )

                # calculate chi2 values to see which is the best fit
                offsets_model = (
                    patmos_model - transitdata['depth'][okPart]
                ) / transitdata['error'][okPart]
                chi2model = np.nansum(offsets_model**2)
                # print('chi2model', chi2model)

                # make an array of some randomly selected walker results
                # fix the random seed for each target/planet, so that results are reproducable
                int_from_target = (
                    1  # arbitrary initialization for the random seed
                )
                for char in trgt + ' ' + p:
                    int_from_target = (
                        runtime_params.randomseed * int_from_target + ord(char)
                    ) % 100000
                np.random.seed(int_from_target)

                chi2best = chi2model
                patmos_best_fit = patmos_model
                param_values_best_fit = param_values_median
                spectrumarray = []
                nwalkersteps = len(np.array(mdptrace)[0, :])
                # print('# of walker steps', nwalkersteps)
                for _ in range(runtime_params.nrandomwalkers):
                    iwalker = int(nwalkersteps * np.random.rand())

                    if fit_cloud_parameters:
                        ctp = ctptrace[iwalker]
                        hazescale = hazescaletrace[iwalker]
                        hazeloc = hazeloctrace[iwalker]
                        hazethick = hazethicktrace[iwalker]
                    if fit_t:
                        tpr = tprtrace[iwalker]
                    mdp = np.array(mdptrace)[:, iwalker]
                    # print('shape mdp',mdp.shape)
                    # if runtime_params.fitCloudParameters:
                    #    print('fit results; CTP:', ctp)
                    #    print('fit results; HScale:', hazescale)
                    #    print('fit results; HLoc:', hazeloc)
                    #    print('fit results; HThick:', hazethick)
                    # print('fit results; T:', tpr)
                    # print('fit results; mdplist:', mdp)

                    if model_name in ['TEC', 'TEA']:
                        mixratio = None
                        tceqdict = {}
                        tceqdict['XtoH'] = float(mdp[0])
                        if fit_c_to_o:
                            tceqdict['CtoO'] = float(mdp[1])
                        else:
                            tceqdict['CtoO'] = cerbatmos[p][
                                'TRUTH_MODELPARAMS'
                            ]['C/O']
                        if fit_n_to_o:
                            tceqdict['NtoO'] = float(mdp[2])
                        else:
                            if ('TRUTH_MODELPARAMS' in cerbatmos[p]) and (
                                'NtoO' in cerbatmos[p]['TRUTH_MODELPARAMS']
                            ):
                                tceqdict['NtoO'] = cerbatmos[p][
                                    'TRUTH_MODELPARAMS'
                                ]['NtoO']
                            else:
                                # log.info('--< NtoO is missing from TRUTH_MODELPARAMS >--')
                                tceqdict['NtoO'] = 0.0

                    elif model_name == 'PHOTOCHEM':
                        tceqdict = None
                        mixratio = {}
                        mixratio['HCN'] = float(mdp[0])
                        mixratio['CH4'] = float(mdp[1])
                        mixratio['C2H2'] = float(mdp[2])
                        mixratio['CO2'] = float(mdp[3])
                        mixratio['H2CO'] = float(mdp[4])

                    fmcrand = crbFM().crbmodel(
                        float(tpr),
                        float(ctp),
                        hazescale=float(hazescale),
                        hazeloc=float(hazeloc),
                        hazethick=float(hazethick),
                        mixratio=mixratio,
                        rp0=rp0,
                        xsecs=cerbxsl[p]['XSECS'],
                        qtgrid=cerbxsl[p]['QTGRID'],
                        wgrid=transitdata['wavelength'],
                        orbp=fin['priors'],
                        hzlib=crbhzlib,
                        cheq=tceqdict,
                        planet=p,
                        knownspecies=runtime_params.knownspecies,
                        cialist=runtime_params.cialist,
                        xmollist=runtime_params.xmollist,
                        lbroadening=runtime_params.lbroadening,
                        lshifting=runtime_params.lshifting,
                        nlevels=runtime_params.nlevels,
                        Hsmax=runtime_params.Hsmax,
                        solrad=runtime_params.solrad,
                    )
                    spectrumrand = fmcrand.spectrum
                    # add offset to match data (i.e. modify Rp)
                    okPart = np.where(np.isfinite(transitdata['depth']))
                    patmos_modelrand = spectrumrand[okPart] + np.average(
                        (transitdata['depth'][okPart] - spectrumrand[okPart]),
                        weights=1 / transitdata['error'][okPart] ** 2,
                    )
                    spectrumarray.append(patmos_modelrand)

                    # check to see if this model is the best one
                    offsets_modelrand = (
                        patmos_modelrand - transitdata['depth'][okPart]
                    ) / transitdata['error'][okPart]
                    chi2modelrand = np.nansum(offsets_modelrand**2)
                    # print('chi2 for a random walker', chi2modelrand)
                    # print('chi2modelrand', chi2modelrand)
                    # print('chi2best', chi2best)
                    if chi2modelrand < chi2best:
                        # print('  using this as best', chi2modelrand)
                        chi2best = chi2modelrand
                        patmos_best_fit = patmos_modelrand
                        param_values_best_fit = (
                            tpr,
                            ctp,
                            hazescale,
                            hazeloc,
                            hazethick,
                            tceqdict,
                            mixratio,
                        )

                # _______________MAKE SOME PLOTS________________
                save_dir = os.path.join(
                    excalibur.context['data_dir'], 'bryden/'
                )

                # _______________BEST-FIT SPECTRUM PLOT________________
                transitdata = rebin_data(transitdata)

                out['data'][p]['plot_spectrum_' + model_name], _ = (
                    plot_spectrumfit(
                        transitdata,
                        patmos_model,
                        patmos_model,
                        patmos_best_fit,
                        spectrumarray,
                        truth_spectrum,
                        fin['priors'],
                        anc['data'][p],
                        cerbatmos[p],
                        filt,
                        model_name,
                        trgt,
                        p,
                        saveDir=save_dir,
                    )
                )

                if verbose:
                    print('paramValues median  ', param_values_median)
                    print('paramValues bestFit ', param_values_best_fit)

                # _______________ML fit-vs-truth PLOT________________
                out['data'][p]['plot_MLfitvstruth'], _ = plot_ML_fits_vs_truths(
                    ML_param_names_forprint,
                    input_params,
                    MLfit_params,
                    verbose=verbose,
                )

                # _______________CORNER PLOT________________
                out['data'][p]['plot_corner_' + model_name], _ = plot_corner(
                    all_keys,
                    all_traces,
                    all_traces,
                    param_values_best_fit,
                    truth_params,
                    prior_ranges,
                    filt,
                    model_name,
                    trgt,
                    p,
                    bins=runtime_params.cornerBins,
                    verbose=verbose,
                    saveDir=save_dir,
                )
            out['target'].append(trgt)
            out['planets'].append(p)
            completed_at_least_one_planet = True
            # print('out-data keys at end of this planet',out['data'][p].keys())

    if completed_at_least_one_planet:
        out['STATUS'].append(True)
    return out['STATUS'][-1]


# ---------------------------------- ---------------------------------


def analysis(aspects, filt, runtime_params, out, verbose=False):
    '''
    Plot out the population analysis (retrieval vs truth, mass-metallicity, etc)
    aspects: cross-target information
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    if verbose:
        print('gemli/analysis...')

    aspecttargets = []
    for a in aspects:
        aspecttargets.append(a)
    log.info(
        '--< GEMLI ANALYSIS: NUMBER OF TARGETS IN ASPECT %s >--',
        len(aspecttargets),
    )

    svname = 'gemli.atmos'

    alltargetlists = get_target_lists()

    # set prior_ranges to avoid possible used-before-assignment problem
    # (ideally it is read in, but possibly not if there's mistake/old formatting)
    # the normal call doesn't work well here actually. and it creates nodes
    # darn.  have to just set something arbitrary
    # _, prior_ranges = addPriors(priorRangeTable, runtime_params, model, modparlbl[model])
    prior_ranges = None

    # allow for analysis of multiple target lists
    analysistargetlists = []
    # optionally specify the specific planets within multi-planet systems
    analysisplanetlist = []

    if filt == 'Ariel-sim':
        if runtime_params.tier == 2:
            #  *** Tier-2 (259 planets) ***
            analysistargetlists.append(
                {
                    'targetlistname': '2-year science time (Tier-2); Thorngren mmw (Nov.2024)',
                    'targets': alltargetlists['ariel_Nov2024_2years'],
                }
            )
            analysisplanetlist = {
                'planetlistname': '2-year science time (Tier-2); Thorngren mmw (Nov.2024)',
                'planets': alltargetlists[
                    'ariel_Nov2024_2years_withPlanetletters'
                ],
            }
        elif runtime_params.tier == 1:
            #  *** Tier-1 (626 planets) ***
            analysistargetlists.append(
                {
                    'targetlistname': '2-year science time (Tier-1); Thorngren mmw (Nov.2024)',
                    'targets': alltargetlists['ariel_Nov2024_2yearsTier1'],
                }
            )
            analysisplanetlist = {
                'planetlistname': '2-year science time (Tier-1); Thorngren mmw (Aug.2024)',
                'planets': alltargetlists[
                    'ariel_Nov2024_2yearsTier1_withPlanetletters'
                ],
            }
        else:
            log.error(
                'ERROR: unknown tier level for mass-metal plot %s',
                runtime_params.tier,
            )
    else:
        analysistargetlists.append(
            {
                'targetlistname': 'All Excalibur targets',
                'targets': alltargetlists['active'],
            }
        )
        # analysistargetlists.append({
        #    'targetlistname':'Roudier+ 2022',
        #    'targets':alltargetlists['roudier62']})
        # analysistargetlists.append({
        #    'targetlistname':'All G141 targets',
        #    'targets':alltargetlists['G141']})

    for targetlist in analysistargetlists:
        # print('  running targetlist=',targetlist['targetlistname'])
        param_names = []
        masses = []
        stellar_fehs = []
        truth_values = defaultdict(list)
        fit_values = defaultdict(list)
        fit_errors = defaultdict(list)
        fit_errors2sided = defaultdict(list)

        # FIXMEE: move to config file and fix this code
        for trgt in targetlist['targets']:
            # print('        cycling through targets',trgt)
            if trgt not in aspecttargets:
                log.warning(
                    '--< GEMLI ANALYSIS: TARGET NOT IN ASPECT %s %s >--',
                    filt,
                    trgt,
                )
            elif svname + '.' + filt not in aspects[trgt]:
                # some targets don't have this filter; no problem
                # log.warning('--< NO CERB.ATMOS for this FILTER+TARGET %s %s >--',filt,trgt)
                pass
            elif 'STATUS' not in aspects[trgt][svname + '.' + filt]:
                log.warning(
                    '--< GEMLI ANALYSIS: FORMAT ERROR - NO STATUS %s %s >--',
                    filt,
                    trgt,
                )
            else:
                # print('target with valid data format for this filter:',filt,trgt)
                atmos_fit = aspects[trgt][svname + '.' + filt]

                # if 'stellar_params' in atmosFit['data']:  # strange. this doesn't work
                if 'stellar_params' in atmos_fit['data'].keys():
                    stellar_feh = atmos_fit['data']['stellar_params']['FEH*']
                else:
                    stellar_feh = 0
                    log.warning('--< GEMLI ANALYSIS: no FEH* for %s >--', trgt)

                # verify SV succeeded for target
                if not atmos_fit['STATUS'][-1]:
                    log.warning(
                        '--< GEMLI ANALYSIS: STATUS IS FALSE FOR CERB.ATMOS %s %s >--',
                        filt,
                        trgt,
                    )
                else:
                    for planet_letter in atmos_fit['data'].keys():
                        # print(trgt,atmosFit['data'][planet_letter]['MODELPARNAMES'])
                        # print(trgt,atmosFit['data'][planet_letter]['planet_params'])

                        # print('   keys:',atmosFit['data'][planet_letter].keys())
                        if (
                            planet_letter == 'stellar_params'
                        ):  # this is not a planet letter
                            pass

                        elif (
                            analysisplanetlist
                            and trgt + ' ' + planet_letter
                            not in analysisplanetlist['planets']
                        ):
                            # print(' DROP: Ariel doesnt observe this planet',trgt+' '+planet_letter)
                            pass

                        elif (
                            'TEC'
                            not in atmos_fit['data'][planet_letter][
                                'MODELPARNAMES'
                            ]
                        ):
                            log.warning(
                                '--< GEMLI ANALYSIS: BIG PROBLEM theres no TEC model! %s %s >--',
                                filt,
                                trgt,
                            )
                        elif (
                            'prior_ranges'
                            not in atmos_fit['data'][planet_letter]['TEC']
                        ):
                            log.warning(
                                '--< GEMLI ANALYSIS: SKIP (no prior info) - %s %s >--',
                                filt,
                                trgt,
                            )
                        else:
                            if (
                                'planet_params'
                                in atmos_fit['data'][planet_letter]
                            ):
                                masses.append(
                                    atmos_fit['data'][planet_letter][
                                        'planet_params'
                                    ]['mass']
                                )
                            else:
                                masses.append(666)

                            stellar_fehs.append(stellar_feh)

                            # (prior range should be the same for all the targets)
                            prior_ranges = atmos_fit['data'][planet_letter][
                                'TEC'
                            ]['prior_ranges']

                            all_traces = []
                            all_keys = []
                            for key in atmos_fit['data'][planet_letter]['TEC'][
                                'MCTRACE'
                            ]:
                                all_traces.append(
                                    atmos_fit['data'][planet_letter]['TEC'][
                                        'MCTRACE'
                                    ][key]
                                )

                                if key in ('TEC[0]', 'TEC'):
                                    all_keys.append('[X/H]')
                                elif key == 'TEC[1]':
                                    all_keys.append('[C/O]')
                                elif key == 'TEC[2]':
                                    all_keys.append('[N/O]')
                                else:
                                    all_keys.append(key)

                            for key, trace in zip(all_keys, all_traces):
                                if key not in param_names:
                                    param_names.append(key)
                                med = np.median(trace)
                                fit_values[key].append(med)
                                lo = np.percentile(np.array(trace), 16)
                                hi = np.percentile(np.array(trace), 84)
                                fit_errors[key].append((hi - lo) / 2)
                                fit_errors2sided[key].append(
                                    [med - lo, hi - med]
                                )
                                if verbose:
                                    if key == '[N/O]' and (hi - lo) / 2 < 2:
                                        print(
                                            'N/O',
                                            trgt,
                                            np.median(trace),
                                            (hi - lo) / 2,
                                        )
                            if (
                                'TRUTH_MODELPARAMS'
                                in atmos_fit['data'][planet_letter].keys()
                            ) and (
                                isinstance(
                                    atmos_fit['data'][planet_letter][
                                        'TRUTH_MODELPARAMS'
                                    ],
                                    dict,
                                )
                            ):
                                truth_params = atmos_fit['data'][planet_letter][
                                    'TRUTH_MODELPARAMS'
                                ].keys()
                                # print('truth keys:',truth_params)
                            else:
                                truth_params = []

                            for trueparam, fitparam in zip(
                                ['Teq', 'metallicity', 'C/O', 'N/O', 'Mp'],
                                ['T', '[X/H]', '[C/O]', '[N/O]', 'Mp'],
                            ):
                                if trueparam in truth_params:
                                    true_value = float(
                                        atmos_fit['data'][planet_letter][
                                            'TRUTH_MODELPARAMS'
                                        ][trueparam]
                                    )
                                    # (metallicity and C/O do not have to be converted to log-solar)
                                    # if trueparam=='metallicity':
                                    #    true_value = np.log10(true_value)
                                    # elif trueparam=='C/O':
                                    #    true_value = np.log10(true_value/0.54951)  # solar is C/O=0.55
                                    # elif trueparam=='N/O':
                                    #     true_value = true_value
                                    if (
                                        fitparam == '[N/O]'
                                        and true_value == 666
                                    ):
                                        truth_values[fitparam].append(0)
                                    else:
                                        truth_values[fitparam].append(
                                            true_value
                                        )

                                    if verbose:
                                        if (
                                            trueparam == 'Teq'
                                            and true_value > 3333
                                        ):
                                            print(
                                                'strangely high T',
                                                trgt,
                                                true_value,
                                            )
                                        if (
                                            trueparam == 'metallicity'
                                            and true_value > 66
                                        ):
                                            print(
                                                'strangely high [X/H]',
                                                trgt,
                                                true_value,
                                            )
                                            print(
                                                'atmosFit',
                                                atmos_fit['data'][
                                                    planet_letter
                                                ],
                                            )
                                        if (
                                            trueparam == 'C/O'
                                            and true_value > 0.5
                                        ):
                                            print(
                                                'strangely high [C/O]',
                                                trgt,
                                                true_value,
                                            )

                                elif trueparam == 'Mp':
                                    # if the planet mass is not in the Truth dictionary, pull it from system
                                    # print(' input keys',atmosFit['data'][planet_letter]['planet_params'])
                                    # print(' planet mass from system params:',
                                    #      atmosFit['data'][planet_letter]['planet_params']['mass'])
                                    truth_values[fitparam].append(
                                        atmos_fit['data'][planet_letter][
                                            'planet_params'
                                        ]['mass']
                                    )
                                elif trueparam == 'N/O':
                                    truth_values[fitparam].append(0)
                                else:
                                    truth_values[fitparam].append(666)

        # plot analysis of the results.  save as png and as state vector for states/view
        save_dir = os.path.join(excalibur.context['data_dir'], 'bryden/')
        fit_co_plot = False
        fit_no_plot = False
        if 'sim' in filt:
            # for simulated data, compare retrieval against the truth
            #  note that the length of plotarray depends on whether N/O and C/O are fit parameters
            # jenkins doesn't like to have a triple-packed return here because it's fussy
            plotarray = plot_fits_vs_truths(
                truth_values,
                fit_values,
                fit_errors,
                prior_ranges,
                filt,
                saveDir=save_dir,
                verbose=verbose,
            )
            # fitTplot, fitMetalplot, fitCOplot, fitNOplot = plotarray[0],plotarray[1],plotarray[2],plotarray[3]
            fit_t_plot = plotarray[0]
            fit_metalplot = plotarray[1]
            if len(plotarray) > 2:
                fit_co_plot = plotarray[2]
            if len(plotarray) > 3:
                fit_no_plot = plotarray[3]
        else:
            # for real data, make a histogram of the retrieved uncertainties
            #  note that the length of plotarray depends on whether N/O and C/O are fit parameters
            plotarray = plot_fit_uncertainties(
                fit_values,
                fit_errors,
                prior_ranges,
                filt,
                saveDir=save_dir,
                verbose=verbose,
            )
            fit_t_plot = plotarray[0]
            fit_metalplot = plotarray[1]
            if len(plotarray) > 2:
                fit_co_plot = plotarray[2]
            if len(plotarray) > 3:
                fit_no_plot = plotarray[3]

        mass_metals_plot, _ = plot_mass_vs_metals(
            truth_values['Mp'],
            stellar_fehs,
            truth_values,
            fit_values,
            fit_errors2sided,
            prior_ranges,
            filt,
            saveDir=save_dir,
            verbose=verbose,
        )

        # save the analysis as .csv file? (in /proj/data/spreadsheets/)
        # savesv(aspects, targetlist)

        # targetlistname = targetlist['targetlistname']

        # Add to SV
        out['data']['truths'] = dict(truth_values)
        out['data']['values'] = dict(fit_values)
        out['data']['errors'] = dict(fit_errors)
        out['data']['plot_mass_v_metals'] = mass_metals_plot
        out['data']['plot_fitT'] = fit_t_plot
        out['data']['plot_fitMetal'] = fit_metalplot
        if fit_co_plot:
            out['data']['plot_fitCO'] = fit_co_plot
        if fit_no_plot:
            out['data']['plot_fitNO'] = fit_no_plot

        out['data']['params'] = param_names
    out['data']['targetlistnames'] = [
        targetlist['targetlistname'] for targetlist in analysistargetlists
    ]

    out['STATUS'].append(True)
    return out['STATUS'][-1]
