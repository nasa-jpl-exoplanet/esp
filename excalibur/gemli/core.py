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
import excalibur.util.cerberus as crbutil
from excalibur.target.targetlists import get_target_lists
from excalibur.cerberus.core import hazelib
from excalibur.cerberus.forward_model import crbFM
from excalibur.cerberus.plotters import (
    rebin_data,
    plot_corner,
    plot_spectrumfit,
)
from excalibur.gemli.plotters import (
    plot_ML_fits_vs_truths,
    plot_ML_spectrumfit,
    plot_overallsample_fits_vs_truths,
)

import logging
import os
import numpy as np

log = logging.getLogger(__name__)


# ---------------------------------- ---------------------------------
def mlfitversion():
    '''
    V1.0.0:
    '''
    return dawgie.VERSION(1, 0, 0)


# ---------------------------------------------------------------------
def features_from_one_spectrum(fluxDepth, Rs, Mp):
    '''
    fluxDepth: 1D array-like, same length/order as training spectra
    Rs: stellar radius, in Rsun
    Mp: planet mass in Mjup
    returns: (1, n_features) array ready for scaler.transform()
    '''
    x = np.asarray(fluxDepth, dtype=float)
    x_mean = x.mean()
    x_std = x.std()

    # per-spectrum normalization (row-wise)
    x_norm = (x - x_mean) / x_std

    Rp_proxy = Rs * np.sqrt(x_mean)

    # stack features
    X_features = np.concatenate([x_norm, [Rp_proxy, x_std, Rs, Mp]]).reshape(
        1, -1
    )

    return X_features


# ---------------------------------------------------------------------
def mlfit(
    trgt,
    filt,
    runtime_params,
    sysfin,
    ancillary,
    cerbxsl,
    cerbatmos,
    arielsim,
    out,
    only_these_planets=None,
    verbose=False,
):
    '''
    Plot out the results from gemli.atmos()
    trgt [INPUT]: target name
    filt [INPUT]: filter
    sysfin [INPUT]: system.finalize.parameters
    ancillary [INPUT]: ancillary.finalize.parameters
    cerbxls [INPUT]: cerberus.xslib.data
    cerbatm [INPUT]: cerberus.atmos.data
    arielsim [INPUT]: ariel.simspectrum.parameters
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    ssc = syscore.ssconstants(mks=True)

    crbhzlib = {'PROFILE': []}
    hazedir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/HAZE')
    hazelib(crbhzlib, hazedir=hazedir, verbose=False)

    if verbose:
        print('starting MLfit for target:', trgt, filt)

    completed_at_least_one_planet = False

    for p in sysfin['priors']['planets']:

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
            #    print('  ', sysfin['priors'].keys())
            #    print('  ', sysfin['priors'][p].keys())
            #    print('  Rp', sysfin['priors'][p]['rp'])
            #    print('  Mp', sysfin['priors'][p]['mass'])
            #    print('  Rs', sysfin['priors']['R*'])
            #    print('  ', arielsim['data'][p].keys())
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
                'R$_p$ (R$_\\oplus$)',
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
            MLtestsample_spectra = data['X_test_2']
            MLtestSample_input_params = data['y_test_2']

            # load scaler and use it to normalize the data
            scaler = joblib.load(os.path.join(ML_inputdir, 'scaler.pkl'))
            MLtestsample_spectra_norm = scaler.transform(MLtestsample_spectra)

            # predict results for the test sets of parameters
            MLfit_results = np.column_stack(
                [
                    model.predict(MLtestsample_spectra_norm)
                    for model in ML_models
                ]
            )

            arielModel = 'cerberus'
            # CAREFUL!  need to use a cloud-free model for self-consistency!!!
            arielModel = 'cerberusNoclouds'
            # arielModel = 'cerberusTEANoClouds'

            # MORE CAREFUL! cerberus.atmos is currently the cloudy model
            # decide here whether to use the cerb.atmos atm or arielsim spc
            # useAtmos = False
            # (currently just using arielsim spc, not cerbatmos

            realspectrum = True
            # realspectrum = False
            if realspectrum:
                Rs = sysfin['priors']['R*']
                Mp = sysfin['priors'][p]['mass']

                # oops actually which spectrum to use?
                # maybe just use the cerberus one?
                # but we want to do ML without cerberus, right? so use ariel-sim

                # print(' uppkeys', arielsim['data'][p][arielModel].keys())
                # print(' truekeys', arielsim['data'][p][arielModel]['true_spectrum'].keys())
                true_spectrum = arielsim['data'][p][arielModel]['true_spectrum']
                # observed_spectrum = true_spectrum['fluxDepth']
                # print('spec,err', np.mean(observed_spectrum), 0)
                # print(' # of waves', len(observed_spectrum))

                observed_spectrum = arielsim['data'][p][arielModel]['ES'] ** 2
                #  not used:
                # observed_spectrum_error = (
                #    2
                #    * arielsim['data'][p][arielModel]['ES']
                #    * arielsim['data'][p][arielModel]['ESerr']
                # )
                # print('dict check', arielsim['data'][p][arielModel].keys())
                # print(
                #    'spec,err',
                #    np.mean(observed_spectrum),
                #    np.mean(observed_spectrum_error),
                # )
                # another_spectrum = cerbatmos[p]['SPECTRUM'] ** 2
                # another_spectrum_error = (
                #    2 * cerbatmos[p]['SPECTRUM'] * cerbatmos[p]['ERRORS']
                # )
                # print(
                #    'spec,err',
                #    np.mean(another_spectrum),
                #    np.mean(another_spectrum_error),
                # )
                # print(' # of waves', len(another_spectrum))
            else:
                # instead of input ariel-sim info, try one of the test spectra
                jtest = 1234

                # number of wavelength bins
                # numwaves = 52
                numwaves = len(MLtestsample_spectra[jtest]) - 4
                # print('raw # of wavelengths', len(exemple))  # 56?!
                observed_spectrum = MLtestsample_spectra[jtest][:numwaves]
                # print('spectrumexample', MLtestsample_spectrum)
                # Rs = MLtestsample_spectra[jtest][numwaves + 2]
                # Mp = MLtestsample_spectra[jtest][numwaves + 3]
                Rs = MLtestsample_spectra[jtest][-2]
                Mp = MLtestsample_spectra[jtest][-1]
            # if verbose:
            #    print(' trying Rs,Mp = ', Rs, Mp)

            # pred = predict_params_from_spectrum(
            #     fluxDepth_ex, Rs, Mp, models, scaler, ML_param_names
            # )
            observed_spectrum_features = features_from_one_spectrum(
                observed_spectrum, Rs, Mp
            )
            observed_spectrum_scaled = scaler.transform(
                observed_spectrum_features
            )
            ML_param_results = [
                model.predict(observed_spectrum_scaled)[0]
                for model in ML_models
            ]
            ML_param_results = dict(zip(ML_param_names, ML_param_results))
            # if verbose:
            #    for k in ML_param_names:
            #        print(
            #            f'MLfit result for {k:7s}: {ML_param_results[k]: .6g}'
            #        )

            # calculate the systematic/ML uncertainty from the test sample
            ML_uncertainties_systematic = {}
            for k, ML_param_name in enumerate(ML_param_names):
                if ML_param_name.startswith('mlp') or (ML_param_name == 'Rp'):
                    # lump together test samples within 0.1 dex
                    thisbin = np.where(
                        np.abs(
                            MLfit_results[:, k]
                            - ML_param_results[ML_param_name]
                        )
                        < 0.1
                    )
                elif ML_param_name == 'Teq':
                    # lump together test samples within 5%
                    thisbin = np.where(
                        np.abs(
                            np.log10(
                                MLfit_results[:, k]
                                / ML_param_results[ML_param_name]
                            )
                        )
                        < np.log(1.05)
                    )
                else:
                    thisbin = [[]]
                    log.warning('Unknown ML parameter!?  %s', ML_param_name)
                # Nsample = len(MLtestSample_input_params)
                # print('binsize', ML_param_name, len(thisbin[0]) / Nsample)
                # determine the range on input values in this bin
                # (the range of input values that can produce the fit result)
                inputvals = MLtestSample_input_params[:, k][thisbin]
                # print(' ML uncertainty', ML_param_name, np.std(inputvals))
                ML_uncertainties_systematic[ML_param_name] = np.std(inputvals)

            # calculate the effect of instrument noise on ML retrieval
            # loop over a bunch of random noise selections on each data point
            # and then see how much the ML param results vary
            # (and then compare this uncertainty against the intrinsic ML uncert)

            # TEST TEST TEST
            #  increase the metallicity by some factor
            # for param in ML_param_names:
            #    if param.startswith('mlp'):
            #        ML_param_results[param] += 1.5
            # TEST TEST TEST

            # calculate the spectrum for the resulting best-fit mixratios
            ML_mixratio = {}
            for param in ML_param_names:
                if param.startswith('mlp'):
                    ML_mixratio[param[3:]] = ML_param_results[param]
            # print('ML_mixratio for forwardmodel call', ML_mixratio)

            true_spectrum = arielsim['data'][p][arielModel]['true_spectrum']
            # print('true keys', true_spectrum.keys())
            truth_spectrum = {
                'depth': true_spectrum['fluxDepth'],
                'wavelength': true_spectrum['wavelength_um'],
            }
            # print('obs keys', arielsim['data'][p][arielModel].keys())
            transitdata = {}
            transitdata['wavelength'] = arielsim['data'][p]['WB']
            transitdata['depth'] = arielsim['data'][p][arielModel]['ES'] ** 2
            transitdata['error'] = (
                2
                * arielsim['data'][p][arielModel]['ES']
                * arielsim['data'][p][arielModel]['ESerr']
            )
            transitdata = rebin_data(transitdata)

            nrandomSample = 1000  # this is fast, so can do a lot itk
            # nrandomSample = 11
            np.random.seed(123)  # make the results reproduciblea
            instrumentNoise = transitdata['error'] * np.random.normal(
                size=(nrandomSample, len(transitdata['wavelength']))
            )
            # print('instrumentNoise shape',instrumentNoise.shape)
            # then add that to the observed spectrum
            #  it seems to do the broadcasting automatically. nice
            noisedspectra = transitdata['depth'] + instrumentNoise
            # print('transit data shape', noisedspectra.shape)
            # print('transit data shape', noisedspectra[:,][0].shape)
            # print('transit data shape', noisedspectra[,:].shape)
            # print(' NOW loop over each of these 1000...')

            ML_param_results_noisedlist = []
            # for noisedspectrum in noisedspectra[:,]:
            for inoise in range(nrandomSample):
                noisedspectrum = noisedspectra[inoise, :]
                observed_spectrum_features = features_from_one_spectrum(
                    noisedspectrum, Rs, Mp
                )
                observed_spectrum_scaled = scaler.transform(
                    observed_spectrum_features
                )
                ML_param_results_noised = [
                    model.predict(observed_spectrum_scaled)[0]
                    for model in ML_models
                ]
                # print('added in', dict(zip(ML_param_names, ML_param_results)))
                ML_param_results_noisedlist.append(
                    dict(zip(ML_param_names, ML_param_results_noised))
                )
            # now check the range of ML results for each parameter
            # print('len check',len(ML_param_results_noised))
            # if 1:
            ML_uncertainties_instrument = {}
            for param in ML_param_names:
                # for instance in ML_param_results_noisedlist:
                #    print('check', instance)
                # print('param', param)
                values = [
                    instance[param] for instance in ML_param_results_noisedlist
                ]
                ML_uncertainties_instrument[param] = np.std(values)

            # print('uncertainty from the instrument', ML_uncertainties_instrument)
            # exit('asdf')

            # calculate the spectrum based on the best-fit ML model
            #  the best-fit model provided Tp, Rp, & mixing ratios
            #  set clouds/haze to zero
            # Rp = 10.**ML_param_results['Rp'] * ssc['Rearth']
            # print('Rp (RJup)', Rp)
            Rp = ML_param_results['Rp'] * ssc['Rjup']
            # print('Rp (cgs) for ML', Rp)
            # Rp = sysfin['priors'][p]['rp'] * ssc['Rjup']
            # print('IGNORING ML Rp; using system.finalize prior', Rp)

            fmc = crbFM().crbmodel(
                ML_param_results['Teq'],
                10.0,
                hazescale=1.0e-10,
                hazeloc=1.0,
                hazethick=1.0,
                mixratio=ML_mixratio,
                cheq=None,
                rp0=Rp,
                xsecs=cerbxsl[p]['XSECS'],
                qtgrid=cerbxsl[p]['QTGRID'],
                wgrid=transitdata['wavelength'],
                orbp=sysfin['priors'],
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
            simulated_spectrum = fmc.spectrum

            # add offset to match data (i.e. modify Rp)
            okPart = np.where(np.isfinite(transitdata['depth']))
            # print(
            #    'average depth before normalizing',
            #    np.mean(simulated_spectrum[okPart]),
            # )
            # OLD normalization method (+- offset)
            ML_best_fit = simulated_spectrum[okPart] + np.average(
                (transitdata['depth'][okPart] - simulated_spectrum[okPart]),
                weights=1 / transitdata['error'][okPart] ** 2,
            )
            # NEW normalization method (*/ factor)
            # ML_best_fit = simulated_spectrum[okPart] * np.average(
            #    (transitdata['depth'][okPart] / simulated_spectrum[okPart]),
            #    weights=1 / transitdata['error'][okPart] ** 2,
            # )
            # print('average depth after  normalizing', np.average(ML_best_fit))

            # print('arielsim keys', arielsim['data'][p][arielModel].keys())
            # print('cerb keys', cerbatmos[p].keys())
            # print('arielsim stuff', arielsim['data'][p][arielModel]['model_params'])
            # print('cerb stuff', cerbatmos[p]['MODELPARNAMES']) TEC and TEA; XtoH,CtoO
            # print('cerb stuff', cerbatmos[p]['TRUTH_MODELPARAMS'])
            # print('cerb stuff', cerbatmos[p]['TEC'].keys()) prior_ranges MCTRACE

            # get the true mixing ratio values, for comparison
            pgrid = np.exp(
                np.arange(
                    np.log(runtime_params.solrad) - runtime_params.Hsmax,
                    np.log(runtime_params.solrad)
                    + runtime_params.Hsmax / runtime_params.nlevels,
                    runtime_params.Hsmax / (runtime_params.nlevels - 1),
                )
            )
            pressure = pgrid[::-1]
            # ok I confirm that these are the same thing
            model_params = arielsim['data'][p][arielModel]['model_params']
            model_params = cerbatmos[p]['TRUTH_MODELPARAMS']
            # print('model_params', model_params)
            if 'TEA' in arielModel:
                tempCoeffs = [0, model_params['Teq'], 0, 1, 0, -1, 1, 0, -1, 1]
                mixratioprofiles = crbutil.calcTEA(
                    tempCoeffs,
                    pressure,
                    metallicity=10.0 ** model_params['metallicity'],
                    C_O=0.55 * 10.0 ** model_params['C/O'],
                )
                truth_params = {}
                for molecule in mixratioprofiles:
                    truth_params[molecule] = np.log10(
                        np.mean(10.0 ** mixratioprofiles[molecule])
                    )
            else:
                truth_params, _, _, _ = crbutil.crbce(
                    pressure,
                    model_params['Teq'],
                    X2Hr=model_params['metallicity'],
                    C2Or=model_params['C/O'],
                )
            # include Teq and Rp truths in with the mixing ratios
            truth_params['Teq'] = model_params['Teq']
            truth_params['Rp'] = model_params['Rp']

            out['data'][p]['ML_param_names'] = ML_param_names
            out['data'][p]['ML_param_results'] = ML_param_results
            out['data'][p][
                'ML_uncertainties_systematic'
            ] = ML_uncertainties_systematic
            out['data'][p][
                'ML_uncertainties_instrument'
            ] = ML_uncertainties_instrument
            out['data'][p]['ML_spectrum_bestfit'] = ML_best_fit
            out['data'][p]['truth_params'] = truth_params

            # plot the best-fit-by-ML model vs the data / truth
            out['data'][p]['plot_MLspectrum'], _ = plot_ML_spectrumfit(
                transitdata,
                truth_spectrum,
                truth_params,
                ML_best_fit,
                ML_param_names,
                ML_param_names_forprint,
                ML_param_results,
                ML_uncertainties_systematic,
                ML_uncertainties_instrument,
                sysfin['priors'],
                ancillary['data'][p],
                filt,
                trgt,
                p,
                verbose=verbose,
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

                rp0 = sysfin['priors'][p]['rp'] * ssc['Rjup']
                # print('Rp (cgs) cerberus', rp0)

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
                    orbp=sysfin['priors'],
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
                        orbp=sysfin['priors'],
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
                        sysfin['priors'],
                        ancillary['data'][p],
                        cerbatmos[p],
                        filt,
                        model_name,
                        trgt,
                        p,
                        saveDir=save_dir,
                    )
                )

                # if verbose:
                #    print('paramValues median  ', param_values_median)
                #    print('paramValues bestFit ', param_values_best_fit)

                # _______________ML fit-vs-truth PLOT________________
                out['data'][p]['plot_MLfitvstruth'], _ = plot_ML_fits_vs_truths(
                    ML_param_names_forprint,
                    MLtestSample_input_params,
                    MLfit_results,
                    # verbose=False,
                    verbose=verbose,  # asdf
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
                    verbose=False,
                    # verbose=verbose,  # asdf
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

    svname = 'gemli.mlfit'

    alltargetlists = get_target_lists()

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
        if verbose:
            print('  running targetlist=', targetlist['targetlistname'])
        MLparams = []
        MLresults = []
        MLtruths = []
        MLerrors = []
        MLerrorssys = []

        for trgt in targetlist['targets']:
            if verbose:
                print('        cycling through targets', trgt)
            if trgt not in aspecttargets:
                log.warning(
                    '--< GEMLI ANALYSIS: TARGET NOT IN ASPECT %s %s >--',
                    filt,
                    trgt,
                )
            elif svname + '.' + filt not in aspects[trgt]:
                # some targets don't have this filter; no problem
                # log.warning('--< NO GEMLI.MLFIT for this FILTER+TARGET %s %s >--',filt,trgt)
                pass
            elif 'STATUS' not in aspects[trgt][svname + '.' + filt]:
                log.warning(
                    '--< GEMLI ANALYSIS: FORMAT ERROR - NO STATUS %s %s >--',
                    filt,
                    trgt,
                )
            else:
                print(
                    'target with valid data format for this filter:', filt, trgt
                )
                mlfit_result = aspects[trgt][svname + '.' + filt]

                # verify SV succeeded for target
                if not mlfit_result['STATUS'][-1]:
                    log.warning(
                        '--< GEMLI ANALYSIS: STATUS IS FALSE FOR GEMLI.MLFIT %s %s >--',
                        filt,
                        trgt,
                    )
                else:
                    for planet_letter in mlfit_result['data'].keys():
                        # print('   keys:',mlfit_result['data'][planet_letter].keys())
                        if (
                            analysisplanetlist
                            and trgt + ' ' + planet_letter
                            not in analysisplanetlist['planets']
                        ):
                            # print(' DROP: Ariel doesnt observe this planet',trgt+' '+planet_letter)
                            pass
                        else:
                            mlfitresult = mlfit_result['data'][planet_letter]
                            MLparams = mlfitresult['ML_param_names']
                            MLresults.append(mlfitresult['ML_param_results'])
                            MLtruths.append(mlfitresult['truth_params'])
                            MLerrors.append(
                                mlfitresult['ML_uncertainties_instrument']
                            )
                            MLerrorssys.append(
                                mlfitresult['ML_uncertainties_systematic']
                            )
        # rearrange to standard plotting format as in cerberus
        # (a dict of lists, rather than list of dicts)
        reformatMLtruths = {}
        reformatMLresults = {}
        reformatMLerrors = {}
        reformatMLerrorssys = {}
        for param in MLparams:
            reformatMLtruths[param] = np.array(
                [MLtruth[param] for MLtruth in MLtruths]
            )
            reformatMLresults[param] = np.array(
                [MLresult[param] for MLresult in MLresults]
            )
            reformatMLerrors[param] = np.array(
                [MLerror[param] for MLerror in MLerrors]
            )
            reformatMLerrorssys[param] = np.array(
                [MLerrorsys[param] for MLerrorsys in MLerrorssys]
            )

        # set path for optional saving plot to disk
        save_dir = os.path.join(excalibur.context['data_dir'], 'bryden/')

        # for simulated data, compare retrieval against the truth
        plotarray = plot_overallsample_fits_vs_truths(
            MLparams,
            reformatMLtruths,
            reformatMLresults,
            reformatMLerrors,
            reformatMLerrorssys,
            filt,
            saveDir=save_dir,
            verbose=verbose,
        )

        # Add to SV
        out['data']['params'] = dict(MLparams)
        out['data']['truths'] = dict(MLtruths)
        out['data']['values'] = dict(MLresults)
        out['data']['errors'] = dict(MLerrors)
        out['data']['errorssys'] = dict(MLerrorssys)
        for param, plot in zip(MLparams, plotarray):
            out['data']['plot_' + param] = plot

    out['data']['targetlistnames'] = [
        targetlist['targetlistname'] for targetlist in analysistargetlists
    ]

    out['STATUS'].append(True)
    return out['STATUS'][-1]
