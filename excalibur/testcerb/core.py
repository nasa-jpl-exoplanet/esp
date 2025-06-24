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

# --------------------------------------------------------------------
def analysis(aspects, filt, out, verbose=False):
    '''
    Plot out the population analysis (retrieval vs truth, mass-metallicity, etc)
    aspects: cross-target information
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    if verbose:
        print('testcerb/analysis...')

    aspecttargets = []
    for a in aspects:
        aspecttargets.append(a)
    log.warning(
        '--< TESTCERB ANALYSIS: NUMBER OF TARGETS IN ASPECT %s >--',
        len(aspecttargets),
    )

    svname = 'testcerb.atmos'

    alltargetlists = get_target_lists()

    target = 'testJup'
    targetlist = []
    for i in range(5):  # asdf
        targetlist.append(f'{target}{i+1:03d}')
    print('targetlist',targetlist)

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

    analysistargetlists = [
        {
            'targetlistname': 'testJup',
            'targets': targetlist,
        }
    ]
    analysisplanetlist = []

    for targetlist in analysistargetlists:
        print('  running targetlist=',targetlist['targetlistname'])
        param_names = []
        masses = []
        stellar_fehs = []
        truth_values = defaultdict(list)
        fit_values = defaultdict(list)
        fit_errors = defaultdict(list)
        fit_errors2sided = defaultdict(list)

        for trgt in targetlist['targets']:
            print('        cycling through targets',trgt)
            if trgt not in aspecttargets:
                log.warning(
                    '--< TESTCERB ANALYSIS: TARGET NOT IN ASPECT %s %s >--',
                    filt,
                    trgt,
                )
            elif svname + '.' + filt not in aspects[trgt]:
                # some targets don't have this filter; no problem
                # log.warning('--< NO CERB.ATMOS for this FILTER+TARGET %s %s >--',filt,trgt)
                pass
            elif 'STATUS' not in aspects[trgt][svname + '.' + filt]:
                log.warning(
                    '--< TESTCERB ANALYSIS: FORMAT ERROR - NO STATUS %s %s >--',
                    filt,
                    trgt,
                )
            else:
                print('target with valid data format for this filter:',filt,trgt)
                atmos_fit = aspects[trgt][svname + '.' + filt]

                # if 'stellar_params' in atmosFit['data']:  # strange. this doesn't work
                if 'stellar_params' in atmos_fit['data'].keys():
                    stellar_feh = atmos_fit['data']['stellar_params']['FEH*']
                else:
                    stellar_feh = 0
                    log.warning(
                        '--< TESTCERB ANALYSIS: no FEH* for %s >--', trgt
                    )

                # verify SV succeeded for target
                if not atmos_fit['STATUS'][-1]:
                    log.warning(
                        '--< TESTCERB ANALYSIS: STATUS IS FALSE FOR CERB.ATMOS %s %s >--',
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
                                '--< TESTCERB ANALYSIS: BIG PROBLEM theres no TEC model! %s %s >--',
                                filt,
                                trgt,
                            )
                        elif (
                            'prior_ranges'
                            not in atmos_fit['data'][planet_letter]['TEC']
                        ):
                            log.warning(
                                '--< TESTCERB ANALYSIS: SKIP (no prior info) - %s %s >--',
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

                                if key == 'TEC[0]':
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
                                    true_value = atmos_fit['data'][
                                        planet_letter
                                    ]['TRUTH_MODELPARAMS'][trueparam]
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
                save_dir,
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
                fit_values, fit_errors, prior_ranges, filt, save_dir
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
            save_dir,
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

# ---------------------------------- ---------------------------------
