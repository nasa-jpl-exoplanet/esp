'''testcerb core ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements
# pylint: disable=duplicate-code


# -- IMPORTS -- ------------------------------------------------------

import os
import numpy as np

from collections import defaultdict

import excalibur
from excalibur.cerberus.bounds import set_prior_bound

from excalibur.testcerb.plotters import plot_fits_vs_truth

from collections import namedtuple
import logging


log = logging.getLogger(__name__)

TestcerbAnalysisParams = namedtuple(
    'cerberus_analysis_params_from_runtime',
    [
        'tier',
        'boundTeq',
        'boundAbundances',
        'boundCTP',
        'boundHLoc',
        'boundHScale',
        'boundHThick',
    ],
)

# --------------------------------------------------------------------
def analysis(aspects, filt, runtime_params, out, verbose=False):
    '''
    Plot out the analysis of the overall sample of test targets
    aspects: cross-target information
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    if verbose:
        print('testcerb/analysis...')

    aspecttargets = []
    for a in aspects:
        aspecttargets.append(a)
    # print('asptargs',aspecttargets)
    # aspecttargets = ['testJup001','testJup002','testJup003']
    # print('asptargs',aspecttargets)
    log.warning(
        '--< TESTCERB ANALYSIS: NUMBER OF TARGETS IN ASPECT %s >--',
        len(aspecttargets),
    )

    svname = 'testcerb.atmos'

    targetbase = 'testJup'
    targetlist = []
    for a in aspects:
        # targetlist.append(f'{target}{i+1:03d}')
        if a.startswith(targetbase):
            targetlist.append(a)
    print('targetlist', targetlist)

    # print('runtime', runtime_params)
    # prior_ranges = set_prior_bound(eqtemp, runtime_params)
    # print('using this prior range:', prior_ranges)
    # prior_ranges = None

    analysistargetlists = [
        {
            'targetlistname': 'testJup',
            'targets': targetlist,
        }
    ]

    for targetlist in analysistargetlists:
        print('  running targetlist=', targetlist['targetlistname'])
        param_names = []
        masses = []
        stellar_fehs = []
        truth_values = defaultdict(list)
        fit_values = defaultdict(list)
        fit_errors = defaultdict(list)
        fit_errors2sided = defaultdict(list)
        fit_errors2sided2sigma = defaultdict(list)

        all_traces = []
        for trgt in targetlist['targets']:
            # print('        cycling through targets', trgt)
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
                print(
                    'target with valid data format for this filter:', filt, trgt
                )
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
                            # print('prior range is passed in, not calculated:', prior_ranges)

                            traces = []
                            all_keys = []
                            for key in atmos_fit['data'][planet_letter]['TEC'][
                                'MCTRACE'
                            ]:
                                traces.append(
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

                            # careful
                            # all_traces is a list of each target iteration
                            # traces is a list of each parameter (matching all_keys)
                            all_traces.append(traces)

                            for key, trace in zip(all_keys, traces):
                                if key not in param_names:
                                    param_names.append(key)
                                med = np.median(trace)
                                fit_values[key].append(med)
                                lo = np.percentile(np.array(trace), 15.9)
                                hi = np.percentile(np.array(trace), 84.1)
                                fit_errors[key].append((hi - lo) / 2)
                                fit_errors2sided[key].append(
                                    [med - lo, hi - med]
                                )
                                lo2 = np.percentile(np.array(trace), 2.3)
                                hi2 = np.percentile(np.array(trace), 97.7)
                                fit_errors2sided2sigma[key].append(
                                    [med - lo2, hi2 - med]
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

        # compare retrieval against the truth
        plotparams, plotarray = plot_fits_vs_truth(
            all_keys,
            all_traces,
            truth_values,
            fit_values,
            # fit_errors,
            fit_errors2sided,
            fit_errors2sided2sigma,
            prior_ranges,
            filt,
        )
        # print('all_keys',all_keys)
        # print('fitvalues keys',fit_values.keys(),'len check',len(plotarray))

        # out['data']['plot_fitT'] = plotarray[0]
        # out['data']['plot_fitMetal'] = plotarray[1]
        # if len(plotarray) > 2:
        #    out['data']['plot_fitCO'] = plotarray[2]
        # if len(plotarray) > 3:
        #    out['data']['plot_fitNO'] = plotarray[3]

        for param, plot in zip(plotparams, plotarray):
            if param == 'T':
                out['data']['plot_fitT'] = plot
            elif param == '[X/H]':
                out['data']['plot_fitMetal'] = plot
            elif param == '[C/O]':
                out['data']['plot_fitCO'] = plot
            elif param == '[N/O]':
                out['data']['plot_fitNO'] = plot
            else:
                print('PROBLEM: unknown parameter has been plotted', param)

        # Add to SV
        out['data']['truths'] = dict(truth_values)
        out['data']['values'] = dict(fit_values)
        out['data']['errors'] = dict(fit_errors)
        out['data']['params'] = param_names
        out['data']['targetlistnames'] = [
            targetlist['targetlistname'] for targetlist in analysistargetlists
        ]

        out['STATUS'].append(True)

    return out['STATUS'][-1]


# ---------------------------------- ---------------------------------
