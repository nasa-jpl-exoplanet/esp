'''cerberus plotting ds'''

# Heritage code shame:
#  no-member is for scipy.special.erfinv
# pylint: disable=invalid-name,no-member
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import numpy as np
import scipy
import matplotlib.pyplot as plt

from excalibur.system.core import ssconstants
from excalibur.util.plotters import save_plot_tosv

import logging


log = logging.getLogger(__name__)

ssc = ssconstants(mks=True)


# --------------------------------------------------------------------
def plot_fits_vs_truth(
    all_keys,
    all_traces,
    truth_values,
    fit_values,
    fit_errors,
    fit_errors2sigma,
    prior_ranges,
    filt,
    saveDir='./',
    savetodisk=False,
):
    '''
    Compare the retrieved values against the original inputs
    '''

    paramlist = []
    for param in ['T', '[X/H]', '[C/O]', '[N/O]']:
        if (
            len(fit_values[param]) == 0
            or len(np.where(fit_values[param] != fit_values[param][0])[0]) == 0
        ):
            # print('drop a blank truth parameter',param)  # (N/O sometimes dropped)
            pass
        else:
            paramlist.append(param)

    plot_statevectors = []
    for param in paramlist:
        # print('PARAM:',param)

        figure = plt.figure(figsize=(13, 4))
        ax = figure.add_subplot(1, 3, 2)

        Niter = len(fit_values[param])

        for testIndex, (truth, fit, error, error2) in enumerate(
            zip(
                truth_values[param],
                fit_values[param],
                fit_errors[param],
                fit_errors2sigma[param],
            )
        ):
            ax.plot([0, 1000], [truth, truth], 'k--')
            ax.scatter(
                testIndex + 1,
                fit,
                facecolor='k',
                edgecolor='k',
                s=30,
                zorder=5,
            )
            # switch to 2-sided (asymmetric) error bars
            ax.errorbar(
                testIndex + 1,
                fit,
                yerr=np.array([[error[0]], [error[1]]]),
                # yerr=error,
                fmt='.',
                color='k',
                lw=2,
                zorder=4,
            )
            ax.errorbar(
                testIndex + 1,
                fit,
                yerr=np.array([[error2[0]], [error2[1]]]),
                fmt='.',
                color='grey',
                lw=0.5,
                zorder=4,
            )
        ax.set_xlim(0, Niter + 1)

        # ax.set_xlabel('test index', fontsize=14)
        ax.set_xlabel('', fontsize=14)
        ax.set_xticklabels([])
        ax.set_ylabel(param + ' fit', fontsize=14)
        # plt.text(
        #    Niter - 3,
        #    truth + 0.03,
        #    param + ' truth',
        #    c='black',
        #    fontsize=12,
        # )

        # plot C/O=1 as a dotted line
        if param == '[C/O]':
            solarCO = np.log10(0.55)
            ax.plot([-100, 100], [-solarCO, -solarCO], 'k:', lw=1, zorder=3)
            # plt.text(
            #    Niter - 3,
            #    -solarCO + 0.03,
            #    'C/O=1',
            #    c='black',
            #    fontsize=12,
            # )

        if param == 'T':  # the prior for T varies between targets
            # Teq = truth_values[param][0]
            # ax.set_ylim(Teq * prior_ranges[param][0], Teq * prior_ranges[param][1])
            # oops this prior_range has already been multiplied by Teq. don't do it again
            ax.set_ylim(prior_ranges[param][0], prior_ranges[param][1])
        elif param not in prior_ranges:
            log.warning(
                '--< Cerb.analysis: Parameter missing from prior_range 2: %s >--',
                param,
            )
        else:
            ax.set_ylim(prior_ranges[param][0], prior_ranges[param][1])

        if param == 'T':
            plt.title('temperature retrieval for ' + str(Niter) + ' planets')
        elif param == '[X/H]':
            plt.title('metallicity retrieval for ' + str(Niter) + ' planets')
        elif param == '[C/O]':
            plt.title('C/O retrieval for ' + str(Niter) + ' planets')
        else:
            plt.title(param + ' retrieval for ' + str(Niter) + ' planets')

        # CHI-OFFSET HISTOGRAMS IN SECOND PANEL
        ax2 = figure.add_subplot(1, 3, 1)

        # old method - simple gaussian-assuming sigma definition
        offsets = np.array(fit_values[param]) - np.array(truth_values[param])
        errors = np.median(np.array(fit_errors[param]), axis=1)
        chis = offsets / errors
        print('old chis', chis)
        # middle method - use two-sided uncertainty, rather than average uncertainty
        offsets = np.array(fit_values[param]) - np.array(truth_values[param])
        errorslo = np.array([errtuple[0] for errtuple in fit_errors[param]])
        errorshi = np.array([errtuple[1] for errtuple in fit_errors[param]])
        print(' lo', errorslo)
        print(' hi', errorshi)
        lowside = np.where(offsets <= 0)  # case where fit is below truth
        highside = np.where(offsets > 0)  # case where fit is above truth
        print(' low', lowside)
        print(' hi', highside)
        # print( offsets[lowside])
        # print(errorshi[lowside])
        chis[lowside] = offsets[lowside] / errorshi[lowside]
        chis[highside] = offsets[highside] / errorslo[highside]
        print('mid chis', chis)
        # print()
        # correct method - calulcate percentile and convert to equivalent sigma
        # print('trace shape',len(all_traces))
        chis = []
        # print('number of targets in all_traces',len(all_traces))
        for truth, traces in zip(
            truth_values[param], all_traces
        ):  # loop over each target
            # print(' number of parameters in this trace',len(traces))
            # print('   index',all_keys.index(param),len(all_keys))
            trace = traces[all_keys.index(param)]  # select the parameter index
            percentile = len(np.where(trace > truth)[0]) / len(trace)
            # percentile = len(trace < truth) / len(trace)  fails
            # print('percentil',percentile*100)
            if percentile <= 0:
                chis.append(-123)  # truth is above the top edge of posterior
            elif percentile >= 1:
                chis.append(123)  # truth is below the bottom edge of posterior
            else:
                chis.append(scipy.special.erfinv(2 * percentile - 1))
        # error function is defined strangely.  have to multiply by sqrt(2)
        # (it's doesn't have the 2 in the exp, like a normal gaussian)
        chis = np.array(chis) * np.sqrt(2)
        print('new chis', chis)

        ax2.set_xlabel(param + ' vs truth (sigma)', fontsize=14)
        lower = -5
        upper = 5
        if len(errors) > 0:
            plt.hist(
                chis,
                range=(lower, upper),
                bins=int(2 * (upper - lower)),
                cumulative=False,
                density=True,
                histtype='step',
                color='black',
                lw=2,
                zorder=2,
                label='retrieved',
            )
            plt.title(
                'cumulative histogram of ' + str(len(errors)) + ' planets'
            )
            ax2.set_xlim(lower, upper)
        #  add a normal dist, for comparison
        # xdist = np.linspace(-10,10,1000)
        # normaldist = np.exp(xdist**2 / 2)
        # normaldist = scipy.stats.norm.pdf(np.random.randn(100000))
        normaldist = scipy.stats.norm.rvs(size=100000)
        # print('average',np.mean(normaldist))
        # print('    std',np.std(normaldist))
        plt.hist(
            normaldist,
            range=(lower, upper),
            bins=int(2 * (upper - lower)),
            cumulative=False,
            density=True,
            histtype='stepfilled',
            color='lightgrey',
            zorder=1,
            label='ideal',
        )
        plt.legend()
        # ax2.set_ylim(0, 1)
        ax2.set_ylabel('fraction of runs', fontsize=14)

        # Posterior distributions
        ax3 = figure.add_subplot(1, 3, 3)

        ax3.set_xlabel(param, fontsize=14)

        # extend the range just a bit, so that the hist edges don't show
        spread = prior_ranges[param][1] - prior_ranges[param][0]
        lower = prior_ranges[param][0] - spread / 10
        upper = prior_ranges[param][1] + spread / 10

        # plot the posteriors for each noise instantiation
        tracetotal = []
        for traces in all_traces:
            trace = traces[all_keys.index(param)]
            tracetotal.extend(trace)

            plt.hist(
                trace,
                range=(lower, upper),
                bins=120,
                cumulative=False,
                density=True,
                histtype='step',
                color='black',
                lw=0.5,
                zorder=1,
            )
        plt.hist(
            tracetotal,
            range=(lower, upper),
            bins=120,
            cumulative=False,
            density=True,
            histtype='step',
            color='black',
            lw=2,
            zorder=2,
        )
        # vertical line for the truth
        ylim = ax3.get_ylim()
        plt.plot(
            [truth_values[param][0], truth_values[param][0]],
            [0, 10],
            c='k',
            ls=':',
            lw=1,
        )
        ax3.set_ylim(ylim)
        ax3.set_xlim(prior_ranges[param][0], prior_ranges[param][1])

        ax3.set_ylabel('posterior probability', fontsize=14)

        figure.tight_layout()

        if savetodisk:
            plt.savefig(
                saveDir
                + 'fitVStruth_'
                + filt
                + '_'
                + param.replace('/', ':')
                + '.png'
            )
        plot_statevectors.append(save_plot_tosv(figure))
        plt.close(figure)
    return paramlist, plot_statevectors
