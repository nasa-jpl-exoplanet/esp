'''cerberus plotting ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# import excalibur
# from excalibur.ariel.metallicity import massMetalRelation
from excalibur.system.core import ssconstants
from excalibur.util.plotters import save_plot_tosv

import logging


log = logging.getLogger(__name__)

ssc = ssconstants(mks=True)


# --------------------------------------------------------------------
def plot_fits_vs_truths(
    truth_values,
    fit_values,
    fit_errors,
    prior_ranges,
    filt,
    saveDir,
    savetodisk=False,
):
    '''
    Compare the retrieved values against the original inputs
    Also (optionally) show a histogram of the uncertainty values
    '''

    # for ppt, stack the two panels on top of each other
    switch_to_vert_stack = False

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

        if switch_to_vert_stack:
            figure = plt.figure(figsize=(5, 9))
            ax = figure.add_subplot(2, 1, 1)
        else:
            figure = plt.figure(figsize=(11, 5))
            ax = figure.add_subplot(1, 2, 1)

        for truth, fit, error in zip(
            truth_values[param], fit_values[param], fit_errors[param]
        ):
            # check whether there is any real information beyond the prior
            # let's say you have to improve uncertainty by a factor of 2
            # but note that the original 1-sigma uncertainty is ~2/3 of prior range
            # oh wait also note that errorbar is one-sided, so another factor of 2
            # 8/1/24 make the criteria more liberal; it's excluding more than just prior-only guys
            #  let's say 80% of prior range, rather than 20%
            minInfo = 0.8 * 0.68 * 0.5
            newInfo = False
            if param not in prior_ranges:
                if param == '[N/O]':
                    priorRangeDiff = 12
                    if error < minInfo * priorRangeDiff:
                        newInfo = True
                else:
                    log.warning(
                        '--< Cerb.analysis: Parameter missing from prior_range 1: %s >--',
                        param,
                    )
            elif param == 'T':
                priorRangeFactor = (
                    prior_ranges[param][1] / prior_ranges[param][0]
                )
                # prior is normally set to 0.75-1.5 times Teq
                if error < minInfo * (priorRangeFactor - 1) * 0.75 * truth:
                    newInfo = True
            else:
                priorRangeDiff = prior_ranges[param][1] - prior_ranges[param][0]
                if error < minInfo * priorRangeDiff:
                    newInfo = True
            if newInfo:
                clr = 'k'
                lwid = 1
                zord = 4
                ptsiz = 40 / 2
            else:
                clr = 'grey'
                lwid = 0.5
                zord = 2
                ptsiz = 10 / 2

            ax.scatter(
                truth,
                fit,
                facecolor=clr,
                edgecolor=clr,
                s=ptsiz,
                zorder=zord + 1,
            )
            ax.errorbar(
                truth, fit, yerr=error, fmt='.', color=clr, lw=lwid, zorder=zord
            )

            # if param=='T' and truth > 3333:
            #    print('strangely high T in plot',truth)
            # if param=='[X/H]' and truth > 66:
            #    print('strangely high [X/H] in plot',truth)
            # if param=='[C/O]' and truth > 0.5:
            #    print('strangely high [C/O] in plot',truth)

        ax.set_xlabel(param + ' truth', fontsize=14)
        ax.set_ylabel(param + ' fit', fontsize=14)

        xrange = ax.get_xlim()
        # overallmin = min(ax.get_xlim()[0],ax.get_ylim()[0])
        overallmax = max(ax.get_xlim()[1], ax.get_ylim()[1])

        # plot equality as a dashed diagonal line
        ax.plot([-10, 10000], [-10, 10000], 'k--', lw=1, zorder=1)
        if param == 'T':  # show T prior (from 0.75 to 1.5 times Teq)
            ax.plot(
                [-10, 10000], [-10 * 0.75, 10000 * 0.75], 'k:', lw=1, zorder=1
            )
            ax.plot(
                [-10, 10000], [-10 * 1.5, 10000 * 1.5], 'k:', lw=1, zorder=1
            )

        # plot C/O=1 as a dotted vertical line
        if param == '[C/O]':
            solarCO = np.log10(0.55)
            ax.plot([-solarCO, -solarCO], [-100, 100], 'k--', lw=1, zorder=3)
            plt.text(
                -solarCO + 0.03,
                -5,
                'C/O=1',
                c='black',
                rotation='vertical',
                va='center',
                fontsize=12,
            )

        # ax.set_xlim(overallmin,overallmax)
        # ax.set_ylim(overallmin,overallmax)
        if param == 'T':  # the prior for T varies between targets
            ax.set_xlim(0, overallmax)
            ax.set_ylim(0, overallmax)
        elif param not in prior_ranges:
            ax.set_xlim(xrange)
            if param == '[N/O]':
                ax.set_ylim(-6, 6)
            else:
                log.warning(
                    '--< Cerb.analysis: Parameter missing from prior_range 2: %s >--',
                    param,
                )
        else:
            # actually, don't use prior range for X/H and X/O on x-axis
            # ax.set_xlim(prior_ranges[param][0],prior_ranges[param][1])
            ax.set_xlim(xrange)
            ax.set_ylim(prior_ranges[param][0], prior_ranges[param][1])

        if param == 'T':
            plt.title(
                'temperature retrieval for '
                + str(len(fit_errors[param]))
                + ' planets'
            )
        elif param == '[X/H]':
            plt.title(
                'metallicity retrieval for '
                + str(len(fit_errors[param]))
                + ' planets'
            )
        elif param == '[C/O]':
            plt.title(
                'C/O retrieval for ' + str(len(fit_errors[param])) + ' planets'
            )
        else:
            plt.title(
                param
                + ' retrieval for '
                + str(len(fit_errors[param]))
                + ' planets'
            )

        # UNCERTAINTY HISTOGRAMS IN SECOND PANEL
        if switch_to_vert_stack:
            ax2 = figure.add_subplot(2, 1, 2)
        else:
            ax2 = figure.add_subplot(1, 2, 2)
        if param == 'T':
            errors = np.array(fit_errors[param]) / np.array(fit_values[param])
            ax2.set_xlabel(param + ' fractional uncertainty', fontsize=14)
        else:
            errors = np.array(fit_errors[param])
            ax2.set_xlabel(param + ' uncertainty', fontsize=14)
        if len(errors) > 0:
            # the histogram range has to go past the data range or you get a vertical line on the right
            lower = errors.min() / 1.5
            upper = errors.max() * 1.5
            # print('uncertainty range (logged)',param,lower,upper)
            plt.hist(
                errors,
                range=(lower, upper),
                bins=1000,
                cumulative=True,
                density=True,
                histtype='step',
                color='black',
                zorder=1,
                label='',
            )
            plt.title(
                'cumulative histogram of ' + str(len(errors)) + ' planets'
            )
            ax2.semilogx()
            ax2.set_xlim(lower, upper)
        ax2.set_ylim(0, 1)
        ax2.set_ylabel('fraction of planets', fontsize=14)

        figure.tight_layout()

        # ('display' doesn't work for pdf files)
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
    return plot_statevectors


# --------------------------------------------------------------------
def plot_fit_uncertainties(
    fit_values, fit_errors, prior_ranges, filt, saveDir, savetodisk=False
):
    '''
    Plot uncertainty as a function of the fit value
    And show a histogram of the uncertainty values
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
    plot_statevectors = []
    for param in paramlist:

        figure = plt.figure(figsize=(11, 5))
        ax = figure.add_subplot(1, 2, 1)

        for fitvalue, error in zip(fit_values[param], fit_errors[param]):

            # check whether there is any real information beyond the prior
            # let's say you have to improve uncertainty by a factor of 2
            # but note that the original 1-sigma uncertainty is ~2/3 of prior range
            # oh wait also note that errorbar is one-sided, so another factor of 2
            # 8/1/24 make the criteria more liberal; it's excluding more than just prior-only guys
            #  let's say 80% of prior range, rather than 20%
            minInfo = 0.8 * 0.68 * 0.5
            newInfo = False
            if param not in prior_ranges:
                if param == '[N/O]':
                    priorRangeDiff = 12
                    if error < minInfo * priorRangeDiff:
                        newInfo = True
                else:
                    log.warning(
                        '--< Cerb.analysis: Parameter missing from prior_range 3: %s >--',
                        param,
                    )
            elif param == 'T':
                priorRangeFactor = (
                    prior_ranges[param][1] / prior_ranges[param][0]
                )
                # prior is normally set to 0.75-1.5 times Teq
                # if error < minInfo * (priorRangeFactor-1) * 0.75*truth:
                # asdf: this needs work maybe.  should pass in Teq as truth?
                if error < minInfo * (priorRangeFactor - 1) * 0.75 * fitvalue:
                    newInfo = True
            else:
                priorRangeDiff = prior_ranges[param][1] - prior_ranges[param][0]
                if error < minInfo * priorRangeDiff:
                    newInfo = True
            if newInfo:
                clr = 'k'
                # lwid = 1
                zord = 4
                ptsiz = 40
            else:
                clr = 'grey'
                # lwid = 0.5
                zord = 2
                ptsiz = 10
            ax.scatter(
                fitvalue,
                error,
                facecolor=clr,
                edgecolor=clr,
                s=ptsiz,
                zorder=zord + 1,
            )

        ax.set_xlabel(param + ' fit value', fontsize=14)
        ax.set_ylabel(param + ' fit uncertainty', fontsize=14)

        if param == 'T':
            plt.title(
                'temperature retrieval for '
                + str(len(fit_errors[param]))
                + ' planets'
            )
        elif param == '[X/H]':
            plt.title(
                'metallicity retrieval for '
                + str(len(fit_errors[param]))
                + ' planets'
            )
        elif param == '[C/O]':
            plt.title(
                'C/O retrieval for ' + str(len(fit_errors[param])) + ' planets'
            )
        else:
            plt.title(
                param
                + ' retrieval for '
                + str(len(fit_errors[param]))
                + ' planets'
            )

        # plot C/O=1 as a dotted vertical line
        if param == '[C/O]':
            yrange = ax.get_ylim()
            solarCO = np.log10(0.55)
            ax.plot([-solarCO, -solarCO], [-100, 100], 'k--', lw=1, zorder=3)
            ax.set_ylim(yrange)
            plt.text(
                -solarCO + 0.03,
                ax.get_ylim()[1] - 0.4,
                'C/O=1',
                c='black',
                rotation='vertical',
                va='center',
                fontsize=12,
            )

        # UNCERTAINTY HISTOGRAMS IN SECOND PANEL
        ax2 = figure.add_subplot(1, 2, 2)
        if param == 'T':
            errors = np.array(fit_errors[param]) / np.array(fit_values[param])
            ax2.set_xlabel(param + ' fractional uncertainty', fontsize=14)
        else:
            errors = np.array(fit_errors[param])
            ax2.set_xlabel(param + ' uncertainty', fontsize=14)
        # the histogram range has to go past the data range or you get a vertical line on the right
        lower = errors.min() / 2.0
        upper = errors.max() * 2.0
        # print('uncertainty range (logged)',param,lower,upper)
        plt.hist(
            errors,
            range=(lower, upper),
            bins=1000,
            cumulative=True,
            density=True,
            histtype='step',
            color='black',
            zorder=1,
            label='',
        )
        plt.title('cumulative histogram of ' + str(len(errors)) + ' planets')
        ax2.semilogx()
        ax2.set_xlim(lower, upper)
        ax2.set_ylim(0, 1)
        ax2.set_ylabel('fraction of planets', fontsize=14)

        figure.tight_layout()

        # ('display' doesn't work for pdf files)
        if savetodisk:
            plt.savefig(
                saveDir
                + 'fitUncertainties_'
                + filt
                + '_'
                + param.replace('/', ':')
                + '.png'
            )
        plot_statevectors.append(save_plot_tosv(figure))
        plt.close(figure)
    return plot_statevectors


# --------------------------------------------------------------------
def plot_massFits(
    masses,
    stellarFEHs,
    truth_values,
    fit_values,
    fit_errors,
    prior_ranges,
    filt,
    saveDir,
    plot_truths=False,  # for Ariel-sims, include truth as open circles?
    savetodisk=False,
):
    '''how well do we retrieve the input mass-metallicity relation?'''

    MEarth = 5.972e27 / 1.898e30

    figure = plt.figure(figsize=(11, 5))
    # with text labels hanging off the right side, need to stretch the figure size
    # figure = plt.figure(figsize=(17,5))
    ax = figure.add_subplot(1, 2, 1)

    # Note: masses_true is not actually used, just masses. (they should be the same thing)
    masses_true = truth_values['Mp']
    metals_true = truth_values['[X/H]']
    metals_fit = fit_values['[X/H]']
    metals_fiterr = fit_errors['[X/H]']
    # switch to 2-sided (asymmetric) error bars
    metals_fiterr = [(f[1] + f[0]) / 2 for f in fit_errors['[X/H]']]
    metals_fiterrhi = [f[1] for f in fit_errors['[X/H]']]
    metals_fiterrlo = [f[0] for f in fit_errors['[X/H]']]

    ilab1 = False
    ilab2 = False
    for (
        mass,
        masstrue,
        metaltrue,
        metalfit,
        metalerror,
        metalerrorlo,
        metalerrorhi,
    ) in zip(
        masses,
        masses_true,
        metals_true,
        metals_fit,
        metals_fiterr,
        metals_fiterrlo,
        metals_fiterrhi,
    ):
        # check whether there is any real information beyond the prior
        # let's say you have to improve uncertainty by a factor of 2
        # but note that the original 1-sigma uncertainty is ~2/3 of prior range
        # oh wait also note that errorbar is one-sided, so another factor of 2
        # 8/1/24 make the criteria more liberal; it's excluding more than just prior-only guys
        #  let's say 80% of prior range, rather than 20%
        minInfo = 0.8 * 0.68 * 0.5
        priorRangeDiff = prior_ranges['[X/H]'][1] - prior_ranges['[X/H]'][0]
        if metalerror < minInfo * priorRangeDiff:
            clr = 'k'
            lwid = 1
            ptsiz = 40
            zord = 5
        else:
            clr = 'grey'
            lwid = 0.5
            ptsiz = 10
            zord = 2
        if 'sim' in filt:
            ptsiz /= 2  # make smaller points for the large Ariel sample

        if 1:
            if plot_truths and metaltrue not in (666, 666666):
                if not ilab1 and clr == 'k':
                    ilab1 = True
                    ax.scatter(
                        masstrue,
                        metaltrue,
                        label='true value',
                        facecolor='w',
                        edgecolor=clr,
                        s=ptsiz,
                        zorder=zord + 1,
                    )
                else:
                    ax.scatter(
                        masstrue,
                        metaltrue,
                        facecolor='w',
                        edgecolor=clr,
                        s=ptsiz,
                        zorder=zord + 1,
                    )
            if not ilab2 and clr == 'k':
                ilab2 = True
                ax.scatter(
                    mass,
                    metalfit,
                    label='retrieved metallicity',
                    facecolor=clr,
                    edgecolor=clr,
                    s=ptsiz,
                    zorder=zord + 2,
                )
            else:
                ax.scatter(
                    mass,
                    metalfit,
                    facecolor=clr,
                    edgecolor=clr,
                    s=ptsiz,
                    zorder=zord + 2,
                )
            # allow for asymmetric error bars!!!
            # ax.errorbar(mass, metalfit, yerr=metalerror,
            #             fmt='.', color=clr, lw=lwid, zorder=zord)
            # if clr=='k':
            #    print('x y',mass,metalfit)
            #    print('  old:',metalerror)
            #    print('  new:',metalerrorlo,metalerrorhi)
            ax.errorbar(
                mass,
                metalfit,
                yerr=np.array([[metalerrorlo], [metalerrorhi]]),
                fmt='.',
                color=clr,
                lw=lwid,
                zorder=zord,
            )
    ax.semilogx()
    ax.set_xlabel('$M_p \\, (M_{\\rm Jup})$', fontsize=14)
    ax.set_ylabel('[X/H]$_p$', fontsize=14)
    xrange = ax.get_xlim()
    yrange = ax.get_ylim()
    # xrange = (0.01, 10.0)

    # plot the underlying distribution (only if this is a simulation)
    # actually, also plot Thorngren relationship for real HST data
    massesThorngren = np.logspace(-5, 3, 100)
    metalsThorngren = massMetalRelation(0, massesThorngren, thorngren=True)
    if 'sim' in filt:
        # ax.plot(massesThorngren,metalsThorngren, 'k:', lw=1, zorder=1, label='true relationship')
        ax.plot(
            massesThorngren,
            metalsThorngren,
            'k:',
            lw=1,
            zorder=1,
            label='Thorngren+ 2016',
        )
    else:
        ax.plot(
            massesThorngren,
            metalsThorngren,
            'k:',
            lw=1,
            zorder=1,
            label='Thorngren+ 2016',
        )
        
    plt.legend()

    # use the prior range for the y-axis
    yrange = (prior_ranges['[X/H]'][0], prior_ranges['[X/H]'][1])

    # display the fit parameters (and Thorgran too)
    fit_mass_exp = polynomialCoeffs[0]
    fit_mass_exp_err = np.sqrt(covariance[0, 0])
    fit_metal_dex = polynomialCoeffs[1]
    fit_metal_dex_err = np.sqrt(covariance[1, 1])
    fit_metal_linear = 10.0**fit_metal_dex
    fit_metal_linear_err = fit_metal_linear * (10.0**fit_metal_dex_err - 1)
    # print('fit result  mass-exp',fit_mass_exp,fit_mass_exp_err)
    # print('fit result  metal(dex)',fit_metal_dex,fit_metal_dex_err)
    # print('fit result  metal(lin)',fit_metal_linear,fit_metal_linear_err)
    # resultsstring = '[X/H]$_p$ = (' + f'{fit_metal_dex:3.2f}' + \
    #    '$\\pm$' + f'{fit_metal_dex_err]):3.2f}' + \
    resultsstring = (
        'Z$_p$ = ('
        + f'{fit_metal_linear:3.2f}'
        + '$\\pm$'
        + f'{fit_metal_linear_err:3.2f}'
        + ') $M_p^{'
        + f'{fit_mass_exp:3.2f}'
        + '\\pm'
        + f'{fit_mass_exp_err:3.2f}'
        + '}$ (fit)'
    )
    # print('resultsstring',resultsstring)
    plt.text(
        xrange[0] * 1.2,
        yrange[0] + 0.8,
        resultsstring,
        c='black',
        ha='left',
        fontsize=10,
    )
    plt.text(
        xrange[0] * 1.2,
        yrange[0] + 0.2,
        'Z$_p$ = (9.7$\\pm$1.3) $M_p^{-0.45\\pm0.09}$ (Thorngren)',
        # '[X/H]$_p$ = (0.97$\\pm$0.05) $M_p^{-0.45\\pm0.09}$ (Thorngren)',
        # '[X/H]$_p$ = (9.7$\\pm$1.3) $M_p^{-0.45\\pm0.09}$ (Thorngren)',
        c='black',
        ha='left',
        fontsize=10,
    )

    ax.set_xlim(xrange)
    ax.set_ylim(yrange)

    # SECOND PANEL - same thing but subtract off the stellar metallicity
    ax2 = figure.add_subplot(1, 2, 2)

    ilab1 = False
    ilab2 = False
    for (
        stellarFEH,
        mass,
        masstrue,
        metaltrue,
        metalfit,
        metalerror,
        metalerrorlo,
        metalerrorhi,
    ) in zip(
        stellarFEHs,
        masses,
        masses_true,
        metals_true,
        metals_fit,
        metals_fiterr,
        metals_fiterrlo,
        metals_fiterrhi,
    ):
        if metalerror < minInfo * priorRangeDiff:
            clr = 'k'
            lwid = 1
            ptsiz = 40
            zord = 5
        else:
            clr = 'grey'
            lwid = 0.5
            ptsiz = 10
            zord = 2
        if 'sim' in filt:
            ptsiz /= 2  # make smaller points for the large Ariel sample
        if 1:
            if plot_truths and metaltrue not in (666, 666666):
                if not ilab1:
                    ilab1 = True
                    ax2.scatter(
                        masstrue,
                        metaltrue,
                        label='true value',
                        facecolor='w',
                        edgecolor=clr,
                        s=ptsiz,
                        zorder=zord + 1,
                    )
                else:
                    ax2.scatter(
                        masstrue,
                        metaltrue,
                        facecolor='w',
                        edgecolor=clr,
                        s=ptsiz,
                        zorder=zord + 1,
                    )
            if not ilab2 and clr == 'k':
                ilab2 = True
                ax2.scatter(
                    mass,
                    metalfit - stellarFEH,
                    label='retrieved metallicity',
                    facecolor=clr,
                    edgecolor=clr,
                    s=ptsiz,
                    zorder=zord + 2,
                )
            else:
                ax2.scatter(
                    mass,
                    metalfit - stellarFEH,
                    facecolor=clr,
                    edgecolor=clr,
                    s=ptsiz,
                    zorder=zord + 2,
                )
            # ax2.errorbar(mass, metalfit - stellarFEH, yerr=metalerror,
            #             fmt='.', color=clr, lw=lwid, zorder=zord)
            # allow for asymmetric error bars!!!
            ax2.errorbar(
                mass,
                metalfit - stellarFEH,
                yerr=np.array([[metalerrorlo], [metalerrorhi]]),
                fmt='.',
                color=clr,
                lw=lwid,
                zorder=zord,
            )
    ax2.semilogx()
    ax2.set_xlabel('$M_p (M_{\\rm Jup})$', fontsize=14)
    ax2.set_ylabel('[X/H]$_p$ - [X/H]$_\\star$', fontsize=14)
    xrange = ax2.get_xlim()
    yrange = ax2.get_ylim()
    # xrange = (0.01, 10.0)

    # plot the underlying distribution (only if this is a simulation)
    # actually, also plot Thorngren relationship for real HST data
    if 'sim' in filt:
        # ax2.plot(massesThorngren,metalsThorngren, 'k:', lw=1, zorder=1, label='true relationship')
        #  careful here. 'true relation' might imply that all planets fall on that dotted line
        ax2.plot(
            massesThorngren,
            metalsThorngren,
            'k:',
            lw=1,
            zorder=1,
            label='Thorngren+ 2016',
        )
    else:
        ax2.plot(
            massesThorngren,
            metalsThorngren,
            'k:',
            lw=1,
            zorder=1,
            label='Thorngren+ 2016',
        )

    plt.legend()

    ax2.set_xlim(xrange)
    ax2.set_ylim(yrange)
    figure.tight_layout()  # some trouble with this, with bigger figure maybe?

    # ('display' doesn't work for pdf files)
    if savetodisk:
        plt.savefig(saveDir + 'massVSmetals_' + filt + '.png')
    return save_plot_tosv(figure), figure
