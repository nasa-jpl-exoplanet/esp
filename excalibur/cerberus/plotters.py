'''cerberus plotting ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import logging
import corner
import numpy as np
import matplotlib.pyplot as plt

# import excalibur
from excalibur.ariel.metallicity import massMetalRelation
from excalibur.system.core import ssconstants
from excalibur.util.plotters import save_plot_tosv

log = logging.getLogger(__name__)

ssc = ssconstants(mks=True)


# --------------------------------------------------------------------
def rebin_data(transitdata, binsize=4):
    '''
    rebin a spectrum (blue points calculated from grey ones)
    set uncertainty to median/sqrt(N); ignore correlated/systematic errors
    this is ok, since it's just for visualization; not used in the analysis
    '''

    nspec = int(transitdata['wavelength'].size / binsize)
    minspec = np.nanmin(transitdata['wavelength'])
    maxspec = np.nanmax(transitdata['wavelength'])
    # add a small bit (1.e-10) to max edge, such that "< wavehi" includes last point
    wavebinedges = np.linspace(minspec, maxspec + 1.0e-10, nspec + 1)

    specbin = []
    errbin = []
    wavebin = []
    for wavelo, wavehi in zip(wavebinedges[:-1], wavebinedges[1:]):
        select = transitdata['wavelength'] < wavehi
        select = select & (transitdata['wavelength'] >= wavelo)
        select = select & np.isfinite(transitdata['depth'])

        if np.sum(np.isfinite(transitdata['depth'][select])) > 0:
            specbin.append(
                np.nansum(
                    transitdata['depth'][select]
                    / (transitdata['error'][select] ** 2)
                )
                / np.nansum(1.0 / (transitdata['error'][select] ** 2))
            )
            errbin.append(
                np.nanmedian(transitdata['error'][select])
                / np.sqrt(np.sum(select))
            )
            wavebin.append(
                np.nansum(
                    transitdata['wavelength'][select]
                    / transitdata['error'][select] ** 2
                )
                / np.nansum(1.0 / transitdata['error'][select] ** 2)
            )
    transitdata['binned_wavelength'] = np.array(wavebin)
    transitdata['binned_depth'] = np.array(specbin)
    transitdata['binned_error'] = np.array(errbin)

    return transitdata


# --------------------------------------------------------------------
# --------------------------------------------------------------------
def plot_spectrumfit(
    transitdata,
    patmos_model,
    patmos_modelProfiled,
    patmos_bestfit,
    fmcarray,
    truth_spectrum,
    system_data,
    ancillary_data,
    atmos_data,
    filt,
    modelName,
    trgt,
    p,
    saveDir='./',
    savetodisk=False,
):
    '''plot the best fit to the data'''

    include_fit_range_as_grey_lines = False

    # figure, ax = plt.subplots(figsize=(8,4))
    figgy = plt.figure(figsize=(20, 4))
    # figgy = plt.figure(figsize=(8,4))
    figgy.subplots_adjust(
        left=0.05, right=0.7, bottom=0.15, top=0.93, wspace=0.8
    )
    ax = plt.subplot(1, 2, 1)

    okPart = np.where(np.isfinite(transitdata['depth']))

    # 1) plot the data
    plt.errorbar(
        transitdata['wavelength'],
        transitdata['depth'] * 100,
        yerr=transitdata['error'] * 100,
        fmt='.',
        color='lightgray',
        zorder=1,
        label='raw data',
    )
    # 2) also plot the rebinned data points
    plt.errorbar(
        transitdata['binned_wavelength'],
        transitdata['binned_depth'] * 100,
        yerr=transitdata['binned_error'] * 100,
        # fmt='o', color='k', markeredgecolor='k', markerfacecolor='w', zorder=5,
        # fmt='^', color='blue', zorder=5,
        # fmt='o', color='royalblue', zorder=5,
        fmt='o',
        markeredgecolor='k',
        color='None',
        ecolor='k',
        zorder=5,
        label='rebinned data',
    )
    # 3a) plot the best-fit model - old posterior-median method
    # plt.plot(transitdata['wavelength'],
    #         patmos_modelProfiled * 100,
    #         # c='k', lw=2, zorder=4,
    #         c='orange', lw=2, ls='--', zorder=4,
    #         label='best fit (old method)')
    # 3b) plot the best-fit model - new random selection parameter-set checking
    plt.plot(
        transitdata['wavelength'][okPart],
        patmos_bestfit * 100,
        # c='k', lw=2, zorder=4,
        c='orange',
        lw=2,
        zorder=4,
        label='best fit',
    )
    # 4) plot a selection of walkers, to see spread
    xlims = plt.xlim()
    ylims = plt.ylim()
    # print('median pmodel',np.nanmedian(patmos_model))
    if not include_fit_range_as_grey_lines:
        fmcarray = []
    for patmos_modeli in fmcarray:
        plt.plot(
            transitdata['wavelength'][okPart],
            patmos_modeli * 100,
            c='grey',
            lw=0.2,
            zorder=2,
        )
    plt.ylim(ylims)  # revert to the original y-bounds, in case messed up
    # 5) plot the true spectrum, if it is a simulation
    if truth_spectrum is not None:
        plt.plot(
            truth_spectrum['wavelength'],
            truth_spectrum['depth'] * 100,
            c='k',
            lw=1.5,
            zorder=3,
            # c='orange', lw=2, zorder=3,
            label='truth',
        )

    offsets_model = (patmos_model - transitdata['depth'][okPart]) / transitdata[
        'error'
    ][okPart]
    offsets_modelProfiled = (
        patmos_modelProfiled - transitdata['depth'][okPart]
    ) / transitdata['error'][okPart]

    # the 'average' function (which allows for weights) doesn't have a NaN version,
    #  so mask out any NaN regions by hand
    flatlineFit = np.average(
        transitdata['depth'][okPart],
        weights=1 / transitdata['error'][okPart] ** 2,
    )
    # print('flatline average',flatlineFit)
    offsets_flat = (flatlineFit - transitdata['depth']) / transitdata['error']
    # print('median chi2, flat',np.nanmedian(offsets_flat**2))

    chi2model = np.nansum(offsets_model**2)
    chi2modelProfiled = np.nansum(offsets_modelProfiled**2)
    chi2flat = np.nansum(offsets_flat**2)

    numParam_model = 8
    numParam_truth = 0
    numParam_flat = 1

    numPoints = len(patmos_model)
    # print('numpoints',numPoints)
    chi2model_red = chi2model / (numPoints - numParam_model)
    chi2modelProfiled_red = chi2modelProfiled / (numPoints - numParam_model)
    chi2flat_red = chi2flat / (numPoints - numParam_flat)

    # add some labels off to the right side
    xoffset = 1.2
    if truth_spectrum is not None:
        offsets_truth = (
            truth_spectrum['depth'] - transitdata['depth']
        ) / transitdata['error']
        chi2truth = np.nansum(offsets_truth**2)
        chi2truth_red = chi2truth / (numPoints - numParam_truth)
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.6,
            '$\\chi^2$-truth=' + f"{chi2truth:5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.36,
            '$\\chi^2_{red}$-truth=' + f"{chi2truth_red:5.2f}",
            fontsize=12,
        )
    if 'tier' in atmos_data:
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 1.00,
            'Tier=' + str(atmos_data['tier']),
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.93,
            '# of Visits=' + str(atmos_data['visits']),
            fontsize=12,
        )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.8,
        'ZFOM='
        + f"{ancillary_data['ZFOM']:4.1f}"
        + '-'
        + f"{ancillary_data['ZFOM_max']:4.1f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.73,
        'TSM=' + f"{ancillary_data['TSM']:5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.53,
        '$\\chi^2$-model=' + f"{chi2modelProfiled:5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.46,
        '$\\chi^2$-flat=' + f"{chi2flat:5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.29,
        '$\\chi^2_{red}$-model=' + f"{chi2modelProfiled_red:5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.22,
        '$\\chi^2_{red}$-flat=' + f"{chi2flat_red:5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.12,
        '$\\chi^2_{red}$(model/flat)='
        + f"{(chi2modelProfiled_red / chi2flat_red):5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.05,
        '$\\Delta\\chi^2$(flat-model)='
        + f"{(chi2flat - chi2modelProfiled):5.1f}",
        fontsize=12,
    )

    if filt == 'Ariel-sim':
        plt.xlim(0, 8)
    plt.title(trgt + ' ' + p, fontsize=16)
    plt.xlabel(str('Wavelength [$\\mu m$]'), fontsize=14)
    plt.ylabel(str('$(R_p/R_*)^2$ [%]'), fontsize=14)
    plt.legend()
    # add the scale-height comparison back in (on the righthand y-axis)
    if 'H_max' in ancillary_data.keys():
        axtwin = ax.twinx()
        axtwin.set_ylabel('$\\Delta$ [H$_{\\rm s}$]')
        axmin, axmax = ax.get_ylim()
        rp0rs = np.sqrt(np.nanmedian(transitdata['depth']))
        # awkward. H is in km but R* is converted to meters
        # Hs0rs = (ancillary_data['H'] *1.e3) / (system_data['R*'] * ssc['Rsun'])
        Hs0rs = (ancillary_data['H_max'] * 1.0e3) / (
            system_data['R*'] * ssc['Rsun']
        )
        # print('rp0rs',rp0rs)
        # print('hsors',Hs0rs)
        # print('mmw,mmwmin', ancillary_data['mmw'] ,ancillary_data['mmw_min'])
        # print('H,Hmax', ancillary_data['H'] ,ancillary_data['H_max'])
        # print('H R*', ancillary_data['H'] , system_data['R*'])
        if axmin >= 0:
            axtwin.set_ylim(
                (np.sqrt(1e-2 * axmin) - rp0rs) / Hs0rs,
                (np.sqrt(1e-2 * axmax) - rp0rs) / Hs0rs,
            )
        else:
            axtwin.set_ylim(
                (-np.sqrt(-1e-2 * axmin) - rp0rs) / Hs0rs,
                (np.sqrt(1e-2 * axmax) - rp0rs) / Hs0rs,
            )

    #  Plot the pre-profiled result off to the right, if there's been profiling
    if np.any(patmos_model != patmos_modelProfiled):
        # ax2 = plt.subplot(1,2,2)
        plt.subplot(1, 2, 2)

        # 1) plot the data
        plt.errorbar(
            transitdata['wavelength'],
            transitdata['depth'] * 100,
            yerr=transitdata['error'] * 100,
            fmt='.',
            color='lightgray',
            zorder=1,
            label='raw data',
        )
        # 2) also plot the rebinned data points
        plt.errorbar(
            transitdata['binned_wavelength'],
            transitdata['binned_depth'] * 100,
            yerr=transitdata['binned_error'] * 100,
            fmt='o',
            markeredgecolor='k',
            color='None',
            ecolor='k',
            zorder=5,
            label='rebinned data',
        )
        # 3a) plot the best-fit model - old posterior-median method
        # plt.plot(transitdata['wavelength'],
        #         patmos_model * 100,
        #         c='orange', lw=2, ls='--', zorder=4,
        #         label='best fit (old method)')
        # 3b) plot the best-fit model - new random selection parameter-set checking
        plt.plot(
            transitdata['wavelength'][okPart],
            patmos_bestfit * 100,
            c='orange',
            lw=2,
            zorder=4,
            label='best fit',
        )
        xlims = plt.xlim()
        ylims = plt.ylim()
        xoffset = (
            0.03  # only a small offset is needed if Hscale twin axis not used
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.53,
            '$\\chi^2$-model=' + f"{chi2model:5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.46,
            '$\\chi^2$-flat=' + f"{chi2flat:5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.29,
            '$\\chi^2_{red}$-model=' + f"{chi2model_red:5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.22,
            '$\\chi^2_{red}$-flat=' + f"{chi2flat_red:5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.12,
            '$\\chi^2_{red}$(model/flat)='
            + f"{(chi2model_red / chi2flat_red):5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.05,
            '$\\Delta\\chi^2$(flat-model)=' + f"{(chi2flat - chi2model):5.1f}",
            fontsize=12,
        )
        if filt == 'Ariel-sim':
            plt.xlim(0, 8)
        plt.title(trgt + ' ' + p + ' (no profiling)', fontsize=16)
        plt.xlabel(str('Wavelength [$\\mu m$]'), fontsize=14)
        plt.ylabel(str('$(R_p/R_*)^2$ [%]'), fontsize=14)
        plt.legend()

    # figgy.tight_layout()
    # plt.show()
    if savetodisk:
        # pdf is so much better, but xv gives error (stick with png for debugging)
        plt.savefig(
            saveDir
            + 'bestFit_'
            + filt
            + '_'
            + modelName
            + '_'
            + trgt
            + ' '
            + p
            + '.png'
        )

    return save_plot_tosv(figgy), figgy


# --------------------------------------------------------------------


def plot_corner(
    allkeys,
    alltraces,
    profiletraces,
    modelParams_bestFit,
    truth_params,
    prior_ranges,
    filt,
    modelName,
    trgt,
    p,
    bins=10,
    saveDir='./',
    savetodisk=False,
    verbose=False,
):
    '''corner plot showing posterior distributions'''

    truthcolor = 'darkgreen'
    fitcolor = 'firebrick'

    paramValues_bestFit = []
    if modelParams_bestFit:
        tpr, ctp, hazescale, hazeloc, hazethick, tceqdict, mixratio = (
            modelParams_bestFit
        )
        # print('model param in corner plot',modelParams_bestFit)

        for param in allkeys:
            if param == 'T':
                paramValues_bestFit.append(tpr)
            elif param == 'CTP':
                paramValues_bestFit.append(ctp)
            elif param == 'HScale':
                paramValues_bestFit.append(hazescale)
            elif param == 'HLoc':
                paramValues_bestFit.append(hazeloc)
            elif param == 'HThick':
                paramValues_bestFit.append(hazethick)
            elif param == '[X/H]':
                paramValues_bestFit.append(tceqdict['XtoH'])
            elif param == '[C/O]':
                paramValues_bestFit.append(tceqdict['CtoO'])
            elif param == '[N/O]':
                paramValues_bestFit.append(tceqdict['NtoO'])
            elif mixratio and param in mixratio:
                paramValues_bestFit.append(mixratio[param])
            else:
                log.warning('--< ERROR: param not in list: %s >--', param)
        # print('best fit values in corner plot',paramValues_bestFit)

    mcmcMedian = np.nanmedian(np.array(profiletraces), axis=1)
    # print(' params inside of corner plotting',allkeys)
    # print('medians inside of corner plotting',mcmcMedian)
    # print('bestfit inside of corner plotting',paramValues_bestFit)
    lo = np.nanpercentile(np.array(alltraces), 16, axis=1)
    hi = np.nanpercentile(np.array(alltraces), 84, axis=1)
    # span = hi - lo
    # Careful!  These are not actually the prior ranges; they're the range of walker values
    priorlo = np.nanmin(np.array(alltraces), axis=1)
    priorhi = np.nanmax(np.array(alltraces), axis=1)
    # OK fixed now. prior ranges are saved as output from atmos and then passed in here
    for ikey, key in enumerate(allkeys):
        # print('param:',key)
        # print(' old prior range:',key,priorlo[ikey],priorhi[ikey])

        # special case for older HST state vectors (probably not needed anymore)
        if 'Hscale' in prior_ranges:
            prior_ranges['HScale'] = prior_ranges['Hscale']

        if key in prior_ranges.keys():
            priorlo[ikey], priorhi[ikey] = prior_ranges[key]
        # else:
        #    print('TROUBLE: param not found',prior_ranges.keys())
        # print(' new prior range:',key,priorlo[ikey],priorhi[ikey])
    # priorspan = priorhi - priorlo
    # priormid = (priorhi + priorlo) / 2.

    # put a line showing the equilibrium temperature
    # eqtemp = orbp['T*']*np.sqrt(orbp['R*']*ssc['Rsun/AU']/(2.*orbp[p]['sma']))
    # print(lo[4], mcmcMedian[4], hi[4], 'Teq', eqtemp)

    lodist = np.array(mcmcMedian) - np.array(lo)
    hidist = np.array(hi) - np.array(mcmcMedian)
    lorange = np.array(mcmcMedian) - 3 * lodist
    hirange = np.array(mcmcMedian) + 3 * hidist
    trange = [tuple([x, y]) for x, y in zip(lorange, hirange)]
    # previous lines did 3-sigma range.  better to just use the prior bounds as bounds
    trange = [tuple([x, y]) for x, y in zip(priorlo, priorhi)]
    # print('trange',trange)
    # print('lodist',lodist)
    # print('lorange',lorange)
    truths = None
    if truth_params is not None:
        truths = []
        for thiskey in allkeys:
            if thiskey == 'T':
                truths.append(truth_params['Teq'])
            elif thiskey == '[X/H]':
                # truths.append(np.log10(truth_params['metallicity']))
                truths.append(truth_params['metallicity'])
            elif thiskey == '[C/O]':
                # truths.append(np.log10(truth_params['C/O'] / 0.54951))
                truths.append(truth_params['C/O'])
            elif thiskey == '[N/O]':
                truths.append(0)
            elif thiskey in truth_params.keys():
                truths.append(truth_params[thiskey])
            else:
                truths.append(666666)
                log.warning(
                    '--< PROBLEM: no truth value for this key: %s >--', thiskey
                )
        # print('truths in corner plot',truths)
    if len(allkeys) <= 1:
        # print('corner plot needs at least 2 parameters; doing histogram instead')

        figure = plt.figure(figsize=(5, 4))
        ax = figure.add_subplot(1, 1, 1)

        plt.hist(
            np.vstack(np.array(profiletraces)).T,
            range=trange[0],
            bins=bins,
            cumulative=False,
            density=True,
            histtype='step',
            color='black',
            lw=2,
            zorder=2,
            label='retrieved',
        )
        ylim = ax.get_ylim()
        plt.plot(
            [truths[0], truths[0]],
            [0, 100],
            c='k',
            ls='--',
            lw=1,
            label='truth',
        )
        ax.set_ylim(ylim)

        ax.set_xlim(trange[0])
        ax.set_xlabel(allkeys[0])
        ax.set_ylabel('posterior probability')

    else:
        figure = corner.corner(
            np.vstack(np.array(profiletraces)).T,
            bins=bins,
            labels=allkeys,
            range=trange,
            truths=truths,
            truth_color=truthcolor,
            show_titles=True,
            quantiles=[0.16, 0.50, 0.84],
        )
        # smaller size for corner plot might fit better, but this creates a bit of checkerboarding
        # figure.set_size_inches(16,16)  # this actually makes it smaller

        ndim = len(alltraces)
        axes = np.array(figure.axes).reshape((ndim, ndim))
        # use larger font size for the axis labels
        for i in range(ndim):
            ax = axes[ndim - 1, i]
            ax.set_xlabel(allkeys[i], fontsize=14)
        for i in range(ndim - 1):
            # skipping the first one on the y side (it's a histo, not a 2-D plot)
            ax = axes[i + 1, 0]
            ax.set_ylabel(allkeys[i + 1], fontsize=14)
        # draw a point and crosshair for the best fit in each subpanel
        # the best fit is not defined for the debug call inside atmos, so skip
        if paramValues_bestFit:
            print('troot', truth_params)
            print('allkeys', allkeys)
            print('ndim', ndim)
            print('paramvals', paramValues_bestFit)
            for yi in range(ndim):
                for xi in range(yi):
                    ax = axes[yi, xi]
                    # ax.axvline(mcmcMedian[xi], color=fitcolor)
                    # ax.axhline(mcmcMedian[yi], color=fitcolor)
                    # ax.plot(mcmcMedian[xi], mcmcMedian[yi], marker='s', c=fitcolor)
                    ax.axvline(paramValues_bestFit[xi], color=fitcolor)
                    ax.axhline(paramValues_bestFit[yi], color=fitcolor)
                    ax.plot(
                        paramValues_bestFit[xi],
                        paramValues_bestFit[yi],
                        marker='s',
                        c=fitcolor,
                    )
            for i in range(ndim):
                ax = axes[i, i]
                # draw light-colored vertical lines in each hisogram for the prior
                #  drop this. adds clutter and is redundant with the vsPrior plot following
                # ax.axvline(priorlo[i] + 0.5*priorspan[i], color='grey', zorder=1)
                # ax.axvline(priorlo[i] + 0.16*priorspan[i], color='grey', zorder=1, ls='--')
                # ax.axvline(priorlo[i] + 0.84*priorspan[i], color='grey', zorder=1, ls='--')

                # darken the lines for the fit results
                # ax.axvline(mcmcMedian[i], color='k', lw=2, zorder=2)
                # ax.axvline(lo[i], color='k', lw=2, zorder=2, ls='--')
                # ax.axvline(hi[i], color='k', lw=2, zorder=2, ls='--')
                # actually the fit result lines are o.k. as is, except all three are dashed
                #  make the median fit a solid line
                #  and maybe change the color to match the central panels
                #  hmm, it's not covering up the dashed line; increase lw and maybe zorder
                # ax.axvline(mcmcMedian[i], color=fitcolor, lw=2, zorder=12)
                ax.axvline(
                    paramValues_bestFit[i], color=fitcolor, lw=2, zorder=12
                )

    if savetodisk:
        plt.savefig(
            saveDir
            + 'corner_'
            + filt
            + '_'
            + modelName
            + '_'
            + trgt
            + ' '
            + p
            + '.png'
        )

    if verbose:
        plt.show()

    return save_plot_tosv(figure), figure


# --------------------------------------------------------------------
def plot_vs_prior(
    allkeys,
    alltraces,
    profiledtraces,
    truth_params,
    prior_ranges,
    appliedLimits,
    filt,
    modelName,
    trgt,
    p,
    saveDir='./',
    savetodisk=False,
):
    '''compare the fit results against the original prior information'''

    mcmcMedian = np.nanmedian(np.array(alltraces), axis=1)
    lo = np.nanpercentile(np.array(alltraces), 16, axis=1)
    hi = np.nanpercentile(np.array(alltraces), 84, axis=1)
    span = hi - lo
    priorlo = np.nanmin(np.array(alltraces), axis=1)
    priorhi = np.nanmax(np.array(alltraces), axis=1)

    mcmcMedianProfiled = np.nanmedian(np.array(profiledtraces), axis=1)
    loProfiled = np.nanpercentile(np.array(profiledtraces), 16, axis=1)
    hiProfiled = np.nanpercentile(np.array(profiledtraces), 84, axis=1)
    spanProfiled = hiProfiled - loProfiled
    priorloProfiled = np.nanmin(np.array(profiledtraces), axis=1)
    priorhiProfiled = np.nanmax(np.array(profiledtraces), axis=1)

    for ikey, key in enumerate(allkeys):
        # print('param:',key)
        # print(' old prior range:',priorlo[ikey],priorhi[ikey])
        if key in prior_ranges.keys():
            priorlo[ikey], priorhi[ikey] = prior_ranges[key]
            # *** this needs work; use appliedLimits somehow ***
            priorloProfiled[ikey], priorhiProfiled[ikey] = prior_ranges[key]
        # print(' new prior range:',priorlo[ikey],priorhi[ikey])
    priorspan = priorhi - priorlo
    priormid = (priorhi + priorlo) / 2.0
    priorspanProfiled = priorhiProfiled - priorloProfiled
    priormidProfiled = (priorhiProfiled + priorloProfiled) / 2.0

    figure = plt.figure(figsize=(12, 6))
    Nparam = len(mcmcMedian)
    for iparam in range(Nparam):
        ax = figure.add_subplot(
            2, int((len(mcmcMedian) + 1.0) / 2.0), iparam + 1
        )
        ax.scatter(
            priormid[iparam],
            priorspan[iparam] * 0.34,
            facecolor='None',
            edgecolor='black',
            s=30,
            zorder=3,
        )
        ax.scatter(
            priormidProfiled[iparam],
            priorspanProfiled[iparam] * 0.34,
            facecolor='black',
            edgecolor='black',
            s=30,
            zorder=3,
        )
        ax.scatter(
            mcmcMedian[iparam],
            span[iparam] / 2.0,
            facecolor='None',
            edgecolor='firebrick',
            s=50,
            zorder=4,
        )
        ax.scatter(
            mcmcMedianProfiled[iparam],
            spanProfiled[iparam] / 2.0,
            facecolor='purple',
            edgecolor='firebrick',
            s=50,
            zorder=4,
        )
        ax.plot(
            [priorlo[iparam], priormid[iparam], priorhi[iparam]],
            [0, priorspan[iparam] * 0.34, 0],
            c='k',
            ls='--',
            lw=0.5,
            zorder=2,
        )
        # ax.plot([priorloProfiled[iparam],priormidProfiled[iparam],priorhiProfiled[iparam]],
        #        [0,priorspanProfiled[iparam]*0.34,0],
        #        c='purple',ls=':',lw=1.5,zorder=2)

        # highlight the background of successful fits (x2 better precision)
        if span[iparam] / 2.0 < priorspan[iparam] * 0.34 / 2:
            # print(' this one has reduced uncertainty',allkeys[iparam])
            ax.set_facecolor('lightyellow')

        # add a dashed line for the true value
        if truth_params is not None:
            keyMatch = {'T': 'Teq', '[X/H]': 'metallicity', '[C/O]': 'C/O'}
            if allkeys[iparam] in keyMatch:
                truthparam = keyMatch[allkeys[iparam]]
            else:
                truthparam = allkeys[iparam]
            if truthparam in truth_params:
                truthvalue = float(truth_params[truthparam])
                ax.plot(
                    [truthvalue, truthvalue],
                    [0, priorspan[iparam]],
                    c='k',
                    ls='--',
                    zorder=5,
                )
            elif truthparam == '[N/O]':
                ax.plot(
                    [0, 0], [0, priorspan[iparam]], c='k', ls='--', zorder=5
                )

        # show if there is some profiling for this parameter
        profiledParams = [limit[0] for limit in appliedLimits]
        profiledValues = [limit[1] for limit in appliedLimits]
        if allkeys[iparam] in profiledParams:
            # print('profiling limit!',allkeys[iparam])
            iprof = profiledParams.index(allkeys[iparam])
            # print('  profiled value',profiledValues[iprof],iprof)
            ax.plot(
                [profiledValues[iprof], profiledValues[iprof]],
                [0, priorspan[iparam]],
                c='purple',
                ls=':',
                lw=2,
                zorder=4,
            )

        ax.set_xlim(priorlo[iparam], priorhi[iparam])
        ax.set_ylim(0, priorspan[iparam] * 0.4)
        ax.set_xlabel(allkeys[iparam], fontsize=14)
        ax.set_ylabel('uncertainty', fontsize=14)
    figure.tight_layout()
    if savetodisk:
        plt.savefig(
            saveDir
            + 'vsprior_'
            + filt
            + '_'
            + modelName
            + '_'
            + trgt
            + ' '
            + p
            + '.png'
        )

    return save_plot_tosv(figure), figure


# --------------------------------------------------------------------
def plot_walker_evolution(
    allkeys,
    alltraces,
    profiledtraces,
    truth_params,
    prior_ranges,
    appliedLimits,
    filt,
    modelName,
    trgt,
    p,
    saveDir='./',
    savetodisk=False,
    Nchains=4,
    verbose=False,
):
    '''trace whether or not the MCMC walkers converge'''

    mcmcMedian = np.nanmedian(np.array(alltraces), axis=1)
    # mcmcMedianProfiled = np.nanmedian(np.array(profiledtraces), axis=1)
    Nparam = len(mcmcMedian)
    priorlo = np.nanmin(np.array(alltraces), axis=1)
    priorhi = np.nanmax(np.array(alltraces), axis=1)
    priorloProfiled = np.nanmin(np.array(profiledtraces), axis=1)
    priorhiProfiled = np.nanmax(np.array(profiledtraces), axis=1)
    for ikey, key in enumerate(allkeys):
        # print('param:',key)
        # print(' old prior range:',priorlo[ikey],priorhi[ikey])
        if key in prior_ranges.keys():
            priorlo[ikey], priorhi[ikey] = prior_ranges[key]
            priorloProfiled[ikey], priorhiProfiled[ikey] = prior_ranges[key]
        # print(' new prior range:',priorlo[ikey],priorhi[ikey]

    figure = plt.figure(figsize=(12, 6))
    linecolors = [
        'crimson',
        'seagreen',
        'royalblue',
        'darkorange',
        'coral',
        'deepskyblue',
        'gold',
        'blueviolet',
    ]
    chainLength = int(len(alltraces[0]) / Nchains)
    # chainLengthProfiled = int(len(profiledtraces[0]) / Nchains)
    for iparam in range(Nparam):
        ax = figure.add_subplot(2, int((Nparam + 1.0) / 2.0), iparam + 1)
        for ic in range(Nchains):
            # jump = ic * chainLength
            # pymc trace has different format than before. no more jumping
            ax.plot(
                np.arange(chainLength) + 1,
                # alltraces[iparam][jump : jump + chainLength],
                alltraces[iparam][ic::Nchains],
                c=linecolors[ic % len(linecolors)],
                alpha=1.0 - (ic / float(Nchains) / 2.0),
                ls='-',
                lw=0.5,
                zorder=3,
            )
        # add a dashed line for the true value
        if truth_params is not None:
            keyMatch = {'T': 'Teq', '[X/H]': 'metallicity', '[C/O]': 'C/O'}
            if allkeys[iparam] in keyMatch:
                truthparam = keyMatch[allkeys[iparam]]
            else:
                truthparam = allkeys[iparam]
            if truthparam in truth_params:
                truthvalue = float(truth_params[truthparam])
                ax.plot(
                    [0, 2 * chainLength],
                    [truthvalue, truthvalue],
                    c='k',
                    ls='--',
                    zorder=5,
                )
            elif truthparam == '[N/O]':
                ax.plot([0, 2 * chainLength], [0, 0], c='k', ls='--', zorder=5)

        # show if there is some profiling for this parameter
        profiledParams = [limit[0] for limit in appliedLimits]
        profiledValues = [limit[1] for limit in appliedLimits]
        profiledSigns = [limit[2] for limit in appliedLimits]
        # print('appliedLimits',appliedLimits)
        if allkeys[iparam] in profiledParams:
            # print('profiling limit!',allkeys[iparam])
            iprof = profiledParams.index(allkeys[iparam])
            # print('  profiled value',profiledValues[iprof],iprof)
            ax.plot(
                [0, 2 * chainLength],
                [profiledValues[iprof], profiledValues[iprof]],
                c='purple',
                ls=':',
                lw=2,
                zorder=4,
            )
            # add grey shading to the area that's dropped by profiling
            greyscale = 0.6
            if profiledSigns[iprof] == '>':
                ax.axhspan(
                    -1.0e5,
                    profiledValues[iprof],
                    alpha=greyscale,
                    color='k',
                    zorder=7,
                )
            elif profiledSigns[iprof] == '<':
                ax.axhspan(
                    1.0e5,
                    profiledValues[iprof],
                    alpha=greyscale,
                    color='k',
                    zorder=7,
                )
            # else:
            #    print(profiledSigns)
            #    exit('ERROR: limit has to be < or >')

        ax.set_xlim(0, chainLength + 1)
        # ax.set_ylim(priorlo[iparam],priorhi[iparam])
        ax.set_ylim(priorloProfiled[iparam], priorhiProfiled[iparam])
        ax.set_xlabel('MCMC step #', fontsize=14)
        ax.set_ylabel(allkeys[iparam], fontsize=14)
    figure.tight_layout()
    if savetodisk:
        plt.savefig(
            saveDir
            + 'walkerevol_'
            + filt
            + '_'
            + modelName
            + '_'
            + trgt
            + ' '
            + p
            + '.png'
        )

    if verbose:
        plt.show()

    return save_plot_tosv(figure), figure


# --------------------------------------------------------------------
def plot_fits_vs_truths(
    truth_values,
    fit_values,
    fit_errors,
    prior_ranges,
    filt,
    saveDir='./',
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
    fit_values, fit_errors, prior_ranges, filt, saveDir='./', savetodisk=False
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
def plot_mass_vs_metals(
    masses,
    stellarFEHs,
    truth_values,
    fit_values,
    fit_errors,
    prior_ranges,
    filt,
    saveDir='./',
    plot_truths=False,  # for Ariel-sims, include truth as open circles?
    savetodisk=False,
):
    '''how well do we retrieve the input mass-metallicity relation?'''

    onlyFitAbove10MEarth = True
    onlyPlotAbove10MEarth = True

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

        if (
            onlyFitAbove10MEarth
            and onlyPlotAbove10MEarth
            and (mass < 10 * MEarth)
        ):
            pass
        else:
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
    if onlyFitAbove10MEarth and onlyPlotAbove10MEarth:
        xrange = (0.01, 10.0)

    # plot the underlying distribution (only if this is a simulation)
    # actually, also plot Thorngren relationship for real HST data
    # and also include the newer Chachan relationship
    massesAnalytic = np.logspace(-5, 3, 100)
    metalsThorngren = massMetalRelation(0, massesAnalytic, thorngren=True)
    metalsChachan = massMetalRelation(0, massesAnalytic, chachan=True)
    # if 'sim' in filt:
    #     label='true relationship')
    ax.plot(
        massesAnalytic,
        metalsThorngren,
        'k:',
        lw=1,
        zorder=1,
        label='Thorngren+ 2016',
    )
    ax.plot(
        massesAnalytic,
        metalsChachan,
        'k--',
        lw=1,
        zorder=1,
        label='Chachan+ 2025',
    )

    # some trouble with a 0+-0 point. make sure uncertainty is not zero
    metals_fiterr = np.array(metals_fiterr)
    metals_fiterr[np.where(metals_fiterr <= 0)] = 666

    # plot a linear fit to the data
    if onlyFitAbove10MEarth:
        limitedMasses = np.where(np.array(masses) > 10.0 * MEarth)
        masses_noSmallOnes = np.array(masses)[limitedMasses]
        metals_noSmallOnes = np.array(metals_fit)[limitedMasses]
        metalserr_noSmallOnes = np.array(metals_fiterr)[limitedMasses]
        # print('# of planet masses', len(masses))
        # print('# of planet masses >10MEarth', len(masses_noSmallOnes))
        polynomialCoeffs, covariance = np.polyfit(
            np.log10(masses_noSmallOnes),
            metals_noSmallOnes,
            1,
            w=1.0 / np.array(metalserr_noSmallOnes),
            cov='unscaled',
        )
    else:
        polynomialCoeffs, covariance = np.polyfit(
            np.log10(masses),
            metals_fit,
            1,
            w=1.0 / np.array(metals_fiterr),
            cov='unscaled',
        )
    # print('polynomialCoeffs',polynomialCoeffs)
    lineFunction = np.poly1d(polynomialCoeffs)
    massrange = np.linspace(xrange[0], xrange[1], 10)
    plt.plot(
        massrange,
        lineFunction(np.log10(massrange)),
        'k--',
        lw=1.5,
        zorder=3,
        label='linear fit',
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
        yrange[0] + 0.2 + 0.06 * (yrange[1] - yrange[0]),
        resultsstring,
        c='black',
        ha='left',
        fontsize=10,
    )
    plt.text(
        xrange[0] * 1.2,
        yrange[0] + 0.2,
        'Z$_p$ = (9.7$\\pm$1.3) $M_p^{-0.45\\pm0.09}$ (Thorngren)',
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
        if (
            onlyFitAbove10MEarth
            and onlyPlotAbove10MEarth
            and (mass < 10 * MEarth)
        ):
            pass
        else:
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
    if onlyFitAbove10MEarth and onlyPlotAbove10MEarth:
        xrange = (0.01, 10.0)

    # plot the underlying distribution (only if this is a simulation)
    # actually, also plot Thorngren relationship for real HST data
    ax2.plot(
        massesAnalytic,
        metalsThorngren,
        'k:',
        lw=1,
        zorder=1,
        label='Thorngren+ 2016',
    )
    ax2.plot(
        massesAnalytic,
        metalsChachan,
        'k--',
        lw=1,
        zorder=1,
        label='Chachan+ 2025',
    )

    # plot a linear fit to the data
    if onlyFitAbove10MEarth:
        MEarth = 5.972e27 / 1.898e30
        limitedMasses = np.where(np.array(masses) > 10.0 * MEarth)
        masses_noSmallOnes = np.array(masses)[limitedMasses]
        metals_noSmallOnes = np.array(metals_fit)[limitedMasses]
        metalserr_noSmallOnes = np.array(metals_fiterr)[limitedMasses]
        stellarFEHs_noSmallOnes = np.array(stellarFEHs)[limitedMasses]
        # print('# of planet masses',len(masses))
        # print('# of planet masses >10MEarth',len(masses_noSmallOnes))
        polynomialCoeffs, covariance = np.polyfit(
            np.log10(masses_noSmallOnes),
            metals_noSmallOnes - stellarFEHs_noSmallOnes,
            1,
            w=1.0 / np.array(metalserr_noSmallOnes),
            cov='unscaled',
        )
    else:
        polynomialCoeffs, covariance = np.polyfit(
            np.log10(masses),
            np.array(metals_fit) - np.array(stellarFEHs),
            1,
            w=1.0 / np.array(metals_fiterr),
            cov='unscaled',
        )
    # print('polynomialCoeffs',polynomialCoeffs)
    # print(' covariance',covariance)
    lineFunction = np.poly1d(polynomialCoeffs)
    massrange = np.linspace(xrange[0], xrange[1], 10)
    plt.plot(
        massrange,
        lineFunction(np.log10(massrange)),
        'k--',
        lw=1.5,
        zorder=3,
        label='linear fit',
    )

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
        yrange[0] + 0.2 + 0.06 * (yrange[1] - yrange[0]),
        resultsstring,
        c='black',
        ha='left',
        fontsize=10,
    )
    plt.text(
        xrange[0] * 1.2,
        yrange[0] + 0.2,
        'Z$_p$ = (9.7$\\pm$1.3) $M_p^{-0.45\\pm0.09}$ (Thorngren)',
        c='black',
        ha='left',
        fontsize=10,
    )

    plt.legend()

    ax2.set_xlim(xrange)
    ax2.set_ylim(yrange)
    figure.tight_layout()  # some trouble with this, with bigger figure maybe?

    # ('display' doesn't work for pdf files)
    if savetodisk:
        plt.savefig(saveDir + 'massVSmetals_' + filt + '.png')
    return save_plot_tosv(figure), figure


# --------------------------------------------------------------------
