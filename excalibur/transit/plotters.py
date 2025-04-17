'''transit plotting ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import os
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


def plot_corner(
    allkeys,
    alltraces,
    modelParams_bestFit,
    prior_ranges,
    filt,
    modelName,
    trgt,
    p,
    saveDir=os.path.join(excalibur.context['data_dir'], 'bryden/'),
    savetodisk=False,
):
    '''corner plot showing posterior distributions'''

    fitcolor = 'firebrick'

    tpr, ctp, hza, hloc, hthc, tceqdict, mixratio = modelParams_bestFit
    # print('model param in corner plot',modelParams_bestFit)

    paramValues_bestFit = []
    for param in allkeys:
        if param == 'T':
            paramValues_bestFit.append(tpr)
        elif param == 'CTP':
            paramValues_bestFit.append(ctp)
        elif param == 'HScale':
            paramValues_bestFit.append(hza)
        elif param == 'HLoc':
            paramValues_bestFit.append(hloc)
        elif param == 'HThick':
            paramValues_bestFit.append(hthc)
        elif param == '[X/H]':
            paramValues_bestFit.append(tceqdict['XtoH'])
        elif param == '[C/O]':
            paramValues_bestFit.append(tceqdict['CtoO'])
        elif param == '[N/O]':
            paramValues_bestFit.append(tceqdict['NtoO'])
        elif param in mixratio:
            paramValues_bestFit.append(mixratio[param])
        else:
            log.warning('--< ERROR: param not in list: %s >--', param)

    # print('best fit values in corner plot',paramValues_bestFit)

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
    figure = corner.corner(
        np.vstack(np.array(alltraces)).T,
        bins=10,
        labels=allkeys,
        range=trange,
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
    #  draw a point and crosshair for the medians in each subpanel
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
        ax.axvline(paramValues_bestFit[i], color=fitcolor, lw=2, zorder=12)

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

    return save_plot_tosv(figure), figure


# --------------------------------------------------------------------
