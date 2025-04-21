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

import excalibur
from excalibur.util.plotters import save_plot_tosv

# from excalibur.system.core import ssconstants
# ssc = ssconstants(mks=True)

log = logging.getLogger(__name__)


# --------------------------------------------------------------------


def plot_corner(
    # allkeys,
    alltraces,
    # modelParams_bestFit,
    prior_ranges,
    p,
    saveDir=os.path.join(excalibur.context['data_dir'], 'bryden/'),
    savetodisk=False,
):
    '''corner plot showing posterior distributions (transit, not cerberus)'''

    fitcolor = 'firebrick'

    allkeys = []
    modelParams_bestFit = []
    for key, values in alltraces.items():
        allkeys.append(key)
        modelParams_bestFit.append(np.nanmedian(values))
        print(
            'WHITELIGHT trace median,std,min,max',
            key,
            np.nanmedian(values),
            np.nanstd(values),
            np.nanmin(values),
            np.nanmax(values),
        )

    print(' params inside of corner plotting', allkeys)

    print('mctrace params', alltraces.keys())

    print('bestfit values passed in (currently median)', modelParams_bestFit)
    print()

    lo = []
    hi = []
    priorlo = []
    priorhi = []
    mcmcMedian = []
    # awkward because alltraces is a dictionary
    #  convert to 2-d array first (needed for corner() anyway)
    paramnames = []
    tracearray = []
    for param, values in alltraces.items():
        paramnames.append(param)
        tracearray.append(values)
        print('shape added', param, values.shape)
    print('tracearray', tracearray)
    tracearray = np.array(tracearray)

    for param, values in alltraces.items():
        mcmcMedian.append(np.nanmedian(values))
        lo.append(np.nanpercentile(np.array(values), 16))
        hi.append(np.nanpercentile(np.array(values), 84))
        # span = hi - lo
        # Careful! these are not actually the prior ranges;
        #  they're the range of walker values (unless set below)
        priorlo.append(np.nanmin(values))
        priorhi.append(np.nanmax(values))
    print('medians inside of corner plotting', mcmcMedian)

    mcmcMedian = np.nanmedian(tracearray, axis=1)
    lo = np.nanpercentile(np.array(tracearray), 16, axis=1)
    hi = np.nanpercentile(np.array(tracearray), 84, axis=1)
    # span = hi - lo
    # Careful! these are not actually the prior ranges;
    #  they're the range of walker values (unless set below)
    priorlo = np.nanmin(tracearray, axis=1)
    priorhi = np.nanmax(tracearray, axis=1)
    print('medians inside of corner plotting', mcmcMedian)
    print()
    print()

    for ikey, key in enumerate(allkeys):
        print('param:', ikey, key)
        print('  old prior range:', key, priorlo[ikey], priorhi[ikey])
        if key in prior_ranges.keys():
            priorlo[ikey] = prior_ranges[key][0]
            priorhi[ikey] = prior_ranges[key][1]
        else:
            print('  TROUBLE: ', key, 'not found in', prior_ranges.keys())
        print('  new prior range:', key, priorlo[ikey], priorhi[ikey])
    # priorspan = priorhi - priorlo
    # priormid = (priorhi + priorlo) / 2.

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
        np.vstack(np.array(tracearray)).T,
        bins=10,
        labels=allkeys,
        range=trange,
        show_titles=True,
        quantiles=[0.16, 0.50, 0.84],
    )

    ndim = len(allkeys)
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
            ax.axvline(modelParams_bestFit[xi], color=fitcolor)
            ax.axhline(modelParams_bestFit[yi], color=fitcolor)
            ax.plot(
                modelParams_bestFit[xi],
                modelParams_bestFit[yi],
                marker='s',
                c=fitcolor,
            )

    # darken the lines for the fit results
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(modelParams_bestFit[i], color=fitcolor, lw=2, zorder=12)

    if savetodisk:
        plt.savefig(saveDir + 'corner-' + p + '.png')

    return save_plot_tosv(figure), figure


# --------------------------------------------------------------------
