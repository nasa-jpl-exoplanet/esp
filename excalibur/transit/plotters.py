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
    alltraces,
    prior_ranges,
    p,
    saveDir=os.path.join(excalibur.context['data_dir'], 'debug/'),
    savetodisk=False,
):
    '''corner plot showing posterior distributions (transit, not cerberus)'''

    fitcolor = 'firebrick'

    fit_param_names = []
    modelParams_bestFit = []
    chainlen = 0
    numwalkers = 0    
    for key, values in alltraces.items():
        fit_param_names.append(key)
        modelParams_bestFit.append(np.nanmedian(values))
        chainlen = np.max([chainlen, values.shape[1]])
        numwalkers = np.max([numwalkers, values.shape[0]])
        # print(
        #    'WHITELIGHT trace median,std,min,max',
        #    key,
        #    np.nanmedian(values),
        #    np.nanstd(values),
        #    np.nanmin(values),
        #    np.nanmax(values),
        # )
    # print('chainlen,numwalkers', chainlen, numwalkers)

    # print('mctrace params', alltraces.keys())
    # print(' params inside of corner plotting', fit_param_names)
    # place rprs as the first parameters? asdf

    # print('bestfit values passed in (currently median)', modelParams_bestFit)
    # print()

    #  convert to the tracearray to a 2-d array (corner() format)
    tracearray = []
    for _, values in alltraces.items():
        # fixed parameters only have 2 steps. extend them
        if len(values[0]) < chainlen:
            # print('EXTENDING A FIXED PARAM', param, numwalkers, chainlen)
            tracearray.append(
                float(values[0][0]) * np.ones((numwalkers, chainlen))
            )
        else:
            tracearray.append(np.array(values))
    tracearray = np.array(tracearray)

    # Careful! these are not actually the prior ranges;
    #  they're the range of walker values (unless set below)
    priorlo = np.nanmin(tracearray, axis=(1, 2))
    priorhi = np.nanmax(tracearray, axis=(1, 2))
    # print('  priorlo', priorlo)

    # for cases with fixed params, make sure the plots have some range
    #  actually this is not needed. corner() can handle it

    if prior_ranges is not None:
        for ikey, key in enumerate(fit_param_names):
            if key in prior_ranges.keys():
                # print(
                #    '  old prior range (median):', key, priorlo[ikey], priorhi[ikey]
                # )
                priorlo[ikey] = prior_ranges[key][0]
                priorhi[ikey] = prior_ranges[key][1]
                # print(
                #    '  new prior range (true)  :', key, priorlo[ikey], priorhi[ikey]
                # )
                pass
            pass
        pass
    trange = [tuple([x, y]) for x, y in zip(priorlo, priorhi)]

    figure = corner.corner(
        np.vstack(np.array(tracearray).T),
        bins=10,
        labels=fit_param_names,
        range=trange,
        show_titles=True,
        quantiles=[0.16, 0.50, 0.84],
    )

    ndim = len(fit_param_names)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    # use larger font size for the axis labels
    for i in range(ndim):
        ax = axes[ndim - 1, i]
        ax.set_xlabel(fit_param_names[i], fontsize=14)
    for i in range(ndim - 1):
        # skipping the first one on the y side (it's a histo, not a 2-D plot)
        ax = axes[i + 1, 0]
        ax.set_ylabel(fit_param_names[i + 1], fontsize=14)
    #  draw a point and crosshair for the medians in each subpanel
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
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
        saveDir = '/proj/sdp/data/debug/'
        # print('  Saving figure to disk!',saveDir+'corner-'+p+'.png')
        plt.savefig(saveDir + 'corner-' + p + '.png')

    return save_plot_tosv(figure)


# --------------------------------------------------------------------
