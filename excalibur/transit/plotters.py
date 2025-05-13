'''transit plotting ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import os
import pymc
import numpy as np
import corner
import logging
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import excalibur
from excalibur.util.plotters import save_plot_tosv

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
    for key, values in alltraces.items():
        fit_param_names.append(key)
        modelParams_bestFit.append(np.nanmedian(values))
        # print(
        #    'WHITELIGHT trace median,std,min,max',
        #    key,
        #    np.nanmedian(values),
        #    np.nanstd(values),
        #    np.nanmin(values),
        #    np.nanmax(values),
        # )
    # print('mctrace params', alltraces.keys())
    # print(' params inside of corner plotting', fit_param_names)
    # place rprs as the first parameters? asdf

    #  convert to the tracearray to a 2-d array (corner() format)
    tracearray = []
    for _, values in alltraces.items():
        tracearray.append(np.array(values).flatten())
    tracearray = np.array(tracearray)

    # Careful! priorlo,priorhi are not actually the prior ranges;
    #  they're the range of walker values (unless set below)
    priorlo = np.nanmin(tracearray, axis=1)
    priorhi = np.nanmax(tracearray, axis=1)
    # print('  priorlo', priorlo)
    # print('  priorhi', priorhi)
    # print()
    if prior_ranges is not None:
        for ikey, key in enumerate(fit_param_names):
            if key in prior_ranges.keys():
                priorlo[ikey] = np.min((prior_ranges[key][0], priorlo[ikey]))
                priorhi[ikey] = np.max((prior_ranges[key][1], priorhi[ikey]))
                pass
            pass
        pass
    # print('  priorlo', priorlo)
    # print('  priorhi', priorhi)
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

    # Corner.corner creates a figure instance.
    # Better close it before memory leaks out
    plt.close()

    return save_plot_tosv(figure)


def simplecorner(mctrace, verbose=False):
    '''
    GMR: Simple corner plot
    '''
    Medians = [float(np.median(mctrace[k])) for k in mctrace]
    # Lowerr = [float(np.percentile(mctrace[k], 50 - 95.4 / 2)) for k in mctrace]
    # Uperr = [float(np.percentile(mctrace[k], 50 + 95.4 / 2)) for k in mctrace]
    LowBound = [
        float(np.percentile(mctrace[k], 50 - 99.7 / 2)) for k in mctrace
    ]
    UpBound = [float(np.percentile(mctrace[k], 50 + 99.7 / 2)) for k in mctrace]
    Ranges = list(zip(LowBound, UpBound))
    Figure = corner.corner(
        mctrace,
        quantiles=(0.5 - 0.68 / 2, 0.5 + 0.68 / 2),
        levels=(
            0.393,
            0.675,
        ),
        range=Ranges,
    )
    corner.overplot_lines(Figure, Medians, c='b')
    corner.overplot_points(Figure, np.array(Medians)[None], marker='o', c='b')
    out = save_plot_tosv(Figure)
    if verbose:
        plt.show()
        pass
    else:
        plt.close()
        pass
    return out


def postpriors(mctrace, prior_center, nodes, verbose=False):
    '''
    GMR: Plots post versus prior per parameter
    '''
    Medians = [float(np.median(mctrace[k])) for k in mctrace]
    Figure, Axes = plt.subplots(
        nrows=len(mctrace), ncols=2, figsize=(4 * 2, 3 * len(mctrace))
    )
    for index, k in enumerate(mctrace):
        Axes[index, 0].set_title(k)
        Axes[index, 0].plot(mctrace[k].flatten(), alpha=0.7)
        Axes[index, 0].axhline(Medians[index], color='k')
        Axes[index, 0].axhline(prior_center[k], ls='--', color='lightgray')
        Axes[index, 1].hist(mctrace[k].flatten())
        Axes[index, 1].hist(
            pymc.draw(nodes[index], draws=len(mctrace[k].flatten())),
            color='lightgray',
            alpha=0.7,
        )
        pass
    out = save_plot_tosv(Figure)
    if verbose:
        plt.show()
        pass
    else:
        plt.close()
        pass
    return out


def lightcurves(SvDataP, p, mergesv=False, verbose=False):
    visits = SvDataP['visits']
    # phase,allwhite is the data before shifting
    phase = SvDataP['phase']
    allwhite = SvDataP['allwhite']
    # postphase,allwhite/postim is the data after shifting
    postphase = SvDataP['postphase']
    postim = SvDataP['postim']
    # phase,postlc is the model
    postflatphase = np.array(SvDataP['postflatphase'])
    postlc = np.array(SvDataP['postlc'])
    # myfig = plt.figure(figsize=(10, 6))
    # plt.title('Planet '+p, fontsize=14)
    myfig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1]},
    )
    myfig.subplots_adjust(hspace=0.03, right=0.8)
    ax1.set_title('Planet ' + p, fontsize=16)
    # ax1.set_xlabel('Orbital Phase', fontsize=14)
    ax1.set_ylabel('Normalized Post-Whitelight Curve', fontsize=14)
    for index, v in enumerate(visits):
        # plot the normalized/shifted data
        if mergesv:
            vlabel = SvDataP['allfltrs'][index]
        else:
            vlabel = 'visit ' + str(v)
            pass
        ax1.plot(
            np.array(postphase[index]),
            np.array(allwhite[index]) / np.array(postim[index]),
            'o',
            zorder=3,
            label=vlabel,
        )
        # plot the pre-correction data
        if index == len(visits) - 1:
            ax1.plot(
                np.array(phase[index]),
                np.array(allwhite[index]),
                'k+',
                zorder=2,
                label='pre-correction',
            )
        else:
            ax1.plot(
                np.array(phase[index]),
                np.array(allwhite[index]),
                'k+',
                zorder=2,
            )
            pass
        # add a lower panel showing the data-model residuals
        model_interpolator = interp1d(
            postflatphase[np.argsort(postflatphase)],
            postlc[np.argsort(postflatphase)],
            kind='linear',
            fill_value='extrapolate',
        )
        model_at_observed_time = model_interpolator(np.array(postphase[index]))
        residuals = (
            np.array(allwhite[index]) / np.array(postim[index])
            - model_at_observed_time
        )

        ax2.plot(
            np.array(postphase[index]),
            residuals,
            'o',
            color='black',
            markerfacecolor='None',
        )
        ax2.axhline(y=0, color='black', linestyle='dashed', linewidth=1)
        ax2.set_xlabel('Orbital Phase', fontsize=14)
        ax2.set_ylabel('Residuals', fontsize=14)
        pass

    xdatarange = ax1.get_xlim()
    # use full model (not just at some points) if available
    if 'modellc' in SvDataP.keys():
        modelphase = np.array(SvDataP['modelphase'])
        modellc = np.array(SvDataP['modellc'])
        # the points are ordered by time, not by phase
        #  so sorting is needed for multi-visit observations
        # otherwise you get weird wrap-around in the line plots
        ax1.plot(
            modelphase[np.argsort(modelphase)],
            modellc[np.argsort(modelphase)],
            '-',
            c='k',
            marker='None',
            zorder=1,
            label='model',
        )
        # model phases only go from -0.5 to 0.5 (not good for eclipse)
        # plot the model line a second time, but shifting the phases over by 1
        ax1.plot(
            modelphase[np.argsort(modelphase)] + 1,
            modellc[np.argsort(modelphase)],
            '-',
            c='k',
            marker='None',
            zorder=1,
        )
    else:
        ax1.plot(postflatphase, postlc, '^', zorder=1, label='model')
        # '^', c='green', zorder=1, label='model')
        pass
    ax1.set_xlim(xdatarange)
    if len(visits) > 14:
        ncol = 2
        pass
    else:
        ncol = 1
        if mergesv:
            ax1.legend(
                loc='best',
                ncol=ncol,
                mode='expand',
                numpoints=1,
                borderaxespad=0.0,
                frameon=False,
            )
            # myfig.tight_layout()
        else:
            ax1.legend(
                bbox_to_anchor=(1 + 0.1 * (ncol - 0.5), 0.5),
                loc=5,
                ncol=ncol,
                mode='expand',
                numpoints=1,
                borderaxespad=0.0,
                frameon=False,
            )
            # myfig.tight_layout(rect=[0,0,(1 - 0.1*ncol),1])
            pass
        pass
    out = save_plot_tosv(myfig)
    if verbose:
        plt.show()
        pass
    else:
        plt.close()
        pass
    return out


# --------------------------------------------------------------------
