'''cerberus plotting ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import logging
import numpy as np
import matplotlib.pyplot as plt

from excalibur.system.core import ssconstants
from excalibur.util.plotters import save_plot_tosv

log = logging.getLogger(__name__)

ssc = ssconstants(mks=True)


# --------------------------------------------------------------------
def plot_ML_fits_vs_truths(
    param_names,
    input_params,
    MLfit_params,
    saveDir='./',
    savetodisk=False,
    verbose=False,
):
    '''
    compare ML model retrieved values against the input parameters
    (no instrument included here; this is just a measure of ML uncertainty)
    '''

    numColumns = int(np.ceil(len(param_names) / 2.0))
    fig, axs = plt.subplots(2, numColumns, figsize=(3.5 * numColumns, 7))
    fig.subplots_adjust(
        left=0.15, right=0.95, bottom=0.15, top=0.92, wspace=0.4, hspace=0.4
    )
    axs = axs.flatten()

    for k, param_name in enumerate(param_names):
        ax = axs[k]
        ax.set_aspect('equal')
        ax.scatter(input_params[:, k], MLfit_params[:, k], s=5, alpha=0.3)

        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        if verbose:
            print(' prior range for', param_name, '=', xlims)
        vmin = min([xlims[0], ylims[0]])
        vmax = max([xlims[1], ylims[1]])
        # vmin = min(input_params[:, k].min(), MLfit_params[:, k].min())
        # vmax = max(input_params[:, k].max(), MLfit_params[:, k].max())
        # print('vmin,vmax', vmin, vmax)
        ax.plot([vmin, vmax], [vmin, vmax], 'k--', lw=1)

        ax.set_xlim((vmin, vmax))
        ax.set_ylim((vmin, vmax))

        # ax.set_xlabel(f'True {param_name}', fontsize=14)
        # ax.set_ylabel(f'Predicted {param_name}', fontsize=14)
        ax.set_xlabel(f'Input value', fontsize=14)
        ax.set_ylabel(f'Retrieved value', fontsize=14)
        ax.set_title(param_name, fontsize=16)

    # if there's an odd number of parameters, hide the last empty frame
    if len(param_names) % 2 == 1:
        axs[-1].axis("off")

    # plt.tight_layout()
    if savetodisk:
        plt.savefig(saveDir + 'MLfitVStruth.png')
    if verbose:
        plt.show()
    return save_plot_tosv(fig), fig


# --------------------------------------------------------------------
def plot_ML_spectrumfit(
    transitdata,
    truth_spectrum,
    ML_best_fit,
    ML_param_names,
    ML_param_names_forprint,
    ML_param_results,
    system_data,
    ancillary_data,
    filt,
    trgt,
    p,
    saveDir='./',
    savetodisk=False,
):
    '''plot the machine-learning best fit to the data'''

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
    # 3b) plot the best-fit model - new random selection parameter-set checking
    plt.plot(
        transitdata['wavelength'][okPart],
        ML_best_fit * 100,
        # c='k', lw=2, zorder=4,
        c='orange',
        lw=2,
        zorder=4,
        label='best fit',
    )
    # 4) plot the true spectrum, if it is a simulation
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

    xlims = plt.xlim()
    ylims = plt.ylim()

    offsets_model = (ML_best_fit - transitdata['depth'][okPart]) / transitdata[
        'error'
    ][okPart]
    chi2model = np.nansum(offsets_model**2)

    numParam_model = 8
    numParam_truth = 0

    numPoints = len(ML_best_fit)
    # print('numpoints',numPoints)
    chi2model_red = chi2model / (numPoints - numParam_model)

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
            ylims[0] + (ylims[1] - ylims[0]) * 0.9,
            '$\\chi^2$-truth=' + f"{chi2truth:5.2f}",
            fontsize=12,
        )
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * 0.74,
            '$\\chi^2_{red}$-truth=' + f"{chi2truth_red:5.2f}",
            fontsize=12,
        )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.83,
        '$\\chi^2$-model=' + f"{chi2model:5.2f}",
        fontsize=12,
    )
    plt.text(
        xlims[1] + xoffset,
        ylims[0] + (ylims[1] - ylims[0]) * 0.67,
        '$\\chi^2_{red}$-model=' + f"{chi2model_red:5.2f}",
        fontsize=12,
    )

    # print out the best-fit parameters on the right side
    yloc = 0.55
    for param, name in zip(ML_param_names, ML_param_names_forprint):
        # print(param, '=', ML_param_results[param])
        yloc -= 0.08
        if '(' in name:
            namesplit = name.index('(')
            units = name[namesplit:]
            name = name[:namesplit]
        else:
            units = ''
        plt.text(
            xlims[1] + xoffset,
            ylims[0] + (ylims[1] - ylims[0]) * yloc,
            name + '= ' + f"{ML_param_results[param]:5.2f}" + ' ' + units,
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

    # figgy.tight_layout()
    # plt.show()
    if savetodisk:
        # pdf is so much better, but xv gives error (stick with png for debugging)
        plt.savefig(saveDir + 'bestFit_' + filt + '_' + trgt + ' ' + p + '.png')

    return save_plot_tosv(figgy), figgy


# --------------------------------------------------------------------
