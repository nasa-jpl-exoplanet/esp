'''phasecurve plotters ds'''

# Heritage code shame:
# asdfasdf pylint: disable=duplicate-code
# asdfadsf  disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements

# -- IMPORTS ---------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import logging

from excalibur.util import elca

log = logging.getLogger(__name__)


# --------------------------------------------------------------------


def plot_phasecurve(
    whitelightdata,
    # systemparam,
    flarephases,
    bin_dt=10.0 / (60 * 24),
    zoom=False,
    phase=True,
    verbose=False,
):
    '''
    plot a cleaned phasecurve.  show where there are flares
    '''
    # print('keys in the input data',whitelightdata.keys())
    # print('finalpar',whitelightdata['final_pars'])

    # ok the final_param numbers are not the same as sys param
    # use the fit values, not the default system params
    # print('period check',
    #      whitelightdata['final_pars']['per'],
    #      systemparam['period'])
    # print('tmid check',
    #      whitelightdata['final_pars']['tmid'],
    #      systemparam['t0'])

    f = plt.figure(figsize=(10, 5))
    ax_lc = plt.subplot2grid((4, 5), (0, 0), colspan=5, rowspan=3)
    ax_res = plt.subplot2grid((4, 5), (3, 0), colspan=5, rowspan=1)
    axs = [ax_lc, ax_res]

    bt, bf, _ = elca.time_bin(
        whitelightdata['time'], whitelightdata['detrended'], bin_dt
    )
    bp = (bt - whitelightdata['final_pars']['tmid']) / whitelightdata[
        'final_pars'
    ]['per']

    if phase:
        axs[0].plot(bp, bf, 'co', alpha=0.5, zorder=2)
        axs[0].plot(
            whitelightdata['phase'], whitelightdata['transit'], 'r-', zorder=3
        )
        axs[0].set_xlim(
            [min(whitelightdata['phase']), max(whitelightdata['phase'])]
        )
        axs[0].set_xlabel('Phase ')
    else:
        axs[0].plot(bt, bf, 'co', alpha=0.5, zorder=2)
        axs[0].plot(
            whitelightdata['time'], whitelightdata['transit'], 'r-', zorder=3
        )
        axs[0].set_xlim(
            [min(whitelightdata['time']), max(whitelightdata['time'])]
        )
        axs[0].set_xlabel('Time [day]')

    axs[0].set_ylabel('Relative Flux')
    axs[0].grid(True, ls='--')

    if zoom:
        axs[0].set_ylim(
            [
                1 - 1.25 * whitelightdata['final_pars']['rprs'] ** 2,
                1 + 0.5 * whitelightdata['final_pars']['rprs'] ** 2,
            ]
        )
    else:
        if phase:
            axs[0].errorbar(
                whitelightdata['phase'],
                whitelightdata['detrended'],
                yerr=np.std(whitelightdata['residuals'])
                / np.median(whitelightdata['flux']),
                ls='none',
                marker='.',
                color='black',
                zorder=1,
                alpha=0.025,
            )
        else:
            axs[0].errorbar(
                whitelightdata['time'],
                whitelightdata['detrended'],
                yerr=np.std(whitelightdata['residuals'])
                / np.median(whitelightdata['flux']),
                ls='none',
                marker='.',
                color='black',
                zorder=1,
                alpha=0.025,
            )

    # ____________ show when the flares happen ____________
    flareplotloc = np.median(whitelightdata['detrended']) + 3.0 * np.std(
        whitelightdata['detrended']
    )
    axs[0].plot(
        flarephases, [flareplotloc] * len(flarephases), 'rv', label='flares'
    )
    axs[0].legend(loc='center left', bbox_to_anchor=(1.06, 0.48))

    bt, br, _ = elca.time_bin(
        whitelightdata['time'],
        whitelightdata['residuals'] / np.median(whitelightdata['flux']) * 1e6,
        bin_dt,
    )
    bp = (bt - whitelightdata['final_pars']['tmid']) / whitelightdata[
        'final_pars'
    ]['per']

    if phase:
        axs[1].plot(
            whitelightdata['phase'],
            whitelightdata['residuals']
            / np.median(whitelightdata['flux'])
            * 1e6,
            'k.',
            alpha=0.15,
            label=fr'$\sigma$ = {np.std(whitelightdata['residuals']
                                        / np.median(whitelightdata['flux']
                                                    ) * 1e6):.0f} ppm',
        )
        axs[1].plot(
            bp,
            br,
            'c.',
            alpha=0.5,
            zorder=2,
            label=fr'$\sigma$ = {np.std(br):.0f} ppm',
        )
        axs[1].set_xlim(
            [min(whitelightdata['phase']), max(whitelightdata['phase'])]
        )
        axs[1].set_xlabel('Phase')
    else:
        axs[1].plot(
            whitelightdata['time'],
            whitelightdata['residuals']
            / np.median(whitelightdata['flux'])
            * 1e6,
            'k.',
            alpha=0.15,
            label=fr'$\sigma$ = {np.std(whitelightdata['residuals']
                                        / np.median(whitelightdata['flux'])
                                        * 1e6):.0f} ppm',
        )
        axs[1].plot(
            bt,
            br,
            'c.',
            alpha=0.5,
            zorder=2,
            label=fr'$\sigma$ = {np.std(br):.0f} ppm',
        )
        axs[1].set_xlim(
            [min(whitelightdata['time']), max(whitelightdata['time'])]
        )
        axs[1].set_xlabel('Time [day]')

    axs[1].legend(loc='center left', bbox_to_anchor=(1.06, 0.48))
    axs[1].set_ylabel('Residuals [ppm]')
    axs[1].grid(True, ls='--')
    plt.tight_layout()

    if verbose:
        plt.show()

    return f, axs
