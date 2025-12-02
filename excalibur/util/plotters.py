'''util plotters dc'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import io
from scipy.stats import skew, kurtosis
from numpy.fft import fft, fftfreq

import dawgie


# -- PLOTTING FUNCTIONS-- --------------------------------------------
# --------------------------------------------------------------------
def save_plot_toscreen(fig, visitor, headertext=' '):
    # save a plot to the screen
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    visitor.add_image('...', headertext, buf.getvalue())
    plt.close(fig)
    return


def save_plot_tosv(fig):
    # save a plot as a state vector
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    return buf.getvalue()


def save_plot_myfit(plotfn):
    # extract plot data for states.py. incoming data is from pc_fitter
    fig, _ = plotfn()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    plt.close(fig)
    return buf.getvalue()


# --------------------------------------------------------------------
def plot_residual_fft(
    selftype,
    fltr,
    p,
    aper,
    subt,
    myfit,
    tdur_freq=None,
):
    # create plot for residual statistics

    # estimates for photon noise
    photons = aper * 1.0  # already converted to e- in data task

    # noise estimate in transit
    # tmask = myfit.transit < 1  # not used
    photon_noise_timeseries = 1 / np.sqrt(photons.mean())

    # photon noise factor based on timeseries
    res_std = np.round(np.std(myfit.residuals / np.median(aper)), 7)
    nf_timeseries = res_std / photon_noise_timeseries
    raw_residual = aper / np.median(aper) - myfit.transit
    nf_timeseries_raw = np.std(raw_residual) / photon_noise_timeseries
    rel_residuals = myfit.residuals / np.median(aper)

    # print(f"raw photon noise:{nf_timeseries_raw}")
    # print(f"photon noise: {nf_timeseries}")

    fig, ax = plt.subplots(3, figsize=(10, 10))
    binspace = np.linspace(-0.02, 0.02, 201)
    raw_label = (
        f"Mean: {np.mean(raw_residual,):.4f} \n"
        f"Stdev: {np.std(raw_residual):.4f} \n"
        f"Skew: {skew(raw_residual):.4f} \n"
        f"Kurtosis: {kurtosis(raw_residual):.4f}\n"
        f"Photon Noise: {nf_timeseries_raw:.2f}"
    )
    ax[0].hist(
        raw_residual,
        bins=binspace,
        label=raw_label,
        color=plt.cm.jet(0.25),
        alpha=0.5,
    )
    detrend_label = (
        f"Mean: {np.mean(rel_residuals):.4f} \n"
        f"Stdev: {np.std(rel_residuals):.4f} \n"
        f"Skew: {skew(rel_residuals):.4f} \n"
        f"Kurtosis: {kurtosis(rel_residuals):.4f}\n"
        f"Photon Noise: {nf_timeseries:.2f}"
    )
    ax[0].hist(
        rel_residuals,
        bins=binspace,
        label=detrend_label,
        color=plt.cm.jet(0.75),
        alpha=0.5,
    )
    ax[0].set_xlabel('Relative Flux Residuals')
    ax[0].legend(loc='best')
    ax[1].scatter(
        subt,
        raw_residual,
        marker='.',
        label=f"Raw ({np.std(raw_residual, 0) * 100:.2f} %)",
        color=plt.cm.jet(0.25),
        alpha=0.25,
    )
    ax[1].scatter(
        subt,
        rel_residuals,
        marker='.',
        label=f"Detrended ({np.std(rel_residuals, 0) * 100:.2f} %)",
        color=plt.cm.jet(0.75),
        alpha=0.25,
    )
    ax[1].legend(loc='best')
    ax[1].set_xlabel('Time [BJD]')
    ax[0].set_title(f'Residual Statistics: {p} {selftype} {fltr}')
    ax[1].set_ylabel("Relative Flux")

    # compute fourier transform of raw_residual
    N = len(raw_residual)
    fft_raw = fft(raw_residual)
    fft_res = fft(rel_residuals)
    xf = fftfreq(len(raw_residual), d=np.diff(subt).mean() * 24 * 60 * 60)[
        : N // 2
    ]
    # fftraw = 2.0/N * np.abs(fft_raw[0:N//2])
    # future: square + integrate under the curve and normalize such that it equals time series variance
    ax[2].loglog(
        xf,
        2.0 / N * np.abs(fft_raw[0 : N // 2]),
        alpha=0.5,
        label='Raw',
        color=plt.cm.jet(0.25),
    )
    ax[2].loglog(
        xf,
        2.0 / N * np.abs(fft_res[0 : N // 2]),
        alpha=0.5,
        label='Detrended',
        color=plt.cm.jet(0.75),
    )
    if tdur_freq:
        ax[2].axvline(tdur_freq, ls='--', color='black', alpha=0.5)
    ax[2].set_ylabel('Power')
    ax[2].set_xlabel('Frequency [Hz]')
    ax[2].legend()
    ax[2].grid(True, ls='--')
    plt.tight_layout()

    # save plot to state vector
    savedplot = save_plot_tosv(fig)
    plt.close(fig)
    return savedplot


# --------------------------------------------------------------------
def rendertable(data, params, visitor: dawgie.Visitor) -> None:
    '''
    Helper function to render a table using data corresponding to
    the passed parameters
    '''
    clabels = ['name', 'estimate', 'units', 'description', 'ref']
    table = visitor.add_table(clabels=clabels, rows=len(params))
    # display stellar estimates
    for i, param in enumerate(params):
        table.get_cell(i, 0).add_primitive(param)
        table.get_cell(i, 1).add_primitive(data[param])
        param_proc = [
            ('_units', 'N/A'),
            ('_descr', 'No Description'),
            ('_ref', 'N/A'),
        ]
        for idx, (suffix, msg) in enumerate(param_proc):
            if param + suffix in data:
                table.get_cell(i, 2 + idx).add_primitive(data[param + suffix])
                pass
            else:
                table.get_cell(i, 2 + idx).add_primitive(msg)
            pass
        pass
    return


# --------------------------------------------------------------------
def barplot(title, categories, counts, categories2, counts2, visitor):
    '''barplot ds'''
    myfig = plt.figure()
    plt.title(title.replace('log(', '').replace(')', ''), fontsize=18)
    plt.bar(categories, counts, color='khaki', zorder=1, label='everything')
    plt.bar(
        categories2,
        counts2,
        color='olive',
        zorder=2,
        label='Roudier et al. 2021',
    )
    plt.ylabel('# of planets', fontsize=14)
    plt.xlabel('Spectral Type', fontsize=14)
    plt.legend()
    save_plot_toscreen(myfig, visitor)
    return


# --------------------------------------------------------------------
def distrplot(paramName, values1, values2, visitor, units=None):
    '''distrplot ds'''

    # clean up the values so that it's just an array of floats; no string error messages
    cleanvalues1 = []
    for value in values1:
        if len(str(value)) > 0:
            if str(value)[0].isdigit() or str(value)[0] == '-':
                cleanvalues1.append(value)
    cleanvalues1 = np.array(cleanvalues1, dtype=float)
    cleanvalues2 = []
    for value2 in values2:
        if len(str(value2)) > 0:
            if str(value2)[0].isdigit() or str(value2)[0] == '-':
                cleanvalues2.append(value2)
    cleanvalues2 = np.array(cleanvalues2, dtype=float)

    # most histograms are better on a log scale
    if paramName == 'luminosity':
        cleanvalues1 = np.log10(cleanvalues1)
        cleanvalues2 = np.log10(cleanvalues2)
        paramName = 'log(L*)'
    elif paramName in [
        'M*',
        'RHO*',
        'L*',
        'mass',
        'sma',
        'period',
        'density',
        'insolation',
        'H',
        'H_max',
        'modulation',
        'modulation_max',
        'ZFOM',
        'ZFOM_max',
        'v_wind',
        'rho_wind',
        'M_loss_rate_wind',
        'M_loss_rate_evap',
        'Beta_rad',
    ]:
        cleanvalues1 = np.log10(cleanvalues1)
        cleanvalues2 = np.log10(cleanvalues2)
        paramName = 'log(' + paramName + ')'

    # not sure why this is necessary.  why are some fields and entire params blank?
    # I guess it's the spTyp field, which seems to be missing from the resulting histograms
    if len(cleanvalues1) == 0:
        return
    if len(cleanvalues2) == 0:
        return

    myfig = plt.figure()
    plt.title(paramName.replace('log(', '').replace(')', ''), fontsize=18)
    outlier_aware_hist(
        cleanvalues1,
        cleanvalues2,
        *calculate_bounds(cleanvalues1),
        label1='everything',
        label2='Roudier et al. 2021',
    )
    plt.ylabel('# of planets', fontsize=14)
    # if units is None: plt.xlabel('Estimate')
    # else: plt.xlabel(f'Estimate [{units}]')
    if units is None:
        plt.xlabel(paramName, fontsize=14)
    else:
        plt.xlabel(paramName + f' [{units}]', fontsize=14)

    # debugging: save to disk!
    # plt.savefig('/proj/data/debug/histogramTest.png')

    save_plot_toscreen(myfig, visitor)
    return


# --------------------------------------------------------------------
def mad(data):
    '''mad ds'''
    median = np.nanmedian(data)
    diff = np.abs(data - median)
    mad_est = np.nanmedian(diff)
    return mad_est


# --------------------------------------------------------------------
def calculate_bounds(data, z_thresh=3.5):
    '''computes outlier cutoffs'''
    MAD = mad(data)
    median = np.nanmedian(data)
    const = z_thresh * MAD / 0.6745
    return (median - const, median + const)


# --------------------------------------------------------------------
def outlier_aware_hist(
    data,
    data2=None,
    lower=None,
    upper=None,
    color1='khaki',
    color2='olive',
    label1='',
    label2='',
    bins=15,
):
    '''note: code is originally from
    https://stackoverflow.com/questions/15837810/making-pyplot-hist-first-and-last-bins-include-outliers
    '''
    if not lower or lower < data.min():
        lower = data.min()
        lower_outliers = False
    else:
        lower_outliers = True

    if not upper or upper > data.max():
        upper = data.max()
        upper_outliers = False
    else:
        upper_outliers = True

    _, _, patches = plt.hist(
        data,
        range=(lower, upper),
        bins=bins,
        color=color1,
        zorder=1,
        label=label1,
    )

    if data2 is not None:
        plt.hist(
            data2,
            range=(lower, upper),
            bins=bins,
            color=color2,
            zorder=2,
            label=label2,
        )

    if lower_outliers:
        n_lower_outliers = (data < lower).sum()
        patches[0].set_height(patches[0].get_height() + n_lower_outliers)
        patches[0].set_facecolor('gold')
        patches[0].set_label(f'Lower outliers: ({data.min():.2f}, {lower:.2f})')

    if upper_outliers:
        n_upper_outliers = (data > upper).sum()
        patches[-1].set_height(patches[-1].get_height() + n_upper_outliers)
        patches[-1].set_facecolor('yellowgreen')
        patches[-1].set_label(
            f'Upper outliers: ({upper:.2f}, {data.max():.2f})'
        )

    if lower_outliers or upper_outliers:
        plt.legend()
    return


# --------------------------------------------------------------------
def plot_normalized_byvisit(data, vrange, visitor):

    for index, v in enumerate(data['visits']):
        wave = data['wave'][index]
        nspec = data['nspec'][index]

        myfig = plt.figure()
        plt.title('Visit: ' + str(v))
        for w, s in zip(wave, nspec):
            select = (w > np.min(vrange)) & (w < np.max(vrange))
            plt.plot(w[select], s[select], 'o')
        plt.ylabel('Normalized Flux')
        plt.xlabel('Wavelength [$\\mu$m]')
        plt.xlim(np.min(vrange), np.max(vrange))

        save_plot_toscreen(myfig, visitor)


# --------------------------------------------------------------------


def add_scale_height_labels(data, vspectrum, ax, fig):

    if 'Hs' in data:
        rpmed = np.sqrt(np.nanmedian(vspectrum))
        Hs = data['Hs'][0]

        # Retro compatibility for Hs in [m]
        if Hs > 1 and ('RSTAR' in data):
            Hs /= data['RSTAR'][0]

        axtwin = ax.twinx()
        axtwin.set_ylabel('$\\Delta$ [H$_s$]', fontsize=12)
        axmin, axmax = ax.get_ylim()

        if np.isnan(np.nanmax(vspectrum)):
            # log.warning(
            #    '--< PROBLEM: spectrum is all NaN %s %s >--',
            #    target, planet_letter)
            pass
        elif axmin >= 0:
            axtwin.set_ylim(
                (np.sqrt(1e-2 * axmin) - rpmed) / Hs,
                (np.sqrt(1e-2 * axmax) - rpmed) / Hs,
            )
        else:
            axtwin.set_ylim(
                (-np.sqrt(-1e-2 * axmin) - rpmed) / Hs,
                (np.sqrt(1e-2 * axmax) - rpmed) / Hs,
            )

        fig.tight_layout()


# --------------------------------------------------------------------


# --------------------------------------------------------------------


# --------------------------------------------------------------------
