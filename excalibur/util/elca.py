"""
Exoplanet light curve analysis

Fit an exoplanet transit model to time series data.

Heavily inspired by: https://github.com/rzellem/EXOTIC/blob/main/exotic/api/elca.py
"""

import copy
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator, NullLocator
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares, curve_fit, brentq
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from scipy import spatial
#from ultranest import ReactiveNestedSampler

# from astropy import units  # CI problem.  strange
import astropy.units
from astropy.time import Time

try:
    from astropy.modeling.models import BlackBody
except ImportError:
    from astropy.modeling.blackbody import blackbody_lambda as BlackBody


def weightedflux(flux, gw, nearest):
    """
    Calculate the weighted flux for each point in a light curve used for
    Gaussian process regression with nearest neighbors.

    Parameters:
        flux (array-like): flux values.
        gw (array-like): gaussian weights.
        nearest (array-like): nearest neighbors.

    Returns:
        array-like: weighted flux values.
    """
    return np.sum(flux[nearest] * gw, axis=-1)


def gaussian_weights(X, w=1, neighbors=50):
    """
    Calculate the gaussian weights for each point in a light curve used for
    nearest neighbor Gaussian process regression.

    Parameters:
        X (array-like): array of features like x/y centroid, noise pixel, etc
        w (float, optional): weight for each feature. Default is 1.
        neighbors (int, optional): number of nearest neighbors to use. Default is 50.

    Returns:
        tuple: a tuple containing:
            - array-like: gaussian weights.
            - array-like: nearest neighbors indices.
    """
    Xm = (X - np.median(X, 0)) * w
    kdtree = spatial.cKDTree(Xm)
    nearest = np.zeros((X.shape[0], neighbors))
    gw = np.zeros((X.shape[0], neighbors), dtype=float)
    for point in range(X.shape[0]):
        ind = kdtree.query(kdtree.data[point], neighbors + 1)[1][1:]
        dX = Xm[ind] - Xm[point]
        Xstd = np.std(dX, 0)
        gX = np.exp(-(dX**2) / (2 * Xstd**2))
        gwX = np.product(gX, 1)
        gw[point, :] = gwX / gwX.sum()
        nearest[point, :] = ind
    gw[np.isnan(gw)] = 0.01
    return gw, nearest.astype(int)


def get_phase(times, per, tmid):
    """
    Calculate the phase of a light curve with bounds of -0.25 to 0.75.
    Good for visualizing phase curves

    Parameters:
        times (array-like): time values.
        per (float): orbital period.
        tmid (float): mid-transit time.

    Returns:
        array-like: phase values.
    """
    return (times - tmid + 0.25 * per) / per % 1 - 0.25


def mc_a1(m_a2, sig_a2, transit_model, airmass, data, n=10000):
    """
    Calculate the mean and standard deviation of airmass-detrended data using a Monte Carlo approach.
    The airmass detrending function is given by:
        A1 * exp( A2 * airmass )

    Parameters:
        m_a2 (float): mean value of A2 parameter.
        sig_a2 (float): standard deviation of A2 parameter.
        transit_model (array-like): transit model.
        airmass (array-like): airmass values.
        data (array-like): observed data.
        n (int, optional): number of iterations in Monte Carlo. Default is 10000.

    Returns:
        tuple: a tuple containing:
            - mean of median-detrended data using Monte Carlo.
            - standard deviation of median-detrended data using Monte Carlo.
    """
    a2 = np.random.normal(m_a2, sig_a2, n)
    model = transit_model * np.exp(
        np.repeat(np.expand_dims(a2, 0), airmass.shape[0], 0).T * airmass
    )
    detrend = data / model
    return np.mean(np.median(detrend, 0)), np.std(np.median(detrend, 0))


def round_to_2(*args):
    """
    Round a number to 2 significant figures.

    Parameters:
        args (float): number to round.

    Returns:
        float: rounded number.
    """
    x = args[0]
    if len(args) == 1:
        y = args[0]
    else:
        y = args[1]
    if np.floor(y) >= 1.0:
        roundval = 2
    else:
        try:
            roundval = -int(np.floor(np.log10(abs(y)))) + 1
        except OverflowError:
            roundval = 1
    return round(x, roundval)


# average data into bins of dt from start to finish
def time_bin(time, flux, dt=1.0 / (60 * 24)):
    """
    Average data into bins of dt from start to finish.

    Parameters:
        time (array-like): time values.
        flux (array-like): flux values.
        dt (float, optional): bin size in days. Default is 1/60/24.

    Returns:
        tuple: a tuple containing:
            - array-like: binned time values.
            - array-like: binned flux values.
            - array-like: binned standard deviations.
    """
    bins = int(np.floor((max(time) - min(time)) / dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    bstds = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time) + i * dt)) & (
            time < (min(time) + (i + 1) * dt)
        )
        if mask.sum() > 0:
            bflux[i] = np.nanmean(flux[mask])
            btime[i] = np.nanmean(time[mask])
            bstds[i] = np.nanstd(flux[mask]) / (mask.sum() ** 0.5)
    zmask = (bflux == 0) | (btime == 0) | np.isnan(bflux) | np.isnan(btime)
    return btime[~zmask], bflux[~zmask], bstds[~zmask]


class lc_fitter:
    """Fit a transit data to a light curve model"""

    def __init__(
        self,
        time,
        data,
        dataerr,
        airmass,
        prior,
        bounds,
        neighbors=100,
        max_ncalls=1e6,
        mode='ns',
        verbose=True,
    ):
        self.time = time
        self.data = data
        self.dataerr = dataerr
        self.airmass = airmass
        self.prior = prior
        self.bounds = bounds
        self.max_ncalls = max_ncalls
        self.verbose = verbose
        self.mode = mode
        self.neighbors = neighbors
        if self.mode == "lm":
            self.fit_LM()
        elif self.mode == "ns":
            self.fit_nested()
        # variables to keep pylint happy, attribute-defined-outside-init
        self.transit = None
        self.time_upsample = None
        self.transit_upsample = None
        self.phase_upsample = None
        self.wf = None
        self.model = None
        self.detrended = None
        self.detrendederr = None
        self.residuals = None
        self.airmass_model = None
        self.phase = None
        self.chi2 = None
        self.bic = None
        self.res_stdev = None
        self.sdata = None
        self.quality = None
        self.duration_measured = None
        self.duration_expected = None

    def fit_LM(self):
        """Fit a transit model to a light curve using least squares"""
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])

        # trim data around predicted transit/eclipse time
        if np.ndim(self.airmass) == 2:
            print(
                f'Computing nearest neighbors and gaussian weights for {len(self.time)} npts...'
            )
            self.gw, self.nearest = gaussian_weights(
                self.airmass, neighbors=self.neighbors
            )

        def lc2min_nneighbor(pars):
            for i, par in enumerate(pars):
                self.prior[freekeys[i]] = par
            lightcurve = transit(self.time, self.prior)
            detrended = self.data / lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve * wf
            return ((self.data - model) / self.dataerr) ** 2

        def lc2min_airmass(pars):
            # chi-squared
            for i, par in enumerate(pars):
                self.prior[freekeys[i]] = par
            model = transit(self.time, self.prior)
            model *= np.exp(self.prior['a2'] * self.airmass)
            detrend = self.data / model  # used to estimate a1
            model *= np.median(detrend)
            return ((self.data - model) / self.dataerr) ** 2

        try:
            if np.ndim(self.airmass) == 2:
                res = least_squares(
                    lc2min_nneighbor,
                    x0=[self.prior[k] for k in freekeys],
                    bounds=[boundarray[:, 0], boundarray[:, 1]],
                    jac='3-point',
                    loss='linear',
                )
            else:
                res = least_squares(
                    lc2min_airmass,
                    x0=[self.prior[k] for k in freekeys],
                    bounds=[boundarray[:, 0], boundarray[:, 1]],
                    jac='3-point',
                    loss='linear',
                )
        except ValueError as e:
            print(
                f"{e} \nbounded light curve fitting failed...check priors "
                "(e.g. estimated mid-transit time + orbital period)"
            )

            for i, k in enumerate(freekeys):
                if not boundarray[i, 0] < self.prior[k] < boundarray[i, 1]:
                    print(
                        f"bound: [{boundarray[i, 0]}, {boundarray[i, 1]}] prior: {self.prior[k]}"
                    )

            print("removing bounds and trying again...")

            if np.ndim(self.airmass) == 2:
                res = least_squares(
                    lc2min_nneighbor,
                    x0=[self.prior[k] for k in freekeys],
                    method='lm',
                    jac='3-point',
                    loss='linear',
                )
            else:
                res = least_squares(
                    lc2min_airmass,
                    x0=[self.prior[k] for k in freekeys],
                    method='lm',
                    jac='3-point',
                    loss='linear',
                )

        self.parameters = copy.deepcopy(self.prior)
        self.errors = {}

        for i, k in enumerate(freekeys):
            self.parameters[k] = res.x[i]
            cov = res.jac.T @ res.jac
            cov = np.linalg.inv(cov)
            cov = cov * res.cost
            self.errors[k] = np.sqrt(cov[i, i])

        self.create_fit_variables()

    def create_fit_variables(self):
        """Create variables for the best fit model"""
        self.phase = get_phase(
            self.time, self.parameters['per'], self.parameters['tmid']
        )
        self.transit = transit(self.time, self.parameters)
        self.time_upsample = np.linspace(min(self.time), max(self.time), 1000)
        self.transit_upsample = transit(self.time_upsample, self.parameters)
        self.phase_upsample = get_phase(
            self.time_upsample, self.parameters['per'], self.parameters['tmid']
        )
        if self.mode == "ns" and np.ndim(self.airmass) == 1:
            self.parameters['a1'], self.errors['a1'] = mc_a1(
                self.parameters.get('a2', 0),
                self.errors.get('a2', 1e-6),
                self.transit,
                self.airmass,
                self.data,
            )
        if np.ndim(self.airmass) == 2:
            detrended = self.data / self.transit
            self.wf = weightedflux(detrended, self.gw, self.nearest)
            self.model = self.transit * self.wf
            self.detrended = self.data / self.wf
            self.detrendederr = self.dataerr / self.wf
        else:
            self.airmass_model = self.parameters['a1'] * np.exp(
                self.parameters.get('a2', 0) * self.airmass
            )
            self.model = self.transit * self.airmass_model
            self.detrended = self.data / self.airmass_model
            self.detrendederr = self.dataerr / self.airmass_model

        self.residuals = self.data - self.model
        self.res_stdev = np.std(self.residuals) / np.median(
            self.data
        )  # relative stdev
        self.chi2 = np.sum(self.residuals**2 / self.dataerr**2)
        self.bic = len(self.bounds) * np.log(len(self.time)) - 2 * np.log(
            self.chi2
        )

        # compare fit chi2 to smoothed data chi2
        dt = np.diff(np.sort(self.time)).mean()
        si = np.argsort(self.time)
        try:
            self.sdata = savgol_filter(
                self.data[si], 1 + 2 * int(0.5 / 24 / dt), 2
            )
        except ValueError:
            self.sdata = np.ones(len(self.time))

        schi2 = np.sum(
            (self.data[si] - self.sdata) ** 2 / self.dataerr[si] ** 2
        )
        self.quality = schi2 / self.chi2

        # measured duration
        tdur = (self.transit < 1).sum() * np.median(np.diff(np.sort(self.time)))

        # test for partial transit
        newtime = np.linspace(
            self.parameters['tmid'] - 0.2, self.parameters['tmid'] + 0.2, 10000
        )
        newtran = transit(newtime, self.parameters)
        masktran = newtran < 1
        newdur = np.diff(newtime).mean() * masktran.sum()

        self.duration_measured = tdur
        self.duration_expected = newdur

    def fit_nested(self):
        """Fit a transit model to a light curve using nested sampling"""
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray, 1).reshape(-1)

        # alloc data for best fit + error
        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        # trim data around predicted transit/eclipse time
        if np.ndim(self.airmass) == 2:
            print(
                f'Computing nearest neighbors and gaussian weights for {len(self.time)} npts...'
            )
            self.gw, self.nearest = gaussian_weights(
                self.airmass, neighbors=self.neighbors
            )

        def loglike_nneighbor(pars):
            """likelihood function for nested sampling using nearest neighbors regression for systematic detrending"""
            for i, par in enumerate(pars):
                self.prior[freekeys[i]] = par
            lightcurve = transit(self.time, self.prior)
            detrended = self.data / lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve * wf
            return -np.sum(((self.data - model) / self.dataerr) ** 2)

        def loglike_airmass(pars):
            """likelihood function for nested sampling using airmass detrending"""
            for i, par in enumerate(pars):
                self.prior[freekeys[i]] = par
            model = transit(self.time, self.prior)
            model *= np.exp(self.prior['a2'] * self.airmass)
            detrend = self.data / model  # used to estimate a1
            model *= np.median(detrend)
            return -np.sum(((self.data - model) / self.dataerr) ** 2)

        def loglike_phasecurve(pars):
            '''lc2min_phasecurve ds'''

            for i, _ in enumerate(pars):
                self.prior[freekeys[i]] = pars[i]
            lightcurve = phasecurve(self.time, self.prior)
            detrended = self.data / lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve * wf
            return -np.sum(((self.data - model) / self.dataerr) ** 2)

        def prior_transform(upars):
            # transform unit cube to prior volume
            return boundarray[:, 0] + bounddiff * upars

        def prior_transform_phasecurve(upars):
            '''prior_transform_phasecurve ds'''
            vals = boundarray[:, 0] + bounddiff * upars

            # set limits of phase amplitude to be less than eclipse depth or user bound
            edepth = (
                vals[freekeys.index('rprs')] ** 2 * vals[freekeys.index('fpfs')]
            )
            for k in ['c1', 'c2']:
                if k in freekeys:
                    # conditional prior needed to conserve energy
                    if k == 'c1':
                        ki = freekeys.index(k)
                        vals[ki] = upars[ki] * 0.4 * edepth + 0.1 * edepth
                    if k == 'c2':
                        ki = freekeys.index(k)
                        vals[ki] = upars[ki] * 0.25 * edepth - 0.125 * edepth
            return vals

        self.ns_type = 'ultranest'
        if np.ndim(self.airmass) == 2:
            if 'fpfs' in freekeys:
                test_sampler = ReactiveNestedSampler(
                    freekeys, loglike_phasecurve, prior_transform_phasecurve
                )
            else:
                test_sampler = ReactiveNestedSampler(
                    freekeys, loglike_nneighbor, prior_transform
                )
        else:
            test_sampler = ReactiveNestedSampler(
                freekeys, loglike_airmass, prior_transform
            )

        def noop(*args, **kwargs):

            pass

        if self.verbose is True:
            self.results = test_sampler.run(max_ncalls=int(self.max_ncalls))
        else:
            self.results = test_sampler.run(
                max_ncalls=int(self.max_ncalls),
                show_status=False,
                viz_callback=noop,
            )

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][
                i
            ]
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i],
            ]

        # final model
        self.create_fit_variables()

    def plot_bestfit(
        self, title="", bin_dt=30.0 / (60 * 24), zoom=False, phase=True
    ):
        """Plot the best fit model"""
        self.create_fit_variables()

        f = plt.figure(figsize=(9, 6))
        f.subplots_adjust(
            top=0.92, bottom=0.09, left=0.14, right=0.98, hspace=0
        )
        ax_lc = plt.subplot2grid((4, 5), (0, 0), colspan=5, rowspan=3)
        ax_res = plt.subplot2grid((4, 5), (3, 0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        axs[0].set_title(title)
        axs[0].set_ylabel("Relative Flux", fontsize=14)
        axs[0].grid(True, ls='--')

        rprs2 = self.parameters['rprs'] ** 2
        rprs2err = 2 * self.parameters['rprs'] * self.errors['rprs']

        lclabel1 = r"$R^{2}_{p}/R^{2}_{s}$ = %s $\pm$ %s" % (
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err)),
        )

        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" % (
            str(
                round_to_2(self.parameters['tmid'], self.errors.get('tmid', 0))
            ),
            str(round_to_2(self.errors.get('tmid', 0))),
        )

        lclabel = lclabel1 + "\n" + lclabel2

        if zoom:
            axs[0].set_ylim(
                [
                    1 - 1.25 * self.parameters['rprs'] ** 2,
                    1 + 0.5 * self.parameters['rprs'] ** 2,
                ]
            )
        else:
            if phase:
                axs[0].errorbar(
                    self.phase,
                    self.detrended,
                    yerr=np.std(self.residuals) / np.median(self.data),
                    ls='none',
                    marker='.',
                    color='black',
                    zorder=1,
                    alpha=0.2,
                )
            else:
                axs[0].errorbar(
                    self.time,
                    self.detrended,
                    yerr=np.std(self.residuals) / np.median(self.data),
                    ls='none',
                    marker='.',
                    color='black',
                    zorder=1,
                    alpha=0.2,
                )

        if phase:
            si = np.argsort(self.phase)
            bt2, br2, _ = time_bin(
                self.phase[si] * self.parameters['per'],
                self.residuals[si] / np.median(self.data) * 1e2,
                bin_dt,
            )
            axs[1].plot(
                self.phase,
                self.residuals / np.median(self.data) * 1e2,
                'k.',
                alpha=0.2,
                label=r'$\sigma$ = {:.2f} %'.format(
                    np.std(self.residuals / np.median(self.data) * 1e2)
                ),
            )
            axs[1].plot(
                bt2 / self.parameters['per'], br2, 'bs', alpha=1, zorder=2
            )
            axs[1].set_xlim([min(self.phase), max(self.phase)])
            axs[1].set_xlabel("Phase", fontsize=14)

            si = np.argsort(self.phase)
            bt2, bf2, bs = time_bin(
                self.phase[si] * self.parameters['per'],
                self.detrended[si],
                bin_dt,
            )
            axs[0].errorbar(
                bt2 / self.parameters['per'],
                bf2,
                yerr=bs,
                alpha=1,
                zorder=2,
                color='blue',
                ls='none',
                marker='s',
            )
            # axs[0].plot(self.phase[si], self.transit[si], 'r-', zorder=3, label=lclabel)
            sii = np.argsort(self.phase_upsample)
            axs[0].plot(
                self.phase_upsample[sii],
                self.transit_upsample[sii],
                'r-',
                zorder=3,
                label=lclabel,
            )
            axs[0].set_xlim([min(self.phase), max(self.phase)])
            axs[0].set_xlabel("Phase ", fontsize=14)
        else:
            bt, br, _ = time_bin(
                self.time, self.residuals / np.median(self.data) * 1e2, bin_dt
            )
            axs[1].plot(
                self.time,
                self.residuals / np.median(self.data) * 1e2,
                'k.',
                alpha=0.2,
                label=r'$\sigma$ = {:.2f} %'.format(
                    np.std(self.residuals / np.median(self.data) * 1e2)
                ),
            )
            axs[1].plot(
                bt,
                br,
                'bs',
                alpha=1,
                zorder=2,
                label=r'$\sigma$ = {:.2f} %'.format(np.std(br)),
            )
            axs[1].set_xlim([min(self.time), max(self.time)])
            axs[1].set_xlabel("Time [day]", fontsize=14)

            bt, bf, bs = time_bin(self.time, self.detrended, bin_dt)
            si = np.argsort(self.time)
            sii = np.argsort(self.time_upsample)
            axs[0].errorbar(
                bt,
                bf,
                yerr=bs,
                alpha=1,
                zorder=2,
                color='blue',
                ls='none',
                marker='s',
            )
            axs[0].plot(
                self.time_upsample[sii],
                self.transit_upsample[sii],
                'r-',
                zorder=3,
                label=lclabel,
            )
            axs[0].set_xlim([min(self.time), max(self.time)])
            axs[0].set_xlabel("Time [day]", fontsize=14)

        axs[0].get_xaxis().set_visible(False)
        axs[1].legend(loc='best')
        axs[0].legend(loc='best')
        axs[1].set_ylabel("Residuals [%]", fontsize=14)
        axs[1].grid(True, ls='--', axis='y')
        return f, axs

    def plot_posterior(self):
        """Plot the posterior distribution"""
        if self.ns_type == 'ultranest':
            ranges = []
            mask1 = np.ones(
                len(self.results['weighted_samples']['logl']), dtype=bool
            )
            mask2 = np.ones(
                len(self.results['weighted_samples']['logl']), dtype=bool
            )
            mask3 = np.ones(
                len(self.results['weighted_samples']['logl']), dtype=bool
            )
            titles = []
            labels = []
            flabels = {
                'rprs': r'R$_{p}$/R$_{s}$',
                'per': r'Period [day]',
                'tmid': r'T$_{mid}$',
                'ars': r'a/R$_{s}$',
                'inc': r'Inc. [deg]',
                'u1': r'u$_1$',
                'fpfs': r'F$_{p}$/F$_{s}$',
                'omega': r'$\omega$ [deg]',
                'mplanet': r'M$_{p}$ [M$_{\oplus}$]',
                'mstar': r'M$_{s}$ [M$_{\odot}$]',
                'ecc': r'$e$',
                'c0': r'$c_0$',
                'c1': r'$c_1$',
                'c2': r'$c_2$',
                'c3': r'$c_3$',
                'c4': r'$c_4$',
                'a0': r'$a_0$',
                'a1': r'$a_1$',
                'a2': r'$a_2$',
            }
            for i, key in enumerate(self.quantiles):
                labels.append(flabels.get(key, key))
                titles.append(
                    f"{self.parameters[key]:.5f} +- {self.errors[key]:.5f}"
                )
                ranges.append(
                    [
                        self.parameters[key] - 5 * self.errors[key],
                        self.parameters[key] + 5 * self.errors[key],
                    ]
                )

                if key in ('a1', 'a2'):
                    continue

                mask3 = (
                    mask3
                    & (
                        self.results['weighted_samples']['points'][:, i]
                        > (self.parameters[key] - 3 * self.errors[key])
                    )
                    & (
                        self.results['weighted_samples']['points'][:, i]
                        < (self.parameters[key] + 3 * self.errors[key])
                    )
                )

                mask1 = (
                    mask1
                    & (
                        self.results['weighted_samples']['points'][:, i]
                        > (self.parameters[key] - self.errors[key])
                    )
                    & (
                        self.results['weighted_samples']['points'][:, i]
                        < (self.parameters[key] + self.errors[key])
                    )
                )

                mask2 = (
                    mask2
                    & (
                        self.results['weighted_samples']['points'][:, i]
                        > (self.parameters[key] - 2 * self.errors[key])
                    )
                    & (
                        self.results['weighted_samples']['points'][:, i]
                        < (self.parameters[key] + 2 * self.errors[key])
                    )
                )

            chi2 = self.results['weighted_samples']['logl'] * -2
            fig = corner(
                self.results['weighted_samples']['points'],
                labels=labels,
                bins=int(np.sqrt(self.results['samples'].shape[0])),
                plot_range=ranges,
                # quantiles=(0.1, 0.84),
                plot_contours=True,
                levels=[
                    np.percentile(chi2[mask1], 95),
                    np.percentile(chi2[mask2], 95),
                    np.percentile(chi2[mask3], 95),
                ],
                titles=titles,
                data_kwargs={
                    'c': chi2,
                    'vmin': np.percentile(chi2[mask3], 1),
                    'vmax': np.percentile(chi2[mask3], 95),
                    'cmap': 'viridis',
                },
                label_kwargs={
                    'labelpad': 15,
                },
                hist_kwargs={
                    'color': 'black',
                },
            )
        else:
            raise ValueError(f"ns_type {self.ns_type} not recognized")
        return fig, None

    def plot_btempcurve(self, bandpass='IRAC 3.6um'):
        '''plot_btempcurve ds'''
        fig = plt.figure(figsize=(13, 7))
        ax_lc = plt.subplot2grid((4, 5), (0, 0), colspan=5, rowspan=3)
        ax_res = plt.subplot2grid((4, 5), (3, 0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        phase = (self.time - self.parameters['tmid']) / self.parameters['per']
        bin_dt = 10.0 / 24.0 / 60.0
        bt, bf, _ = time_bin(self.time, self.detrended, bin_dt)
        bp = (bt - self.parameters['tmid']) / self.parameters['per']
        bt, br, _ = time_bin(self.time, self.residuals, bin_dt)

        bcurve = brightness(bt, self.parameters)
        ogfpfs = self.parameters['fpfs']
        tbcurve = np.ones(bcurve.shape)
        for i, bc in enumerate(bcurve):
            self.parameters['fpfs'] = max(
                (bc - 1) / self.parameters['rprs'] ** 2, 0.00001
            )
            try:
                tbcurve[i] = brightnessTemp(self.parameters, bandpass)
            except ValueError:
                tbcurve[i] = np.nan
        self.parameters['fpfs'] = ogfpfs

        # residuals
        axs[1].plot(
            phase,
            self.residuals / np.median(self.data) * 1e6,
            'k.',
            alpha=0.15,
            label=fr'$\sigma$ = {np.std(self.residuals / np.median(self.data) * 1e6):.0f} ppm',
        )

        axs[1].plot(
            bp,
            1e6 * br / np.median(self.data),
            'w.',
            zorder=2,
            label=fr'$\sigma$ = {np.std(1e6 * br / np.median(self.data)):.0f} ppm',
        )

        axs[1].set_xlim([min(phase), max(phase)])
        axs[1].set_xlabel("Phase")
        axs[1].legend(loc='best')
        axs[1].set_ylabel("Residuals [ppm]")
        axs[1].grid(True, ls='--')

        axs[0].errorbar(
            phase,
            self.detrended,
            yerr=np.std(self.residuals) / np.median(self.data),
            ls='none',
            marker='.',
            color='black',
            alpha=0.1,
            zorder=1,
        )

        # map color to equilibrium temperature
        im = axs[0].scatter(
            bp,
            bf,
            marker='o',
            c=tbcurve,
            vmin=500,
            vmax=2750,
            cmap='jet',
            zorder=2,
            s=20,
        )
        cbar = plt.colorbar(im)
        cbar.ax.set_xlabel("B. Temp. [K]")

        axs[0].plot(phase, self.transit, 'w--', zorder=3)
        axs[0].set_xlim([min(phase), max(phase)])
        axs[0].set_xlabel("Phase ")

        axs[0].set_ylabel("Relative Flux")
        axs[0].grid(True, ls='--')
        axs[0].set_ylim([0.955, 1.03])

        plt.tight_layout()
        return fig, axs

    def plot_pixelmap(self, title='', savedir=None):
        '''plot_pixelmap ds'''
        fig, ax = plt.subplots(1, figsize=(8.5, 7))
        # airmass = system parameters
        xcent = self.airmass[:, 0]  # weighted flux x-cntroid
        ycent = self.airmass[:, 1]  # weighted flux y-centroid
        npp = self.airmass[:, 2]  # noise pixel parameter
        normpp = (npp - npp.min()) / (npp.max() - npp.min())  # norm btwn 0-1
        normpp *= 20
        normpp += 20
        im = ax.scatter(
            xcent,
            ycent,
            c=self.wf / np.median(self.wf),
            marker='.',
            vmin=0.99,
            vmax=1.01,
            alpha=0.5,
            cmap='jet',
            s=normpp,
        )
        ax.set_xlim(
            [
                np.median(xcent) - 3 * np.std(xcent),
                np.median(xcent) + 3 * np.std(xcent),
            ]
        )
        ax.set_ylim(
            [
                np.median(ycent) - 3 * np.std(ycent),
                np.median(ycent) + 3 * np.std(ycent),
            ]
        )

        ax.set_title(title, fontsize=14)
        ax.set_xlabel('X-Centroid [px]', fontsize=14)
        ax.set_ylabel('Y-Centroid [px]', fontsize=14)
        cbar = fig.colorbar(im)
        cbar.set_label(
            'Relative Pixel Response', fontsize=14, rotation=270, labelpad=15
        )

        plt.tight_layout()
        if savedir:
            plt.savefig(savedir + title + ".png")
            plt.close()
        return fig, ax


def brightnessTemp(priors, f='IRAC 3.6um'):
    '''Solve for Tb using Fp/Fs, Ts and a filter bandpass'''
    if '3.6' in f or '36' in f:
        waveset = np.linspace(3.15, 3.9, 1000) * astropy.units.micron
    else:
        waveset = np.linspace(4, 5, 1000) * astropy.units.micron

    def f2min(T, *args):
        fpfs, tstar, waveset = args
        fstar = BlackBody(waveset, tstar * astropy.units.K)
        fplanet = BlackBody(waveset, T * astropy.units.K)
        fp = np.trapz(fplanet, waveset)
        fs = np.trapz(fstar, waveset)
        return (fp / fs) - fpfs

    tb = brentq(f2min, 1, 4500, args=(priors['fpfs'], priors['T*'], waveset))
    return tb


class glc_fitter(lc_fitter):
    """Fit multiple light curves simultaneously with global and local parameters"""

    ns_type = 'ultranest'

    def __init__(
        self,
        input_data,
        global_bounds,
        local_bounds,
        max_ncalls=1e6,
        individual_fit=False,
        stdev_cutoff=0.03,
        verbose=False,
    ):
        # keys for input_data: time, flux, ferr, airmass, priors all numpy arrays
        self.lc_data = copy.deepcopy(input_data)
        self.global_bounds = global_bounds
        self.local_bounds = local_bounds
        self.individual_fit = individual_fit
        self.max_ncalls = max_ncalls

        self.stdev_cutoff = stdev_cutoff
        self.verbose = verbose

        self.fit_nested()

    def fit_nested(self):

        # create bound arrays for generating samples
        nobs = len(self.lc_data)
        gfreekeys = list(self.global_bounds.keys())

        # if isinstance(self.local_bounds, dict):
        #     lfreekeys = list(self.local_bounds.keys())
        #     boundarray = np.vstack([ [self.global_bounds[k] for k in gfreekeys], [self.local_bounds[k] for k in lfreekeys]*nobs ])
        # else:
        #     # if list type
        lfreekeys = []
        boundarray = [self.global_bounds[k] for k in gfreekeys]
        for i in range(nobs):
            lfreekeys.append(list(self.local_bounds[i].keys()))
            boundarray.extend([self.local_bounds[i][k] for k in lfreekeys[-1]])
        boundarray = np.array(boundarray)

        # fit individual light curves to constrain priors
        if self.individual_fit:
            for i in range(nobs):

                print(f"Fitting individual light curve {i + 1} / {nobs}")
                try:
                    mybounds = {**self.local_bounds[i], **self.global_bounds}
                except TypeError:
                    mybounds = {}
                    for k in self.local_bounds[i]:
                        mybounds[k] = self.local_bounds[i][k]
                    for k in self.global_bounds:
                        mybounds[k] = self.global_bounds[k]
                if 'per' in mybounds:
                    del mybounds['per']
                if 'inc' in mybounds and 'rprs' in mybounds:
                    del mybounds['inc']
                if 'tmid' in mybounds:
                    # find the closet mid transit time to the last observation
                    phase = (
                        self.lc_data[i]['time'][-1]
                        - self.lc_data[i]['priors']['tmid']
                    ) / self.lc_data[i]['priors']['per']
                    nepochs = np.round(phase)
                    newtmid = (
                        self.lc_data[i]['priors']['tmid']
                        + nepochs * self.lc_data[i]['priors']['per']
                    )
                    err = np.diff(mybounds['tmid'])[0] / 2.0
                    mybounds['tmid'] = [newtmid - err, newtmid + err]

                # fit individual light curve
                myfit = lc_fitter(
                    self.lc_data[i]['time'],
                    self.lc_data[i]['flux'],
                    self.lc_data[i]['ferr'],
                    self.lc_data[i]['airmass'],
                    self.lc_data[i]['priors'],
                    mybounds,
                )

                # check stdev_cutoff and residuals
                if myfit.res_stdev > self.stdev_cutoff:
                    print(
                        f"WARNING: Stdev of residuals is large! {myfit.res_stdev:.3f} > {self.stdev_cutoff:.3f}"
                    )

                # copy data over for individual fits
                self.lc_data[i]['individual'] = myfit.parameters.copy()
                self.lc_data[i]['individual_err'] = myfit.errors.copy()
                self.lc_data[i]['res_stdev'] = myfit.res_stdev
                self.lc_data[i]['quality'] = myfit.quality

                ti = sum(len(self.local_bounds[k]) for k in range(i))

                # update local priors
                for j, key in enumerate(self.local_bounds[i].keys()):

                    boundarray[j + ti + len(gfreekeys), 0] = (
                        myfit.parameters[key] - 5 * myfit.errors[key]
                    )
                    boundarray[j + ti + len(gfreekeys), 1] = (
                        myfit.parameters[key] + 5 * myfit.errors[key]
                    )
                    self.local_bounds[i][key] = [
                        myfit.parameters[key] - 5 * myfit.errors[key],
                        myfit.parameters[key] + 5 * myfit.errors[key],
                    ]

                    if key == 'rprs':
                        boundarray[j + ti + len(gfreekeys), 0] = max(
                            0, myfit.parameters[key] - 5 * myfit.errors[key]
                        )
                        self.local_bounds[i][key][0] = max(
                            0, myfit.parameters[key] - 5 * myfit.errors[key]
                        )

                # print name and stdev of residuals
                mint = np.min(self.lc_data[i]['time'])
                maxt = np.max(self.lc_data[i]['time'])
                try:
                    print(
                        f"{self.lc_data[i]['name']} & {Time(mint, format='jd').isot} & {Time(maxt, format='jd').isot} & {np.std(myfit.residuals)} & {len(self.lc_data[i]['time'])}"
                    )
                except ValueError:
                    print(
                        f"{self.lc_data[i]['name']} & {mint} & {maxt} & {np.std(myfit.residuals)} & {len(self.lc_data[i]['time'])}"
                    )

                del myfit

        # transform unit cube to prior volume
        bounddiff = np.diff(boundarray, 1).reshape(-1)

        def prior_transform(upars):
            return boundarray[:, 0] + bounddiff * upars

        def loglike(pars):
            chi2 = 0

            # for each light curve
            for i in range(nobs):

                # global keys
                for j, key in enumerate(gfreekeys):
                    self.lc_data[i]['priors'][key] = pars[j]

                # local keys
                ti = sum(len(self.local_bounds[k]) for k in range(i))
                for j, key in enumerate(lfreekeys[i]):
                    self.lc_data[i]['priors'][key] = pars[
                        j + ti + len(gfreekeys)
                    ]

                # compute model
                model = transit(
                    self.lc_data[i]['time'], self.lc_data[i]['priors']
                )
                model *= np.exp(
                    self.lc_data[i]['priors']['a2'] * self.lc_data[i]['airmass']
                )
                detrend = self.lc_data[i]['flux'] / model
                model *= np.mean(detrend)

                chi2 += np.sum(
                    (
                        (self.lc_data[i]['flux'] - model)
                        / self.lc_data[i]['ferr']
                    )
                    ** 2
                )

            # maximization metric for nested sampling
            return -0.5 * chi2

        freekeys = [] + gfreekeys
        for n in range(nobs):
            for k in lfreekeys[n]:
                freekeys.append(f"local_{k}_{n}")

        if self.verbose:
            self.results = ReactiveNestedSampler(
                freekeys, loglike, prior_transform
            ).run(max_ncalls=self.max_ncalls)
        else:
            self.results = ReactiveNestedSampler(
                freekeys, loglike, prior_transform
            ).run(
                max_ncalls=self.max_ncalls,
                show_status=self.verbose,
                viz_callback=self.verbose,
            )

        self.quantiles = {}
        self.errors = {}
        self.parameters = self.lc_data[0]['priors'].copy()

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][
                i
            ]
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i],
            ]

        for n in range(nobs):
            self.lc_data[n]['errors'] = {}
            for k in lfreekeys[n]:
                pkey = f"local_{k}_{n}"

                self.lc_data[n]['priors'][k] = self.parameters[pkey]
                self.lc_data[n]['errors'][k] = self.errors[pkey]

                if k == 'rprs' and 'rprs' not in freekeys:
                    self.parameters[k] = self.lc_data[n]['priors'][k]
                    self.errors[k] = self.lc_data[n]['errors'][k]

            # solve for a1
            model = transit(self.lc_data[n]['time'], self.lc_data[n]['priors'])
            airmass = np.exp(
                self.lc_data[n]['airmass'] * self.lc_data[n]['priors']['a2']
            )
            detrend = self.lc_data[n]['flux'] / (model * airmass)
            self.lc_data[n]['priors']['a1'] = np.mean(detrend)
            self.lc_data[n]['residuals'] = (
                self.lc_data[n]['flux']
                - model * airmass * self.lc_data[n]['priors']['a1']
            )
            self.lc_data[n]['detrend'] = self.lc_data[n]['flux'] / (
                airmass * self.lc_data[n]['priors']['a1']
            )
            # phase
            self.lc_data[n]['phase'] = get_phase(
                self.lc_data[n]['time'],
                self.lc_data[n]['priors']['per'],
                self.lc_data[n]['priors']['tmid'],
            )

    def plot_bestfits(self):
        """plot best fit"""
        nrows = len(self.lc_data) // 4 + 1
        fig, ax = plt.subplots(nrows, 4, figsize=(5 + 5 * nrows, 5 * nrows))

        # turn off all axes
        for i in range(nrows * 4):
            ri = int(i / 4)
            ci = i % 4
            if ax.ndim == 1:
                ax[i].axis('off')
            else:
                ax[ri, ci].axis('off')

        # cycle the colors and markers
        markers = cycle(
            ['o', 'v', '^', '<', '>', 's', '*', 'h', 'H', 'D', 'd', 'P', 'X']
        )
        colors = cycle(
            [
                'black',
                'blue',
                'green',
                'orange',
                'purple',
                'grey',
                'magenta',
                'cyan',
                'lime',
            ]
        )

        # plot observations
        for i, lc_data in enumerate(self.lc_data):
            ri = int(i / 4)
            ci = i % 4
            ncolor = next(colors)
            nmarker = next(markers)

            model = transit(lc_data['time'], lc_data['priors'])
            airmass = np.exp(
                self.lc_data[i]['airmass'] * self.lc_data[i]['priors']['a2']
            )
            detrend = self.lc_data[i]['flux'] / (model * airmass)

            if ax.ndim == 1:
                ax[i].axis('on')
                ax[i].errorbar(
                    self.lc_data[i]['time'],
                    self.lc_data[i]['flux'] / airmass / detrend.mean(),
                    yerr=self.lc_data[i]['ferr'] / airmass / detrend.mean(),
                    ls='none',
                    marker=nmarker,
                    color=ncolor,
                    alpha=0.5,
                )
                ax[i].plot(self.lc_data[i]['time'], model, 'r-', zorder=2)
                ax[i].set_xlabel("Time [BJD]", fontsize=14)
                ax[i].set_ylabel("Relative Flux", fontsize=14)
                ax[i].set_title(
                    f"{self.lc_data[i].get('name', '')}", fontsize=16
                )
            else:
                ax[ri, ci].axis('on')
                ax[ri, ci].errorbar(
                    self.lc_data[i]['time'],
                    self.lc_data[i]['flux'] / airmass / detrend.mean(),
                    yerr=self.lc_data[i]['ferr'] / airmass / detrend.mean(),
                    ls='none',
                    marker=nmarker,
                    color=ncolor,
                    alpha=0.5,
                )
                ax[ri, ci].plot(self.lc_data[i]['time'], model, 'r-', zorder=2)
                ax[ri, ci].set_xlabel("Time[BJD]", fontsize=14)
                ax[ri, ci].set_ylabel("Relative Flux", fontsize=14)
                ax[ri, ci].set_title(
                    f"{self.lc_data[i].get('name', '')}", fontsize=16
                )

        plt.tight_layout()
        return fig

    def plot_bestfit(
        self,
        title="",
        bin_dt=30.0 / (60 * 24),
        alpha=0.05,
        ylim_sigma=6,
        legend_loc='best',
        phase_limits='median',
    ):
        f = plt.figure(figsize=(15, 12))
        f.subplots_adjust(top=0.92, bottom=0.09, left=0.1, right=0.98, hspace=0)
        ax_lc = plt.subplot2grid((4, 5), (0, 0), colspan=5, rowspan=3)
        ax_res = plt.subplot2grid((4, 5), (3, 0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        axs[0].set_title(title)
        axs[0].set_ylabel("Relative Flux", fontsize=14)
        axs[0].grid(True, ls='--')

        try:
            rprs2 = self.parameters['rprs'] ** 2
            rprs2err = 2 * self.parameters['rprs'] * self.errors['rprs']
        except KeyError:
            rprs2 = self.lc_data[0]['priors']['rprs'] ** 2
            rprs2err = (
                2
                * self.lc_data[0]['priors']['rprs']
                * self.lc_data[0]['errors']['rprs']
            )

        lclabel1 = r"$R^{2}_{p}/R^{2}_{s}$ = %s $\pm$ %s" % (
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err)),
        )

        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" % (
            str(
                round_to_2(self.parameters['tmid'], self.errors.get('tmid', 0))
            ),
            str(round_to_2(self.errors.get('tmid', 0))),
        )

        lclabel = lclabel1 + "\n" + lclabel2
        minp = 1
        maxp = 0

        min_std = 1
        # cycle the colors and markers
        markers = cycle(
            ['o', 'v', '^', '<', '>', 's', '*', 'h', 'H', 'D', 'd', 'P', 'X']
        )
        colors = cycle(
            [
                'black',
                'blue',
                'green',
                'orange',
                'purple',
                'grey',
                'magenta',
                'cyan',
                'lime',
            ]
        )

        alldata = {
            'time': [],
            'flux': [],
            'detrend': [],
            'ferr': [],
            'residuals': [],
        }

        for n, lc_data in enumerate(self.lc_data):
            ncolor = next(colors)
            nmarker = next(markers)
            alldata['time'].extend(lc_data['time'].tolist())
            alldata['detrend'].extend(self.lc_data[n]['detrend'].tolist())
            alldata['flux'].extend(self.lc_data[n]['flux'].tolist())
            alldata['ferr'].extend(self.lc_data[n]['ferr'].tolist())
            alldata['residuals'].extend(self.lc_data[n]['residuals'].tolist())

            phase = get_phase(
                self.lc_data[n]['time'],
                self.parameters['per'],
                self.lc_data[n]['priors']['tmid'],
            )
            si = np.argsort(phase)

            # plot data
            axs[0].errorbar(
                phase,
                self.lc_data[n]['detrend'],
                yerr=np.std(self.lc_data[n]['residuals'])
                / np.median(self.lc_data[n]['flux']),
                ls='none',
                marker=nmarker,
                color=ncolor,
                zorder=1,
                alpha=alpha,
            )

            # plot residuals
            axs[1].plot(
                phase,
                self.lc_data[n]['residuals']
                / np.median(self.lc_data[n]['flux'])
                * 1e2,
                color=ncolor,
                marker=nmarker,
                ls='none',
                alpha=0.2,
            )

            # plot binned data
            bt2, bf2, bs = time_bin(
                phase[si] * self.lc_data[n]['priors']['per'],
                self.lc_data[n]['detrend'][si],
                bin_dt,
            )
            axs[0].errorbar(
                bt2 / self.lc_data[n]['priors']['per'],
                bf2,
                yerr=bs,
                alpha=1,
                zorder=2,
                color=ncolor,
                ls='none',
                marker=nmarker,
                label=r'{}: {:.2f} %'.format(
                    self.lc_data[n].get('name', ''),
                    np.std(
                        self.lc_data[n]['residuals']
                        / np.median(self.lc_data[n]['flux'])
                        * 1e2
                    ),
                ),
            )

            # replace min and max for upsampled lc model
            minp = min(minp, phase)
            maxp = max(maxp, phase)
            min_std = min(
                min_std,
                np.std(
                    self.lc_data[n]['residuals']
                    / np.median(self.lc_data[n]['flux'])
                ),
            )

        # create binned plot for all the data

        for k in alldata:
            alldata[k] = np.array(alldata[k])

        phase = get_phase(
            alldata['time'],
            self.parameters['per'],
            self.lc_data[-1]['priors']['tmid'],
        )
        si = np.argsort(phase)
        bt, br, _ = time_bin(
            phase[si] * self.parameters['per'],
            alldata['residuals'][si] / np.median(alldata['flux']),
            2 * bin_dt,
        )
        bt, bf, bs = time_bin(
            phase[si] * self.parameters['per'],
            alldata['detrend'][si],
            2 * bin_dt,
        )

        axs[0].errorbar(
            bt / self.parameters['per'],
            bf,
            yerr=bs,
            alpha=1,
            zorder=2,
            color='white',
            ls='none',
            marker='o',
            ms=15,
            markeredgecolor='black',
            ecolor='black',
            label=r'Binned Data: {:.2f} %'.format(np.std(br) * 1e2),
        )

        axs[1].plot(
            bt / self.parameters['per'],
            br * 1e2,
            color='white',
            ls='none',
            marker='o',
            ms=11,
            markeredgecolor='black',
        )

        # best fit model
        self.time_upsample = np.linspace(
            minp * self.parameters['per'] + self.parameters['tmid'],
            maxp * self.parameters['per'] + self.parameters['tmid'],
            10000,
        )
        self.transit_upsample = transit(
            self.time_upsample, self.lc_data[0]['priors']
        )
        self.phase_upsample = get_phase(
            self.time_upsample, self.parameters['per'], self.parameters['tmid']
        )
        sii = np.argsort(self.phase_upsample)
        axs[0].plot(
            self.phase_upsample[sii],
            self.transit_upsample[sii],
            'r-',
            zorder=3,
            label=lclabel,
            lw=3,
        )

        axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        axs[0].set_xlabel("Phase ", fontsize=14)
        axs[0].set_ylim(
            [
                1 - self.parameters['rprs'] ** 2 - ylim_sigma * min_std,
                1 + ylim_sigma * min_std,
            ]
        )
        axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        axs[1].set_xlabel("Phase", fontsize=14)
        axs[1].set_ylim([-6 * min_std * 1e2, 6 * min_std * 1e2])

        # compute average min and max for all the data
        mins = []
        maxs = []
        for n, lc_data in enumerate(self.lc_data):
            mins.append(min(self.lc_data[n]['phase']))
            maxs.append(max(self.lc_data[n]['phase']))

        # set up phase limits
        if isinstance(phase_limits, str):
            if phase_limits == "minmax":
                axs[0].set_xlim(
                    [min(self.phase_upsample), max(self.phase_upsample)]
                )
                axs[1].set_xlim(
                    [min(self.phase_upsample), max(self.phase_upsample)]
                )
            elif phase_limits == "median":
                axs[0].set_xlim([np.median(mins), np.median(maxs)])
                axs[1].set_xlim([np.median(mins), np.median(maxs)])
            else:
                axs[0].set_xlim(
                    [min(self.phase_upsample), max(self.phase_upsample)]
                )
                axs[1].set_xlim(
                    [min(self.phase_upsample), max(self.phase_upsample)]
                )
        elif isinstance(phase_limits, list):
            axs[0].set_xlim([phase_limits[0], phase_limits[1]])
            axs[1].set_xlim([phase_limits[0], phase_limits[1]])
        elif isinstance(phase_limits, tuple):
            axs[0].set_xlim([phase_limits[0], phase_limits[1]])
            axs[1].set_xlim([phase_limits[0], phase_limits[1]])
        else:
            axs[0].set_xlim(
                [min(self.phase_upsample), max(self.phase_upsample)]
            )
            axs[1].set_xlim(
                [min(self.phase_upsample), max(self.phase_upsample)]
            )

        axs[0].get_xaxis().set_visible(False)
        axs[0].legend(loc=legend_loc, ncol=len(self.lc_data) // 7 + 1)
        axs[1].set_ylabel("Residuals [%]", fontsize=14)
        axs[1].grid(True, ls='--', axis='y')
        return f, axs

    def plot_stack(self, title="", bin_dt=30.0 / (60 * 24), dy=0.02):
        """plot mosaic of light curve fits"""
        f, ax = plt.subplots(1, figsize=(9, 12))

        ax.set_title(title)
        ax.set_ylabel("Relative Flux", fontsize=14)
        ax.grid(True, ls='--')

        minp = 1
        maxp = 0

        min_std = 1
        # cycle the colors and markers
        markers = cycle(
            ['o', 'v', '^', '<', '>', 's', '*', 'h', 'H', 'D', 'd', 'P', 'X']
        )
        colors = cycle(
            [
                'black',
                'blue',
                'green',
                'orange',
                'purple',
                'grey',
                'magenta',
                'cyan',
                'lime',
            ]
        )
        for n, lc_data in enumerate(self.lc_data):
            ncolor = next(colors)
            nmarker = next(markers)

            phase = get_phase(
                lc_data['time'],
                self.parameters['per'],
                lc_data['priors']['tmid'],
            )
            si = np.argsort(phase)
            bt2, _, _ = time_bin(
                phase[si] * self.parameters['per'],
                self.lc_data[n]['residuals'][si]
                / np.median(self.lc_data[n]['flux'])
                * 1e2,
                bin_dt,
            )

            # plot data
            ax.errorbar(
                phase,
                self.lc_data[n]['detrend'] - n * dy,
                yerr=np.std(self.lc_data[n]['residuals'])
                / np.median(self.lc_data[n]['flux']),
                ls='none',
                marker=nmarker,
                color=ncolor,
                zorder=1,
                alpha=0.25,
            )

            # plot binned data
            bt2, bf2, bs = time_bin(
                phase[si] * self.lc_data[n]['priors']['per'],
                self.lc_data[n]['detrend'][si] - n * dy,
                bin_dt,
            )
            ax.errorbar(
                bt2 / self.lc_data[n]['priors']['per'],
                bf2,
                yerr=bs,
                alpha=1,
                zorder=2,
                color=ncolor,
                ls='none',
                marker=nmarker,
            )

            # replace min and max for upsampled lc model
            minp = min(minp, phase)
            maxp = max(maxp, phase)
            min_std = min(
                min_std,
                np.std(
                    self.lc_data[n]['residuals']
                    / np.median(self.lc_data[n]['flux'])
                ),
            )

            # best fit model
            self.time_upsample = np.linspace(
                minp * self.parameters['per'] + self.parameters['tmid'],
                maxp * self.parameters['per'] + self.parameters['tmid'],
                10000,
            )
            self.transit_upsample = transit(self.time_upsample, self.parameters)
            self.phase_upsample = get_phase(
                self.time_upsample,
                self.parameters['per'],
                self.parameters['tmid'],
            )
            sii = np.argsort(self.phase_upsample)
            ax.plot(
                self.phase_upsample[sii],
                self.transit_upsample[sii] - n * dy,
                ls='-',
                color=ncolor,
                zorder=3,
                label=self.lc_data[n].get('name', ''),
            )

        n = len(self.lc_data) - 1
        ax.set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        ax.set_xlabel("Phase ", fontsize=14)
        ax.set_ylim(
            [
                1 - self.parameters['rprs'] ** 2 - 5 * min_std - n * dy,
                1 + 5 * min_std,
            ]
        )
        ax.get_xaxis().set_visible(False)
        ax.legend(loc='best')
        return f, ax


# rip off of corner.py so we can cmap scatter to chi2
def corner(
    xs,
    bins=20,
    plot_range=None,
    weights=None,
    color="k",
    hist_bin_factor=1,
    smooth1d=None,
    levels=None,
    labels=None,
    label_kwargs=None,
    titles=None,
    title_kwargs=None,
    truths=None,
    truth_color="#4682b4",
    scale_hist=False,
    quantiles=None,
    fig=None,
    max_n_ticks=5,
    top_ticks=False,
    use_math_text=False,
    reverse=False,
    hist_kwargs=None,
    **hist2d_kwargs,
):
    """Make a *sick* corner plot"""
    if quantiles is None:
        quantiles = []
    if title_kwargs is None:
        title_kwargs = {}
    if label_kwargs is None:
        label_kwargs = {}

    # Try filling in labels from pandas.DataFrame columns.
    if labels is None:
        try:
            labels = xs.columns
        except AttributeError:
            pass

    # Deal with 1D sample lists.
    xs = np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = np.atleast_2d(xs)
    else:
        xs = xs.T

    # Parse the weight array.
    if weights is not None:
        weights = np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("Weights must be 1-D")
        if xs.shape[1] != weights.shape[0]:
            raise ValueError("Lengths of weights must match number of samples")

    # Parse the parameter ranges.
    if plot_range is None:
        if "extents" in hist2d_kwargs:
            print("Deprecated keyword argument 'extents'. Use 'range' instead.")
            plot_range = hist2d_kwargs.pop("extents")
        else:
            plot_range = [[x.min(), x.max()] for x in xs]
            # Check for parameters that never change.
            m = np.array([e[0] == e[1] for e in plot_range], dtype=bool)
            if np.any(m):

                raise ValueError(
                    (
                        "It looks like the parameter(s) in "
                        "column(s) {0} have no dynamic range. "
                        "Please provide a `plot_range` argument."
                    ).format(", ".join(map("{0}".format, np.arange(len(m))[m])))
                )

    else:
        # If any of the extents are percentiles, convert them to ranges.
        # Also make sure it's a normal list.
        plot_range = list(plot_range)
        # GMR: No comment

        for i, _ in enumerate(plot_range):
            try:

                emin, emax = plot_range[i]
            except TypeError:
                q = [0.5 - 0.5 * plot_range[i], 0.5 + 0.5 * plot_range[i]]
                plot_range[i] = quantile(xs[i], q, weights=weights)
                pass
            pass

        pass

    if len(plot_range) != xs.shape[0]:
        raise ValueError(
            "Dimension mismatch between samples and range"
        ) from None

    # Parse the bin specifications.
    try:
        bins = [int(bins) for _ in plot_range]
    except TypeError:
        if len(bins) != len(plot_range):
            raise ValueError(
                "Dimension mismatch between bins and range"
            ) from None
    try:
        hist_bin_factor = [float(hist_bin_factor) for _ in plot_range]
    except TypeError:
        if len(hist_bin_factor) != len(plot_range):
            raise ValueError(
                "Dimension mismatch between hist_bin_factor and plot_range"
            ) from None

    # Some magic numbers for pretty axis layout.
    K = len(xs)
    factor = 2.0  # size of one side of one panel
    if reverse:
        lbdim = 0.2 * factor  # size of left/bottom margin
        trdim = 0.5 * factor  # size of top/right margin
    else:
        lbdim = 0.5 * factor  # size of left/bottom margin
        trdim = 0.2 * factor  # size of top/right margin
    whspace = 0.05  # w/hspace size
    plotdim = factor * K + factor * (K - 1.0) * whspace
    dim = lbdim + plotdim + trdim

    # Create a new figure if one wasn't provided.
    if fig is None:
        fig, axes = plt.subplots(K, K, figsize=(dim, dim))
    else:
        try:
            axes = np.array(fig.axes).reshape((K, K))
        except:  # fmt: skip # noqa: E722 # pylint: disable=bare-except # because we do not expect exceptions and do not want to crash
            raise ValueError(
                "Provided figure has {0} axes, but data has "
                "dimensions K={1}".format(len(fig.axes), K)
            ) from None

    # Format the figure.
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(
        left=lb, bottom=lb, right=tr, top=tr, wspace=whspace, hspace=whspace
    )

    # Set up the default histogram keywords.
    if hist_kwargs is None:
        hist_kwargs = {}
    hist_kwargs["color"] = hist_kwargs.get("color", color)
    if smooth1d is None:
        hist_kwargs["histtype"] = hist_kwargs.get("histtype", "step")

    for i, x in enumerate(xs):
        # Deal with masked arrays.
        if hasattr(x, "compressed"):
            x = x.compressed()

        if np.shape(xs)[0] == 1:
            ax = axes
        else:
            if reverse:
                ax = axes[K - i - 1, K - i - 1]
            else:
                ax = axes[i, i]
        # Plot the histograms.
        if smooth1d is None:
            bins_1d = int(max(1, np.round(hist_bin_factor[i] * bins[i])))
            n, _, _ = ax.hist(
                x,
                bins=bins_1d,
                weights=weights,
                range=np.sort(plot_range[i]),
                **hist_kwargs,
            )
        else:
            if gaussian_filter is None:
                raise ImportError("Please install scipy for smoothing")
            n, b = np.histogram(
                x, bins=bins[i], weights=weights, range=np.sort(plot_range[i])
            )
            n = gaussian_filter(n, smooth1d)
            x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y0 = np.array(list(zip(n, n))).flatten()
            ax.plot(x0, y0, **hist_kwargs)

        if truths is not None and truths[i] is not None:
            ax.axvline(truths[i], color=truth_color)

        # Plot quantiles if wanted.
        if len(quantiles) > 0:
            qvalues = quantile(x, quantiles, weights=weights)
            for q in qvalues:
                ax.axvline(q, ls="dashed", color=color)

            # if verbose:
            #     print("Quantiles:")
            #     print([item for item in zip(quantiles, qvalues)])

        if len(titles) > 0:
            ax.set_title(titles[i], **title_kwargs)

        # Set up the axes.
        ax.set_xlim(plot_range[i])
        if scale_hist:
            maxn = np.max(n)
            ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
        else:
            ax.set_ylim(0, 1.1 * np.max(n))
        ax.set_yticklabels([])
        if max_n_ticks == 0:
            ax.xaxis.set_major_locator(NullLocator())
            ax.yaxis.set_major_locator(NullLocator())
        else:
            ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
            ax.yaxis.set_major_locator(NullLocator())

        if i < K - 1:
            if top_ticks:
                ax.xaxis.set_ticks_position("top")
                for ticklabel in ax.get_xticklabels():
                    ticklabel.set_rotation(45)
            else:
                ax.set_xticklabels([])
        else:
            if reverse:
                ax.xaxis.tick_top()
            for ticklabel in ax.get_xticklabels():
                ticklabel.set_rotation(45)
            if labels is not None:
                if reverse:
                    ax.set_title(labels[i], y=1.25, **label_kwargs)
                else:
                    ax.set_xlabel(labels[i], **label_kwargs)

            # use MathText for axes ticks
            ax.xaxis.set_major_formatter(
                ScalarFormatter(useMathText=use_math_text)
            )

        for j, y in enumerate(xs):
            if np.shape(xs)[0] == 1:
                ax = axes
            else:
                if reverse:
                    ax = axes[K - i - 1, K - j - 1]
                else:
                    ax = axes[i, j]
            if j > i:
                ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            if j == i:
                continue

            # Deal with masked arrays.
            if hasattr(y, "compressed"):
                y = y.compressed()

            hist2d(
                y,
                x,
                ax=ax,
                plot_range=[plot_range[j], plot_range[i]],
                levels=levels,
                **hist2d_kwargs,
            )

            if truths is not None:
                if truths[i] is not None and truths[j] is not None:
                    ax.plot(truths[j], truths[i], "s", color=truth_color)
                if truths[j] is not None:
                    ax.axvline(truths[j], color=truth_color)
                if truths[i] is not None:
                    ax.axhline(truths[i], color=truth_color)

            if max_n_ticks == 0:
                ax.xaxis.set_major_locator(NullLocator())
                ax.yaxis.set_major_locator(NullLocator())
            else:
                ax.xaxis.set_major_locator(
                    MaxNLocator(max_n_ticks, prune="lower")
                )
                ax.yaxis.set_major_locator(
                    MaxNLocator(max_n_ticks, prune="lower")
                )

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                if reverse:
                    ax.xaxis.tick_top()
                for ticklabel in ax.get_xticklabels():
                    ticklabel.set_rotation(45)
                if labels is not None:
                    ax.set_xlabel(labels[j], **label_kwargs)
                    if reverse:
                        ax.xaxis.set_label_coords(0.5, 1.4)
                    else:
                        ax.xaxis.set_label_coords(0.5, -0.3)

                # use MathText for axes ticks
                ax.xaxis.set_major_formatter(
                    ScalarFormatter(useMathText=use_math_text)
                )

            if j > 0:
                ax.set_yticklabels([])
            else:
                if reverse:
                    ax.yaxis.tick_right()
                for ticklabel in ax.get_yticklabels():
                    ticklabel.set_rotation(45)
                if labels is not None:
                    if reverse:
                        ax.set_ylabel(labels[i], rotation=-90, **label_kwargs)
                        ax.yaxis.set_label_coords(1.3, 0.5)
                    else:
                        ax.set_ylabel(labels[i], **label_kwargs)
                        ax.yaxis.set_label_coords(-0.3, 0.5)

                # use MathText for axes ticks
                ax.yaxis.set_major_formatter(
                    ScalarFormatter(useMathText=use_math_text)
                )

    return fig


def quantile(x, q, weights=None):
    """
    Compute sample quantiles with support for weighted samples.

    Note
    ----
    When ``weights`` is ``None``, this method simply calls numpy's percentile
    function with the values of ``q`` multiplied by 100.

    Parameters
    ----------
    x : array_like[nsamples,]
       The samples.

    q : array_like[nquantiles,]
       The list of quantiles to compute. These should all be in the range
       ``[0, 1]``.

    weights : Optional[array_like[nsamples,]]
        An optional weight corresponding to each sample. These

    Returns
    -------
    quantiles : array_like[nquantiles,]
        The sample quantiles computed at ``q``.

    Raises
    ------
    ValueError
        For invalid quantiles; ``q`` not in ``[0, 1]`` or dimension mismatch
        between ``x`` and ``weights``.

    """
    x = np.atleast_1d(x)
    q = np.atleast_1d(q)

    if np.any(q < 0.0) or np.any(q > 1.0):
        raise ValueError("Quantiles must be between 0 and 1")

    if weights is None:
        return np.percentile(x, list(100.0 * q))

    weights = np.atleast_1d(weights)
    if len(x) != len(weights):
        raise ValueError("Dimension mismatch: len(weights) != len(x)")
    idx = np.argsort(x)
    sw = weights[idx]
    cdf = np.cumsum(sw)[:-1]
    cdf /= cdf[-1]
    cdf = np.append(0, cdf)
    return np.interp(q, cdf, x[idx]).tolist()


def hist2d(
    x,
    y,
    plot_range=None,
    levels=None,
    ax=None,
    plot_datapoints=True,
    plot_contours=True,
    contour_kwargs=None,
    data_kwargs=None,
):
    '''Plot a 2-D histogram of samples for plotting posteriors'''
    if ax is None:
        ax = plt.gca()

    if plot_datapoints:
        if data_kwargs is None:
            data_kwargs = {}
        data_kwargs["s"] = data_kwargs.get("s", 2.0)
        data_kwargs["alpha"] = data_kwargs.get("alpha", 0.2)
        ax.scatter(x, y, marker="o", zorder=-1, rasterized=True, **data_kwargs)

    # Plot the contour edge colors.
    if plot_contours:
        if contour_kwargs is None:
            contour_kwargs = {}

        # mask data in range + chi2
        maskx = (x > plot_range[0][0]) & (x < plot_range[0][1])
        masky = (y > plot_range[1][0]) & (y < plot_range[1][1])
        mask = maskx & masky & (data_kwargs['c'] < data_kwargs['vmax'] * 1.2)

        try:  # contour
            # approx posterior + smooth
            xg, yg = np.meshgrid(
                np.linspace(x[mask].min(), x[mask].max(), 256),
                np.linspace(y[mask].min(), y[mask].max(), 256),
            )
            cg = griddata(
                np.vstack([x[mask], y[mask]]).T,
                data_kwargs['c'][mask],
                (xg, yg),
                method='nearest',
                rescale=True,
            )
            scg = gaussian_filter(cg, sigma=15)

            ax.contour(
                xg,
                yg,
                scg * np.nanmin(cg) / np.nanmin(scg),
                np.sort(levels),
                **contour_kwargs,
                vmin=data_kwargs['vmin'],
                vmax=data_kwargs['vmax'],
            )

        except Exception as err:
            print(err)
            print("contour plotting failed")

    ax.set_xlim(plot_range[0])
    ax.set_ylim(plot_range[1])


################################################
# PyLightcurve "mini" - https://github.com/ucl-exoplanets/pylightcurve

# coefficients from https://pomax.github.io/bezierinfo/legendre-gauss.html

gauss0 = [
    [1.0000000000000000, -0.5773502691896257],
    [1.0000000000000000, 0.5773502691896257],
]

gauss10 = [
    [0.2955242247147529, -0.1488743389816312],
    [0.2955242247147529, 0.1488743389816312],
    [0.2692667193099963, -0.4333953941292472],
    [0.2692667193099963, 0.4333953941292472],
    [0.2190863625159820, -0.6794095682990244],
    [0.2190863625159820, 0.6794095682990244],
    [0.1494513491505806, -0.8650633666889845],
    [0.1494513491505806, 0.8650633666889845],
    [0.0666713443086881, -0.9739065285171717],
    [0.0666713443086881, 0.9739065285171717],
]

gauss20 = [
    [0.1527533871307258, -0.0765265211334973],
    [0.1527533871307258, 0.0765265211334973],
    [0.1491729864726037, -0.2277858511416451],
    [0.1491729864726037, 0.2277858511416451],
    [0.1420961093183820, -0.3737060887154195],
    [0.1420961093183820, 0.3737060887154195],
    [0.1316886384491766, -0.5108670019508271],
    [0.1316886384491766, 0.5108670019508271],
    [0.1181945319615184, -0.6360536807265150],
    [0.1181945319615184, 0.6360536807265150],
    [0.1019301198172404, -0.7463319064601508],
    [0.1019301198172404, 0.7463319064601508],
    [0.0832767415767048, -0.8391169718222188],
    [0.0832767415767048, 0.8391169718222188],
    [0.0626720483341091, -0.9122344282513259],
    [0.0626720483341091, 0.9122344282513259],
    [0.0406014298003869, -0.9639719272779138],
    [0.0406014298003869, 0.9639719272779138],
    [0.0176140071391521, -0.9931285991850949],
    [0.0176140071391521, 0.9931285991850949],
]

gauss30 = [
    [0.1028526528935588, -0.0514718425553177],
    [0.1028526528935588, 0.0514718425553177],
    [0.1017623897484055, -0.1538699136085835],
    [0.1017623897484055, 0.1538699136085835],
    [0.0995934205867953, -0.2546369261678899],
    [0.0995934205867953, 0.2546369261678899],
    [0.0963687371746443, -0.3527047255308781],
    [0.0963687371746443, 0.3527047255308781],
    [0.0921225222377861, -0.4470337695380892],
    [0.0921225222377861, 0.4470337695380892],
    [0.0868997872010830, -0.5366241481420199],
    [0.0868997872010830, 0.5366241481420199],
    [0.0807558952294202, -0.6205261829892429],
    [0.0807558952294202, 0.6205261829892429],
    [0.0737559747377052, -0.6978504947933158],
    [0.0737559747377052, 0.6978504947933158],
    [0.0659742298821805, -0.7677774321048262],
    [0.0659742298821805, 0.7677774321048262],
    [0.0574931562176191, -0.8295657623827684],
    [0.0574931562176191, 0.8295657623827684],
    [0.0484026728305941, -0.8825605357920527],
    [0.0484026728305941, 0.8825605357920527],
    [0.0387991925696271, -0.9262000474292743],
    [0.0387991925696271, 0.9262000474292743],
    [0.0287847078833234, -0.9600218649683075],
    [0.0287847078833234, 0.9600218649683075],
    [0.0184664683110910, -0.9836681232797472],
    [0.0184664683110910, 0.9836681232797472],
    [0.0079681924961666, -0.9968934840746495],
    [0.0079681924961666, 0.9968934840746495],
]

gauss40 = [
    [0.0775059479784248, -0.0387724175060508],
    [0.0775059479784248, 0.0387724175060508],
    [0.0770398181642480, -0.1160840706752552],
    [0.0770398181642480, 0.1160840706752552],
    [0.0761103619006262, -0.1926975807013711],
    [0.0761103619006262, 0.1926975807013711],
    [0.0747231690579683, -0.2681521850072537],
    [0.0747231690579683, 0.2681521850072537],
    [0.0728865823958041, -0.3419940908257585],
    [0.0728865823958041, 0.3419940908257585],
    [0.0706116473912868, -0.4137792043716050],
    [0.0706116473912868, 0.4137792043716050],
    [0.0679120458152339, -0.4830758016861787],
    [0.0679120458152339, 0.4830758016861787],
    [0.0648040134566010, -0.5494671250951282],
    [0.0648040134566010, 0.5494671250951282],
    [0.0613062424929289, -0.6125538896679802],
    [0.0613062424929289, 0.6125538896679802],
    [0.0574397690993916, -0.6719566846141796],
    [0.0574397690993916, 0.6719566846141796],
    [0.0532278469839368, -0.7273182551899271],
    [0.0532278469839368, 0.7273182551899271],
    [0.0486958076350722, -0.7783056514265194],
    [0.0486958076350722, 0.7783056514265194],
    [0.0438709081856733, -0.8246122308333117],
    [0.0438709081856733, 0.8246122308333117],
    [0.0387821679744720, -0.8659595032122595],
    [0.0387821679744720, 0.8659595032122595],
    [0.0334601952825478, -0.9020988069688743],
    [0.0334601952825478, 0.9020988069688743],
    [0.0279370069800234, -0.9328128082786765],
    [0.0279370069800234, 0.9328128082786765],
    [0.0222458491941670, -0.9579168192137917],
    [0.0222458491941670, 0.9579168192137917],
    [0.0164210583819079, -0.9772599499837743],
    [0.0164210583819079, 0.9772599499837743],
    [0.0104982845311528, -0.9907262386994570],
    [0.0104982845311528, 0.9907262386994570],
    [0.0045212770985332, -0.9982377097105593],
    [0.0045212770985332, 0.9982377097105593],
]

gauss50 = [
    [0.0621766166553473, -0.0310983383271889],
    [0.0621766166553473, 0.0310983383271889],
    [0.0619360674206832, -0.0931747015600861],
    [0.0619360674206832, 0.0931747015600861],
    [0.0614558995903167, -0.1548905899981459],
    [0.0614558995903167, 0.1548905899981459],
    [0.0607379708417702, -0.2160072368760418],
    [0.0607379708417702, 0.2160072368760418],
    [0.0597850587042655, -0.2762881937795320],
    [0.0597850587042655, 0.2762881937795320],
    [0.0586008498132224, -0.3355002454194373],
    [0.0586008498132224, 0.3355002454194373],
    [0.0571899256477284, -0.3934143118975651],
    [0.0571899256477284, 0.3934143118975651],
    [0.0555577448062125, -0.4498063349740388],
    [0.0555577448062125, 0.4498063349740388],
    [0.0537106218889962, -0.5044581449074642],
    [0.0537106218889962, 0.5044581449074642],
    [0.0516557030695811, -0.5571583045146501],
    [0.0516557030695811, 0.5571583045146501],
    [0.0494009384494663, -0.6077029271849502],
    [0.0494009384494663, 0.6077029271849502],
    [0.0469550513039484, -0.6558964656854394],
    [0.0469550513039484, 0.6558964656854394],
    [0.0443275043388033, -0.7015524687068222],
    [0.0443275043388033, 0.7015524687068222],
    [0.0415284630901477, -0.7444943022260685],
    [0.0415284630901477, 0.7444943022260685],
    [0.0385687566125877, -0.7845558329003993],
    [0.0385687566125877, 0.7845558329003993],
    [0.0354598356151462, -0.8215820708593360],
    [0.0354598356151462, 0.8215820708593360],
    [0.0322137282235780, -0.8554297694299461],
    [0.0322137282235780, 0.8554297694299461],
    [0.0288429935805352, -0.8859679795236131],
    [0.0288429935805352, 0.8859679795236131],
    [0.0253606735700124, -0.9130785566557919],
    [0.0253606735700124, 0.9130785566557919],
    [0.0217802431701248, -0.9366566189448780],
    [0.0217802431701248, 0.9366566189448780],
    [0.0181155607134894, -0.9566109552428079],
    [0.0181155607134894, 0.9566109552428079],
    [0.0143808227614856, -0.9728643851066920],
    [0.0143808227614856, 0.9728643851066920],
    [0.0105905483836510, -0.9853540840480058],
    [0.0105905483836510, 0.9853540840480058],
    [0.0067597991957454, -0.9940319694320907],
    [0.0067597991957454, 0.9940319694320907],
    [0.0029086225531551, -0.9988664044200710],
    [0.0029086225531551, 0.9988664044200710],
]

gauss60 = [
    [0.0519078776312206, -0.0259597723012478],
    [0.0519078776312206, 0.0259597723012478],
    [0.0517679431749102, -0.0778093339495366],
    [0.0517679431749102, 0.0778093339495366],
    [0.0514884515009809, -0.1294491353969450],
    [0.0514884515009809, 0.1294491353969450],
    [0.0510701560698556, -0.1807399648734254],
    [0.0510701560698556, 0.1807399648734254],
    [0.0505141845325094, -0.2315435513760293],
    [0.0505141845325094, 0.2315435513760293],
    [0.0498220356905502, -0.2817229374232617],
    [0.0498220356905502, 0.2817229374232617],
    [0.0489955754557568, -0.3311428482684482],
    [0.0489955754557568, 0.3311428482684482],
    [0.0480370318199712, -0.3796700565767980],
    [0.0480370318199712, 0.3796700565767980],
    [0.0469489888489122, -0.4271737415830784],
    [0.0469489888489122, 0.4271737415830784],
    [0.0457343797161145, -0.4735258417617071],
    [0.0457343797161145, 0.4735258417617071],
    [0.0443964787957871, -0.5186014000585697],
    [0.0443964787957871, 0.5186014000585697],
    [0.0429388928359356, -0.5622789007539445],
    [0.0429388928359356, 0.5622789007539445],
    [0.0413655512355848, -0.6044405970485104],
    [0.0413655512355848, 0.6044405970485104],
    [0.0396806954523808, -0.6449728284894770],
    [0.0396806954523808, 0.6449728284894770],
    [0.0378888675692434, -0.6837663273813555],
    [0.0378888675692434, 0.6837663273813555],
    [0.0359948980510845, -0.7207165133557304],
    [0.0359948980510845, 0.7207165133557304],
    [0.0340038927249464, -0.7557237753065856],
    [0.0340038927249464, 0.7557237753065856],
    [0.0319212190192963, -0.7886937399322641],
    [0.0319212190192963, 0.7886937399322641],
    [0.0297524915007889, -0.8195375261621458],
    [0.0297524915007889, 0.8195375261621458],
    [0.0275035567499248, -0.8481719847859296],
    [0.0275035567499248, 0.8481719847859296],
    [0.0251804776215212, -0.8745199226468983],
    [0.0251804776215212, 0.8745199226468983],
    [0.0227895169439978, -0.8985103108100460],
    [0.0227895169439978, 0.8985103108100460],
    [0.0203371207294573, -0.9200784761776275],
    [0.0203371207294573, 0.9200784761776275],
    [0.0178299010142077, -0.9391662761164232],
    [0.0178299010142077, 0.9391662761164232],
    [0.0152746185967848, -0.9557222558399961],
    [0.0152746185967848, 0.9557222558399961],
    [0.0126781664768160, -0.9697017887650528],
    [0.0126781664768160, 0.9697017887650528],
    [0.0100475571822880, -0.9810672017525982],
    [0.0100475571822880, 0.9810672017525982],
    [0.0073899311633455, -0.9897878952222218],
    [0.0073899311633455, 0.9897878952222218],
    [0.0047127299269536, -0.9958405251188381],
    [0.0047127299269536, 0.9958405251188381],
    [0.0020268119688738, -0.9992101232274361],
    [0.0020268119688738, 0.9992101232274361],
]

gauss_table = [
    np.swapaxes(gauss0, 0, 1),
    np.swapaxes(gauss10, 0, 1),
    np.swapaxes(gauss20, 0, 1),
    np.swapaxes(gauss30, 0, 1),
    np.swapaxes(gauss40, 0, 1),
    np.swapaxes(gauss50, 0, 1),
    np.swapaxes(gauss60, 0, 1),
]


def gauss_numerical_integration(f, x1, x2, precision, *f_args):
    '''gauss_numerical_integration ds'''
    x1, x2 = (x2 - x1) / 2, (x2 + x1) / 2

    return x1 * np.sum(
        gauss_table[precision][0][:, None]
        * f(
            x1[None, :] * gauss_table[precision][1][:, None] + x2[None, :],
            *f_args,
        ),
        0,
    )


def sample_function(f, precision=3):
    '''sample_function ds'''

    def sampled_function(x12_array, *args):
        '''sampled_function ds'''
        x1_array, x2_array = x12_array

        return gauss_numerical_integration(
            f, x1_array, x2_array, precision, *list(args)
        )

    return sampled_function


# orbit
def planet_orbit(
    period,
    sma_over_rs,
    eccentricity,
    inclination,
    periastron,
    mid_time,
    time_array,
    ww=0,
):
    '''planet_orbit ds'''

    inclination = inclination * np.pi / 180.0
    periastron = periastron * np.pi / 180.0
    ww = ww * np.pi / 180.0

    if eccentricity == 0 and ww == 0:
        vv = 2 * np.pi * (time_array - mid_time) / period
        bb = sma_over_rs * np.cos(vv)
        return [
            bb * np.sin(inclination),
            sma_over_rs * np.sin(vv),
            -bb * np.cos(inclination),
        ]

    if periastron < np.pi / 2:
        aa = 1.0 * np.pi / 2 - periastron
    else:
        aa = 5.0 * np.pi / 2 - periastron
    bb = 2 * np.arctan(
        np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(aa / 2)
    )
    if bb < 0:
        bb += 2 * np.pi
    mid_time = float(mid_time) - (period / 2.0 / np.pi) * (
        bb - eccentricity * np.sin(bb)
    )
    m = (
        (
            time_array
            - mid_time
            - np.int_((time_array - mid_time) / period) * period
        )
        * 2.0
        * np.pi
        / period
    )
    u0 = m
    stop = False
    u1 = 0
    for _ in range(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - eccentricity * np.sin(u0) - m) / (
            1 - eccentricity * np.cos(u0)
        )
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        u0 = u1
        pass

    if not stop:
        raise RuntimeError('Failed to find a solution in 10000 loops')

    vv = 2 * np.arctan(
        np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(u1 / 2)
    )
    #
    rr = (
        sma_over_rs
        * (1 - (eccentricity**2))
        / (np.ones_like(vv) + eccentricity * np.cos(vv))
    )
    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)
    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))

    return [x, y, z]


def planet_star_projected_distance(
    period,
    sma_over_rs,
    eccentricity,
    inclination,
    periastron,
    mid_time,
    time_array,
):
    '''planet_star_projected_distance ds'''
    position_vector = planet_orbit(
        period,
        sma_over_rs,
        eccentricity,
        inclination,
        periastron,
        mid_time,
        time_array,
    )

    return np.sqrt(
        position_vector[1] * position_vector[1]
        + position_vector[2] * position_vector[2]
    )


def planet_phase(period, mid_time, time_array):
    '''planet_phase ds'''
    return (time_array - mid_time) / period


# flux drop


def integral_r_claret(limb_darkening_coefficients, r):
    '''integral_r_claret ds'''
    a1, a2, a3, a4 = limb_darkening_coefficients
    mu44 = 1.0 - r * r
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return (
        -(2.0 * (1.0 - a1 - a2 - a3 - a4) / 4) * mu44
        - (2.0 * a1 / 5) * mu44 * mu14
        - (2.0 * a2 / 6) * mu44 * mu24
        - (2.0 * a3 / 7) * mu44 * mu24 * mu14
        - (2.0 * a4 / 8) * mu44 * mu44
    )


def num_claret(r, limb_darkening_coefficients, rprs, z):
    '''num_claret ds'''

    a1, a2, a3, a4 = limb_darkening_coefficients
    rsq = r * r
    mu44 = 1.0 - rsq
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return (
        (
            (1.0 - a1 - a2 - a3 - a4)
            + a1 * mu14
            + a2 * mu24
            + a3 * mu24 * mu14
            + a4 * mu44
        )
        * r
        * np.arccos(np.minimum((-(rprs**2) + z * z + rsq) / (2.0 * z * r), 1.0))
    )


def integral_r_f_claret(
    limb_darkening_coefficients, rprs, z, r1, r2, precision=3
):
    '''integral_r_f_claret ds'''
    return gauss_numerical_integration(
        num_claret, r1, r2, precision, limb_darkening_coefficients, rprs, z
    )


def integral_r_zero(_, r):
    '''integral definitions for zero method'''
    musq = 1 - r * r
    return (-1.0 / 6) * musq * 3.0


def num_zero(r, _, rprs, z):
    '''num_zero ds'''

    rsq = r * r
    return r * np.arccos(
        np.minimum((-(rprs**2) + z * z + rsq) / (2.0 * z * r), 1.0)
    )


def integral_r_f_zero(
    limb_darkening_coefficients, rprs, z, r1, r2, precision=3
):
    '''integral_r_f_zero ds'''
    return gauss_numerical_integration(
        num_zero, r1, r2, precision, limb_darkening_coefficients, rprs, z
    )


# integral definitions for linear method
def integral_r_linear(limb_darkening_coefficients, r):
    '''integral_r_linear ds'''
    a1 = limb_darkening_coefficients[0]
    musq = 1 - r * r
    return (-1.0 / 6) * musq * (3.0 + a1 * (-3.0 + 2.0 * np.sqrt(musq)))


def num_linear(r, limb_darkening_coefficients, rprs, z):
    '''num_linear ds'''

    a1 = limb_darkening_coefficients[0]
    rsq = r * r
    return (
        (1.0 - a1 * (1.0 - np.sqrt(1.0 - rsq)))
        * r
        * np.arccos(np.minimum((-(rprs**2) + z * z + rsq) / (2.0 * z * r), 1.0))
    )


def integral_r_f_linear(
    limb_darkening_coefficients, rprs, z, r1, r2, precision=3
):
    '''integral_r_f_linear ds'''
    return gauss_numerical_integration(
        num_linear, r1, r2, precision, limb_darkening_coefficients, rprs, z
    )


# integral definitions for quadratic method


def integral_r_quad(limb_darkening_coefficients, r):
    '''integral_r_quad ds'''
    a1, a2 = limb_darkening_coefficients[:2]
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return (1.0 / 12) * (
        -4.0 * (a1 + 2.0 * a2) * mu * musq
        + 6.0 * (-1 + a1 + a2) * musq
        + 3.0 * a2 * musq * musq
    )


def num_quad(r, limb_darkening_coefficients, rprs, z):
    '''num_quad ds'''

    a1, a2 = limb_darkening_coefficients[:2]
    rsq = r * r
    cc = 1.0 - np.sqrt(1.0 - rsq)
    return (
        (1.0 - a1 * cc - a2 * cc * cc)
        * r
        * np.arccos(np.minimum((-(rprs**2) + z * z + rsq) / (2.0 * z * r), 1.0))
    )


def integral_r_f_quad(
    limb_darkening_coefficients, rprs, z, r1, r2, precision=3
):
    '''integral_r_f_quad ds'''
    return gauss_numerical_integration(
        num_quad, r1, r2, precision, limb_darkening_coefficients, rprs, z
    )


# integral definitions for square root method


def integral_r_sqrt(limb_darkening_coefficients, r):
    '''integral_r_sqrt ds'''
    a1, a2 = limb_darkening_coefficients[:2]
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return (
        (-2.0 / 5) * a2 * np.sqrt(mu)
        - (1.0 / 3) * a1 * mu
        + (1.0 / 2) * (-1 + a1 + a2)
    ) * musq


def num_sqrt(r, limb_darkening_coefficients, rprs, z):
    '''num_sqrt ds'''

    a1, a2 = limb_darkening_coefficients[:2]
    rsq = r * r
    mu = np.sqrt(1.0 - rsq)
    return (
        (1.0 - a1 * (1 - mu) - a2 * (1.0 - np.sqrt(mu)))
        * r
        * np.arccos(np.minimum((-(rprs**2) + z * z + rsq) / (2.0 * z * r), 1.0))
    )


def integral_r_f_sqrt(
    limb_darkening_coefficients, rprs, z, r1, r2, precision=3
):
    '''integral_r_f_sqrt ds'''
    return gauss_numerical_integration(
        num_sqrt, r1, r2, precision, limb_darkening_coefficients, rprs, z
    )


# dictionaries containing the different methods,
# if you define a new method, include the functions in the dictionary as well

integral_r = {
    'claret': integral_r_claret,
    'linear': integral_r_linear,
    'quad': integral_r_quad,
    'sqrt': integral_r_sqrt,
    'zero': integral_r_zero,
}

integral_r_f = {
    'claret': integral_r_f_claret,
    'linear': integral_r_f_linear,
    'quad': integral_r_f_quad,
    'sqrt': integral_r_f_sqrt,
    'zero': integral_r_f_zero,
}


def integral_centred(method, limb_darkening_coefficients, rprs, ww1, ww2):
    '''integral_centred ds'''
    return (
        integral_r[method](limb_darkening_coefficients, rprs)
        - integral_r[method](limb_darkening_coefficients, 0.0)
    ) * np.abs(ww2 - ww1)


def integral_plus_core(
    method, limb_darkening_coefficients, rprs, z, ww1, ww2, precision=3
):
    '''integral_plus_core ds'''

    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) + np.sqrt(
        np.maximum(rprs**2 - (z * np.sin(ww1)) ** 2, 0)
    )
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) + np.sqrt(
        np.maximum(rprs**2 - (z * np.sin(ww2)) ** 2, 0)
    )
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * w2
    partc = integral_r[method](limb_darkening_coefficients, r2) * (-w1)
    partd = integral_r_f[method](
        limb_darkening_coefficients, rprs, z, r1, r2, precision=precision
    )
    return parta + partb + partc + partd


def integral_minus_core(
    method, limb_darkening_coefficients, rprs, z, ww1, ww2, precision=3
):
    '''integral_minus_core ds'''

    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) - np.sqrt(
        np.maximum(rprs**2 - (z * np.sin(ww1)) ** 2, 0)
    )
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) - np.sqrt(
        np.maximum(rprs**2 - (z * np.sin(ww2)) ** 2, 0)
    )
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * (-w1)
    partc = integral_r[method](limb_darkening_coefficients, r2) * w2
    partd = integral_r_f[method](
        limb_darkening_coefficients, rprs, z, r1, r2, precision=precision
    )
    return parta + partb + partc - partd


def transit_flux_drop(
    limb_darkening_coefficients,
    rp_over_rs,
    z_over_rs,
    method='claret',
    precision=3,
):
    '''transit_flux_drop ds'''

    z_over_rs = np.where(z_over_rs < 0, 1.0 + 100.0 * rp_over_rs, z_over_rs)
    z_over_rs = np.maximum(z_over_rs, 10 ** (-10))

    # cases
    zsq = z_over_rs * z_over_rs
    sum_z_rprs = z_over_rs + rp_over_rs
    dif_z_rprs = rp_over_rs - z_over_rs
    sqr_dif_z_rprs = zsq - rp_over_rs**2
    case0 = np.where((z_over_rs == 0) & (rp_over_rs <= 1))
    case1 = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs <= 1))
    casea = np.where(
        (z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs < 1)
    )
    caseb = np.where(
        (z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs > 1)
    )
    case2 = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs <= 1))
    casec = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs > 1))
    case3 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs < 1))
    case4 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs == 1))
    case5 = np.where(
        (z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs < 1)
    )
    case6 = np.where(
        (z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs == 1)
    )
    case7 = np.where(
        (z_over_rs > rp_over_rs)
        & (sum_z_rprs > 1)
        & (sqr_dif_z_rprs > 1)
        & (-1 < dif_z_rprs)
    )
    plus_case = np.concatenate(
        (case1[0], case2[0], case3[0], case4[0], case5[0], casea[0], casec[0])
    )
    minus_case = np.concatenate(
        (case3[0], case4[0], case5[0], case6[0], case7[0])
    )
    star_case = np.concatenate(
        (case5[0], case6[0], case7[0], casea[0], casec[0])
    )

    # cross points
    ph = np.arccos(
        np.clip((1.0 - rp_over_rs**2 + zsq) / (2.0 * z_over_rs), -1, 1)
    )
    theta_1 = np.zeros(len(z_over_rs))
    ph_case = np.concatenate((case5[0], casea[0], casec[0]))
    theta_1[ph_case] = ph[ph_case]
    theta_2 = np.arcsin(np.minimum(rp_over_rs / z_over_rs, 1))
    theta_2[case1] = np.pi
    theta_2[case2] = np.pi / 2.0
    theta_2[casea] = np.pi
    theta_2[casec] = np.pi / 2.0
    theta_2[case7] = ph[case7]

    # flux_upper
    plusflux = np.zeros(len(z_over_rs))
    plusflux[plus_case] = integral_plus_core(
        method,
        limb_darkening_coefficients,
        rp_over_rs,
        z_over_rs[plus_case],
        theta_1[plus_case],
        theta_2[plus_case],
        precision=precision,
    )
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(
            method, limb_darkening_coefficients, rp_over_rs, 0.0, np.pi
        )
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(
            method, limb_darkening_coefficients, 1, 0.0, np.pi
        )

    # flux_lower
    minsflux = np.zeros(len(z_over_rs))
    minsflux[minus_case] = integral_minus_core(
        method,
        limb_darkening_coefficients,
        rp_over_rs,
        z_over_rs[minus_case],
        0.0,
        theta_2[minus_case],
        precision=precision,
    )

    # flux_star
    starflux = np.zeros(len(z_over_rs))
    starflux[star_case] = integral_centred(
        method, limb_darkening_coefficients, 1, 0.0, ph[star_case]
    )

    # flux_total
    total_flux = integral_centred(
        method, limb_darkening_coefficients, 1, 0.0, 2.0 * np.pi
    )

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)


# transit
def pytransit(
    limb_darkening_coefficients,
    rp_over_rs,
    period,
    sma_over_rs,
    eccentricity,
    inclination,
    periastron,
    mid_time,
    time_array,
    method='claret',
    precision=3,
):
    '''pytransit ds'''

    position_vector = planet_orbit(
        period,
        sma_over_rs,
        eccentricity,
        inclination,
        periastron,
        mid_time,
        time_array,
    )

    projected_distance = np.where(
        position_vector[0] < 0,
        1.0 + 5.0 * rp_over_rs,
        np.sqrt(
            position_vector[1] * position_vector[1]
            + position_vector[2] * position_vector[2]
        ),
    )

    return transit_flux_drop(
        limb_darkening_coefficients,
        rp_over_rs,
        projected_distance,
        method=method,
        precision=precision,
    )


def eclipse_mid_time(
    period, sma_over_rs, eccentricity, inclination, periastron, mid_time
):
    '''eclipse_mid_time ds'''
    test_array = np.arange(0, period, 0.001)
    xx, yy, _ = planet_orbit(
        period,
        sma_over_rs,
        eccentricity,
        inclination,
        periastron,
        mid_time,
        test_array + mid_time,
    )

    test1 = np.where(xx < 0)
    yy = yy[test1]
    test_array = test_array[test1]

    aprox = test_array[np.argmin(np.abs(yy))]

    def function_to_fit(_, t):
        '''function_to_fit ds'''
        _, yy, _ = planet_orbit(
            period,
            sma_over_rs,
            eccentricity,
            inclination,
            periastron,
            mid_time,
            np.array(mid_time + t),
        )
        return yy

    popt, *_ = curve_fit(function_to_fit, [0], [0], p0=[aprox])

    return mid_time + popt[0]


def transit(times, values):
    '''transit ds'''
    return pytransit(
        [values['u0'], values['u1'], values['u2'], values['u3']],
        values['rprs'],
        values['per'],
        values['ars'],
        values['ecc'],
        values['inc'],
        values['omega'],
        values['tmid'],
        times,
        method='claret',
        precision=3,
    )


def eclipse(times, values):
    '''eclipse ds'''
    tme = eclipse_mid_time(
        values['per'],
        values['ars'],
        values['ecc'],
        values['inc'],
        values['omega'],
        values['tmid'],
    )
    model = pytransit(
        [0, 0, 0, 0],
        values['rprs'] * values['fpfs'] ** 0.5,
        values['per'],
        values['ars'],
        values['ecc'],
        values['inc'],
        values['omega'] + 180,
        tme,
        times,
        method='claret',
        precision=3,
    )
    return model


def brightness(time, values):
    '''brightness ds'''
    # compute mid-eclipse time
    tme = eclipse_mid_time(
        values['per'],
        values['ars'],
        values['ecc'],
        values['inc'],
        values['omega'],
        values['tmid'],
    )

    # compute phase based on mid-eclipse
    phase = (time - tme) / values['per']

    # brightness amplitude variation
    bdata = (
        values['c1'] * np.cos(2 * np.pi * phase)
        + values['c2'] * np.sin(2 * np.pi * phase)
        + values['c3'] * np.cos(4 * np.pi * phase)
        + values['c4'] * np.sin(4 * np.pi * phase)
    )

    # offset so eclipse is around 1 in norm flux
    c0 = values['fpfs'] * values['rprs'] ** 2 - (values['c1'] + values['c3'])
    return 1 + c0 + bdata


def phasecurve(time, values):
    '''phasecurve ds'''
    # transit time series
    tdata = transit(time, values)

    # eclipse (similar to transit but with 0 Ld and 180+omega)
    edata = eclipse(time, values)

    # brightness variation
    bdata = brightness(time, values)

    # combine models
    pdata = tdata * bdata * edata

    # mask in-eclipse data
    emask = ~np.floor(edata).astype(bool)

    # offset data to be at 1 (todo inspect this line)
    pdata[emask] = edata[emask] + values['fpfs'] * values['rprs'] ** 2
    return pdata


################################################


def test():
    '''test ds'''

    prior = {
        'rprs': 0.02,  # Rp/Rs
        'ars': 14.25,  # a/Rs
        'per': 3.33,  # Period [day]
        'inc': 88.5,  # Inclination [deg]
        'ecc': 0.5,  # Eccentricity
        'omega': 120,  # Arg of periastron
        'tmid': 0.75,  # Time of mid transit [day],
        'a1': 50,  # Airmass coefficients
        'a2': 0.0,  # trend = a1 * np.exp(a2 * airmass)
        'u0': 0.489193,
        'u1': -0.017295,
        'u2': 0.3956781,
        'u3': -0.165573,
        'teff': 5000,
        'tefferr': 50,
        'met': 0,
        'meterr': 0,
        'logg': 3.89,
        'loggerr': 0.01,
    }

    # # example generating LD coefficients
    # from pylightcurve import exotethys
    # u0, u1, u2, u3 = exotethys(prior['logg'], prior['teff'], prior['met'], 'TESS', method='claret',
    #                            stellar_model='phoenix')
    # prior['u0'], prior['u1'], prior['u2'], prior['u3'] = u0, u1, u2, u3

    time = np.linspace(0.7, 0.8, 1000)  # [day]

    # simulate extinction from airmass
    # airmass = 1./np.cos(np.deg2rad(90-alt))
    airmass = np.zeros(time.shape[0])

    # GENERATE NOISY DATA
    data = transit(time, prior) * prior['a1'] * np.exp(prior['a2'] * airmass)
    data += np.random.normal(0, prior['a1'] * 250e-6, len(time))
    dataerr = np.random.normal(300e-6, 50e-6, len(time)) + np.random.normal(
        300e-6, 50e-6, len(time)
    )

    # add bounds for free parameters only
    mybounds = {
        'rprs': [0, 0.1],
        'tmid': [prior['tmid'] - 0.1, prior['tmid'] + 0.1],
        # 'ars': [13, 15],
        'inc': [86, 89],
    }

    myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds, mode='ns')

    for k in myfit.bounds.keys():
        print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

    _, _ = myfit.plot_bestfit()
    plt.tight_layout()
    plt.show()

    _ = myfit.plot_posterior()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    test()
