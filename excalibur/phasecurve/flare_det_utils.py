'''phasecurve flare_det_utils ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=broad-exception-caught
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-branches,too-many-statements,too-many-locals,too-many-nested-blocks
# no-member is for pm.plot_trace which is not in pymc

# -- IMPORTS -- ------------------------------------------------------
import numpy as np
from scipy import integrate
from collections import defaultdict

# import os
# import json
# import matplotlib.pyplot as plt

import pymc as pm
import pytensor.graph as pg
import pytensor.tensor as pt

# from altaipony.flarelc import FlareLightCurve
from excalibur.util import elca

# from excalibur.system import algorithms as sysalg
from excalibur.util.time import time2z
import excalibur.system.core as syscore

# from .notebook_utils import *

# def get_flare_times(bt, flcd, N1, N2, N3, sigma, diff, transits):
#     # set higher threshold for all transit windows
#     N1_arr = np.full_like(bt, N1, dtype=int)
#     for start, stop in transits:
#         start -= 0.05
#         stop += 0.05
#         mask = (bt >= start) & (bt <= stop)
#         N1_arr[mask] = N1 + diff
#     N2_arr = N1_arr - 1

#     # find flares and flare data point categorization relative to threshold
#     flcd, isflares = flcd.find_flares(N1=N1_arr, N2=N2_arr, N3=N3, sigma=sigma,
#                                         addtail=True, tailthreshdiff=1)
#     flares = flcd.flares.sort_values(by="tstart", ascending=True)

#     # combine all isflares arrays from each segment
#     all_isflare = []
#     for (start, stop), isflare in isflares.items():
#         print(f"From {start} to {stop}, there are {len(isflare)} flux values.")
#         all_isflare.extend(isflare)

#     return flares, all_isflare, N1_arr


def get_transits(
    whitelight: dict,
    priors: dict,
) -> dict:
    """
    Creates a JSON file of transit timestamps across all visits for each planet in target star.

    Parameters:
    - whitelight (dict): Data loaded from pipeline for target star.
                         Retrieved by retrieving value for 'data' key, resulting in each planet mapped to its visits.
                         Eg. whitelight = {'b': [visit1, visit2, ...], 'c': [visit1, visit2, ...], ...}
    - priors (dict): prior system info

    Returns:
    - (str) Path of created transit file
    """

    ssc = syscore.ssconstants()
    planets = whitelight.keys()

    # Find transits for each planet's light curves (one light curve per visit)
    all_transits = defaultdict(lambda: defaultdict(list))

    # Iterate through planets in target star
    for planet, visits in whitelight.items():
        # Iterate through visits in each planet
        for idx, thisvisit in enumerate(visits):
            print(f"Processing planet {planet} visit {idx}...")

            z_dict = {}
            sft_dict = {}

            # Extract light curve variables
            time = np.array(thisvisit['time'])
            fluxd = np.array(thisvisit['detrended'])
            time, fluxd, _ = elca.time_bin(time, fluxd, 5.0 / (60 * 24))

            # find separation values for each planet during given time samples
            for p in planets:
                # current planet's parameters
                planet_priors = priors[p]

                inclination = planet_priors['inc']
                tknot = planet_priors['t0']
                sma = planet_priors['sma'] / priors['R*'] / ssc['Rsun/AU']
                period = planet_priors['period']
                ecc = planet_priors['ecc']

                z, sft = time2z(time, inclination, tknot, sma, period, ecc)

                z_dict[p] = z
                sft_dict[p] = sft

            # find timestamps where transit/eclipse occurs
            s_rad = priors['R*']
            rjup_to_rsun = (
                69911 / 695700
            )  # Jupiter mean radius [km] / Sun mean radius [km]

            # find indices where planet is transiting/eclipsing
            timestamps = defaultdict(list)
            for p, z_planet in z_dict.items():
                # calibrate z-lim with planet's impact
                p_rad = priors[p]['rp']
                z_lim = 1 + p_rad / s_rad * rjup_to_rsun

                mask = (np.array(z_planet) <= z_lim) & (
                    np.array(z_planet) >= -z_lim
                )
                timestamps[p] = list(np.where(mask)[0])

            # extract start and stop indices of transits/eclipses
            limits = defaultdict(list)
            for p, timestamps_p in timestamps.items():
                if not timestamps_p:
                    continue

                limits[p].append(timestamps_p[0])
                for i in range(len(timestamps_p) - 1):
                    if timestamps_p[i] + 1 != timestamps_p[i + 1]:
                        limits[p].append(timestamps_p[i])
                        limits[p].append(timestamps_p[i + 1])
                limits[p].append(timestamps_p[-1])

            # check sft values to distinguish transits from eclipses
            transits = defaultdict(list)
            for p, limits_p in limits.items():
                for i in range(0, len(limits_p), 2):
                    lim_start, lim_stop = limits_p[i], limits_p[i + 1]
                    sft_start, sft_stop = (
                        sft_dict[p][lim_start],
                        sft_dict[p][lim_stop],
                    )

                    # sft crosses 0 from negative to positive for transits
                    if sft_start < 0 < sft_stop:
                        transits[p].append((lim_start, lim_stop))

            # record transits
            for p, transit_list in transits.items():
                for start, stop in transit_list:
                    t_start, t_stop = time[start], time[stop]
                    all_transits[planet][idx].append((t_start, t_stop))

            # if verify:
            #    # plot light curve smoothed with 5-min bins
            #    plt.figure(figsize=(12, 5))
            #    plt.plot(time, fluxd)
            #    # shade transits
            #    for p, transit_list in transits.items():
            #        for start, stop in transit_list:
            #            t_start, t_stop = time[start], time[stop]
            #            plt.axvspan(t_start, t_stop, color='red', alpha=0.3)
            #    plt.xlabel(f"Time -{thisvisit['time'][0]} [day]")
            #    plt.ylabel('5-min Relative Flux')
            #    plt.title(
            #        f"Smoothed light curve of visit {idx} to planet {planet}"
            #    )
            #    # plt.savefig(
            #    #    os.path.join(img_dir, f"planet_{planet}_visit_{idx}.png")
            #    # )

    plain_transits = {
        planet: dict(visits.items()) for planet, visits in all_transits.items()
    }
    # 'unnecessary use of a comprehension' error
    #    planet: {idx: windows for idx, windows in visits.items()}
    #    for planet, visits in all_transits.items()

    return plain_transits


def get_flare_times(bt, flcd, N1, N2, N3, sigma, diff, transits):
    """
    Increase the effective threshold inside (padded) transit windows by scaling sigma,
    while keeping N1, N2 as scalars (avoids AltaiPony shape/broadcast issues).
    """
    # Build a mask for transit windows padded by ±0.05
    transit_mask = np.zeros_like(bt, dtype=bool)
    for start, stop in transits:
        transit_mask |= (bt >= start - 0.05) & (bt <= stop + 0.05)

    # Scale sigma in transit windows so it emulates N1 -> N1+diff
    # (since thresholds compare T1 > N1, raising sigma tightens detection equivalently)
    scale = (N1 + diff) / max(N1, 1)  # avoid division by zero
    sigma_adj = sigma.copy()
    sigma_adj[transit_mask] *= scale

    # Call find_flares with scalar N1, N2; sigma carries the local tightening
    res = flcd.find_flares(
        N1=N1, N2=N2, N3=N3, sigma=sigma_adj, addtail=True, tailthreshdiff=1
    )

    # Normalize return style:
    # - newer AltaiPony: returns a FlareLightCurve (with .flares)
    # - older: may return a tuple (flcd_out, isflares_out, ...)
    # isflares_out = None
    if isinstance(res, tuple):
        flcd_out = res[0]
        # if len(res) > 1:
        #    isflares_out = res[1]
    else:
        flcd_out = res  # detections live here

    # --- extract flares table and convert to dict of lists ---
    flares_dict = {"tstart": [], "tstop": []}
    flares_tbl = getattr(flcd_out, "flares", None)

    # Convert astropy Table -> pandas DataFrame if needed/available
    if flares_tbl is not None:
        try:
            # pandas DataFrame already?
            if (
                hasattr(flares_tbl, "columns")
                and ("tstart" in flares_tbl.columns)
                and ("tstop" in flares_tbl.columns)
            ):
                flares_df = flares_tbl.sort_values(by="tstart", ascending=True)
                tstart = np.asarray(flares_df["tstart"].values, dtype=float)
                tstop = np.asarray(flares_df["tstop"].values, dtype=float)
            # astropy Table
            elif (
                hasattr(flares_tbl, "colnames")
                and ("tstart" in flares_tbl.colnames)
                and ("tstop" in flares_tbl.colnames)
            ):
                # try to_pandas if present; otherwise read columns directly
                try:
                    flares_df = flares_tbl.to_pandas()
                    flares_df = flares_df.sort_values(
                        by="tstart", ascending=True
                    )
                    tstart = np.asarray(flares_df["tstart"].values, dtype=float)
                    tstop = np.asarray(flares_df["tstop"].values, dtype=float)
                # 'broad-exception-caught' error
                except Exception:
                    tstart = np.asarray(flares_tbl["tstart"], dtype=float)
                    tstop = np.asarray(flares_tbl["tstop"], dtype=float)
                    if tstart.size:
                        order = np.argsort(tstart)
                        tstart, tstop = tstart[order], tstop[order]
            else:
                tstart = np.array([], dtype=float)
                tstop = np.array([], dtype=float)
        # 'broad-exception-caught' error
        except Exception:
            tstart = np.array([], dtype=float)
            tstop = np.array([], dtype=float)
    else:
        tstart = np.array([], dtype=float)
        tstop = np.array([], dtype=float)

    # Pack into dict
    if tstart.size:
        flares_dict["tstart"] = tstart.tolist()
        flares_dict["tstop"] = tstop.tolist()

    # --- build boolean mask aligned to bt for plotting ---
    isflares_mask = np.zeros_like(bt, dtype=bool)
    for a, b in zip(tstart, tstop):
        if b < a:
            a, b = b, a
        isflares_mask |= (bt >= a) & (bt <= b)

    threshold = np.asarray(flcd_out.it_med) + (N1 * sigma_adj)

    return flares_dict, isflares_mask, threshold


def get_area_under_lc(model, tpeak, fwhm, ampl, start, stop):
    """Computes area under the flare model with integration."""
    area, error = integrate.quad(
        lambda t: model(t, tpeak=tpeak, fwhm=fwhm, ampl=ampl), start, stop
    )
    return area, error


def fit_flare_model(masked_time, masked_flux, masked_err, model, start, stop):
    """
    Fit Mendoza 2022's model to the given flare,
    sampling parameters peak time, full width at half-maximum, and amplitude.
    """

    # --< MODEL >--
    def mycall(tpeak=np.float128, fwhm=np.float128, ampl=np.float128):
        params = {
            'tpeak': tpeak,
            'fwhm': fwhm,
            'ampl': ampl,
        }  # parameters to be retrieved
        out = myfavRT(params)  # model
        return out

    def myfavRT(
        arg={'tpeak': np.float128, 'fwhm': np.float128, 'ampl': np.float128}
    ):  # RT = relative transfer
        """
        A continuous flare template whose shape is defined by the convolution of a Gaussian and double exponential
        and can be parameterized by three parameters: center time (tpeak), FWHM, and ampitude
        """
        tpeak = arg['tpeak']
        fwhm = arg['fwhm']
        ampl = arg['ampl']

        # out = arg['a'] * xdata + arg['b']
        # out = gauss(t=masked_time, tpeak=tpeak, dur=fwhm, ampl=ampl)
        out = model(t=masked_time, tpeak=tpeak, fwhm=fwhm, ampl=ampl)
        return out

    # --< >--

    # --< PYMC THINGS >--
    def LL(arg1):  # log likelihood
        fwm = mycall(*arg1)
        norm = np.log(np.sqrt(2e0 * np.pi)) + np.log(pymcsigma)
        out = -(((data - fwm) / pymcsigma) ** 2) / 2e0 - norm
        return out

    # class faketensor(pg.Op):
    #    def make_node(self, flatargs) -> pg.Apply:
    #        inputs = [pt.as_tensor(a) for a in flatargs]
    #        outputs = [pt.vector()]
    #        return pg.Apply(self, inputs, outputs)
    #    def perform(
    #        self,
    #        node: pg.Apply,
    #        inputs: list[np.ndarray],
    #        outputs: list[list[None]],
    #    ) -> None:
    #        outputs[0][0] = np.asarray(LL(inputs))
    #        return
    #    pass
    # fkt = faketensor()

    # def fakeshell(tensordata, flatargs):  # tensordata not used
    # actually fakeshell isn't used (now that unused likelihood commented out)
    # def fakeshell(_, flatargs):
    #    return fkt(flatargs)

    # --< >--

    # ndata = int(1e4)
    # global xdata
    # xdata = np.arange(ndata)
    global pymcsigma
    pymcsigma = masked_err
    # pymcsigma = np.array(thisvisit['err'][mask])
    # pymcsigma = np.std(masked_flux[-500:])
    # pymcsigma = 1e-1
    # noise = pymcsigma*np.random.normal(size=ndata)
    global data
    # data = xdata*5e-1 + 5e-1 + noise
    data = masked_flux

    chlen = int(5e4)
    with pm.Model():
        # define priors for parameters
        # arg1 = pm.Uniform('slope', lower=-1e1, upper=1e1, shape=3) # priors
        # flatargs = []
        # flatargs.extend(arg1)
        tpeak = pm.Uniform('tpeak', lower=start, upper=stop)
        fwhm = pm.Uniform('fwhm', lower=0, upper=stop - start)
        ampl = pm.Uniform('ampl', lower=0, upper=max(masked_flux))
        # fwhm = pm.HalfNormal('fwhm', lower=0, upper=stop-start)
        # ampl = pm.HalfNormal('ampl', sigma=np.std(masked_flux))

        # flatargs = [tpeak, fwhm, ampl]
        # likelihood = pm.CustomDist(
        #    "likelihood", flatargs, observed=data, logp=fakeshell
        # )
        step = pm.Metropolis()
        trace = pm.sample(
            chlen,
            step=step,
            tune=int(chlen / 2),
            chains=1,
            cores=1,
            compute_convergence_checks=False,
            progressbar=True,
        )

        pm.plot_trace(trace)

    fin_params = {}
    for key in trace['posterior']:
        # get median
        flat_samples = trace['posterior'][key].values.flatten()
        median = np.median(flat_samples)
        fin_params[key] = median

    return fin_params


def plot_params(ax, tpeak, fwhm, ampl):
    '''Plot parameters of the flare model fit on the given axis.'''
    ax.axhline(y=0, color='black', linestyle='--')
    ax.axvline(x=tpeak, color='red', label=f'Peak time={tpeak:.3f}')
    ax.hlines(
        y=0,
        xmin=tpeak - fwhm / 2,
        xmax=tpeak + fwhm / 2,
        color='magenta',
        label=f'Full width half maximum={fwhm:.3f}',
    )
    ax.vlines(
        x=tpeak,
        ymin=0,
        ymax=ampl,
        color='yellow',
        label=f'Amplitude={ampl:.3f}',
    )


def plot_threshold(ax, bt, bf, med, threshold, isflares):
    '''Plot the threshold for flare detection on the given axis.'''
    # plot median and threshold
    ax.plot(bt, med, color='orange', label='Local scatter', zorder=2)
    ax.plot(
        bt, threshold, color='magenta', label='Detection threshold', zorder=3
    )

    # plot flux data points with colors relative to threshold
    above = np.array(isflares)
    ax.scatter(bt[~above], bf[~above], color='black', zorder=4)
    ax.scatter(bt[above], bf[above], color='red', zorder=5)


def find_flare_overlaps(all_flares):
    flare_groups = defaultdict(list)

    for planet, visits in all_flares.items():
        for visit, flares in visits.items():
            for flare, data in flares.items():
                flare_id = f"{planet}{visit}.{flare}"
                start, stop = data['start'], data['stop']
                added = False
                print(f"Checking flare {flare_id}...")

                if not flare_groups:  # empty dictionary
                    print(
                        f"Adding first flare {flare_id} with start {start} and stop {stop}."
                    )
                    flare_groups[(start, stop)].append(flare_id)
                else:
                    for start_iter, stop_iter in flare_groups.keys():
                        # overlap detected
                        if start < stop_iter and stop > start_iter:
                            key = start_iter, stop_iter
                            print(
                                f"Overlap detected between {flare} and {flare_groups[key]}."
                            )

                            # recalibrate start and stop
                            if start < start_iter or stop > stop_iter:
                                new_key = (
                                    min(start, start_iter),
                                    max(stop, stop_iter),
                                )
                                flare_groups[new_key] = flare_groups.pop(key)
                                key = new_key

                            # add overlapping flare
                            flare_groups[key].append(flare_id)
                            added = True
                            break

                    if not added:  # no overlap detected
                        print(
                            f"Adding unique flare {flare_id} with start {start} and stop {stop}."
                        )
                        flare_groups[(start, stop)].append(flare_id)

    return flare_groups


def merge_flares(short, long):
    """
    Takes the union of flares detected in short cadence data with those in long cadence data.

    Parameters:
    - flares_short_cad (DataFrame): Short cadence (1-min) flare times, sorted by start time.
    - flares_long_cad (DataFrame): Long cadence (5-min) flare times, sorted by start time.

    Returns:
    - merged (list of tuples): Adjusted flare start and stop times to maximize overlap.
    """
    merged = []

    print(f"# of short cadence flares: {len(short)}")
    print(f"# of long cadence flares: {len(long)}")

    for long_start, long_stop in zip(long['tstart'], long['tstop']):
        added = False

        # adjust start and stop times to maximize partial overlap with short-cadence detections
        for short_start, short_stop in zip(short['tstart'], short['tstop']):
            if short_start < long_start < short_stop:
                print(
                    f"Detected flare overlap: short ({short_start}, {short_stop}) with long ({long_start}, {long_stop})"
                )
                # adjust start and stop times to maximize overlap
                # ensures that the merged flare covers the full duration of both detections
                merged_start = min(long_start, short_start)
                merged_stop = max(long_stop, short_stop)
                merged.append((merged_start, merged_stop))
                added = True
                break  # only add each long cadence flare once

        # if no overlap with short cadence data, keep long cadence detection as is
        if not added:
            merged.append((long_start, long_stop))

    print(f"Merged flares: {merged}\n")
    return merged
