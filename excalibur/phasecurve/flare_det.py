'''phasecurve flare_det ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-branches,too-many-statements,too-many-locals

# -- IMPORTS -- ------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

from altaipony.flarelc import FlareLightCurve
from altaipony.fakeflares import flare_model_mendoza2022 as model

from excalibur.util import elca

from excalibur.util.plotters import save_plot_tosv

from excalibur.phasecurve.flare_det_utils import (
    fit_flare_model,
    get_area_under_lc,
    get_flare_times,
    get_transits,
    plot_params,
    plot_threshold,
)

import logging

log = logging.getLogger(__name__)

PC_TO_M = 3.08567758149e16
C_LIGHT = 2.99792458e8
JY_TO_W_PER_M2_PER_HZ = 1e-26
H_PLANCK = 6.62607015e-34
K_BOLTZMANN = 1.380649e-23
SEC_PER_DAY = 86400.0
MIN_PER_DAY = 1440.0
TWO_MASS_KS_LAMBDA_M = 2.159e-6
TWO_MASS_KS_ZEROPOINT_JY = 666.7


def calculate_quiescent_luminosity(
    flux_density_mjy,
    distance_pc,
    lambda_c_m=4.493e-6,
    bandwidth_m=1.015e-6,
):
    '''
    Estimate quiescent stellar luminosity in the observing bandpass.
    '''
    distance_m = distance_pc * PC_TO_M
    flux_density_jy = flux_density_mjy * 1e-3
    flux_density = flux_density_jy * JY_TO_W_PER_M2_PER_HZ
    delta_nu_hz = C_LIGHT * bandwidth_m / (lambda_c_m**2)
    luminosity = 4.0 * np.pi * (distance_m**2) * flux_density * delta_nu_hz
    return luminosity


def _coerce_target_name(target, whitelight):
    if target:
        return target
    return whitelight.get('target', 'UNKNOWN_TARGET')


def _normalize_whitelight(whitelight):
    if 'data' in whitelight:
        return whitelight['data']
    return whitelight


def _normalize_priors(fin):
    if 'priors' in fin:
        return fin['priors']
    return fin


def _resolve_distance_pc(priors, stellar_params):
    if (
        stellar_params is not None
        and stellar_params.get('distance_pc') is not None
    ):
        return stellar_params['distance_pc']
    return priors.get('dist')


def _resolve_bandpass(fltr, stellar_params):
    params = stellar_params or {}
    if (
        params.get('lambda_c_m') is not None
        and params.get('bandwidth_m') is not None
    ):
        return params['lambda_c_m'], params['bandwidth_m']

    if 'IR-36' in fltr:
        return 3.551e-6, 0.750e-6
    if 'IR-45' in fltr:
        return 4.493e-6, 1.010e-6

    return 4.493e-6, 1.010e-6


def _planck_fnu(lambda_m, temperature_k):
    nu_hz = C_LIGHT / lambda_m
    exponent = H_PLANCK * nu_hz / (K_BOLTZMANN * temperature_k)
    return (2.0 * H_PLANCK * nu_hz**3 / C_LIGHT**2) / np.expm1(exponent)


def _estimate_flux_density_mjy(priors, lambda_c_m):
    k_mag = priors.get('Kmag')
    temperature_k = priors.get('T*')
    if k_mag is None or temperature_k is None:
        return None

    flux_k_jy = TWO_MASS_KS_ZEROPOINT_JY * 10 ** (-0.4 * k_mag)
    color_ratio = _planck_fnu(lambda_c_m, temperature_k) / _planck_fnu(
        TWO_MASS_KS_LAMBDA_M, temperature_k
    )
    return flux_k_jy * color_ratio * 1e3


def detect_flares(
    whitelight,
    fin,
    fltr,
    out,
    target=None,
    stellar_params=None,
    verbose=False,
):
    '''
    Detect and characterize flares in Spitzer phase-curve white-light data.
    '''
    whitelight_data = _normalize_whitelight(whitelight)
    priors = _normalize_priors(fin)
    target_name = _coerce_target_name(target, whitelight)

    transits = get_transits(
        whitelight=whitelight_data,
        priors=priors,
    )

    quiescent_luminosity = None
    c_bol = None
    params = stellar_params or {}
    distance_pc = _resolve_distance_pc(priors, stellar_params)
    lambda_c_m, bandwidth_m = _resolve_bandpass(fltr, stellar_params)
    flux_density_mjy = params.get('flux_density_mjy')
    flux_estimated = False
    if flux_density_mjy is None:
        flux_density_mjy = _estimate_flux_density_mjy(priors, lambda_c_m)
        flux_estimated = flux_density_mjy is not None

    if flux_density_mjy is not None:
        if distance_pc is None:
            log.warning(
                'Absolute luminosity/energy not computed: no distance_pc found in stellar_params or system priors.'
            )
        else:
            if (
                params.get('distance_pc') is None
                and priors.get('dist') is not None
            ):
                print(
                    f'Using distance from system priors: {distance_pc:.3f} pc'
                )
            if flux_estimated:
                if verbose:
                    print(
                        'Estimated flux_density_mjy from priors '
                        f'(Kmag={priors.get("Kmag")}, T*={priors.get("T*")} K): '
                        f'{flux_density_mjy:.3f} mJy'
                    )
            quiescent_luminosity = calculate_quiescent_luminosity(
                flux_density_mjy=flux_density_mjy,
                distance_pc=distance_pc,
                lambda_c_m=lambda_c_m,
                bandwidth_m=bandwidth_m,
            )
            c_bol = params.get('c_bol', 100.0)
            if verbose:
                print(
                    'Using bandpass: '
                    f'lambda_c={lambda_c_m:.3e} m, bandwidth={bandwidth_m:.3e} m'
                )
                print(
                    'Band quiescent luminosity: '
                    f'{quiescent_luminosity:.2e} W'
                )
    else:
        print(
            'Absolute luminosity/energy not computed: '
            'no flux_density_mjy provided and could not estimate it from priors.'
        )

    results = {}
    all_flares_fig, all_flares_ax = plt.subplots(figsize=(20, 8))
    any_figures = False

    for planet, visits_list in whitelight_data.items():
        out['data'][planet] = []
        results[planet] = []

        for idx, thisvisit in enumerate(visits_list):
            if verbose:
                print('-----------------------------------------------------')
                print(f'Planet {planet} visit {idx} in {target_name} ({fltr})')
                print('-----------------------------------------------------')

            base_time = thisvisit['time'][0]
            time = np.array(thisvisit['time']) - base_time
            fluxd = np.array(thisvisit['detrended'])
            fluxd_err = fluxd * (
                np.array(thisvisit['err']) / np.array(thisvisit['flux'])
            )

            bt, bf, bstds = elca.time_bin(
                time,
                fluxd,
                2.0 / (60 * 24),
            )
            flc = FlareLightCurve(time=bt, flux=bf, flux_err=bstds)
            flcd = flc.detrend('savgol')
            ft = np.asarray(
                flcd.time.jd if len(bt) != len(flcd.time.jd) else bt
            )
            ff = np.asarray(flcd.detrended_flux)
            fstds = np.asarray(flcd.detrended_flux_err)
            med = np.asarray(flcd.it_med)

            visit_transits = np.asarray(
                transits.get(planet, {}).get(idx, []),
                dtype=float,
            )
            if visit_transits.size:
                visit_transits = visit_transits - base_time

            if verbose:
                print('Finding flares...')
            flares, isflares, threshold = get_flare_times(
                bt=ft,
                flcd=flcd,
                N1=3,
                N2=2,
                N3=3,
                sigma=fstds,
                diff=3,
                transits=visit_transits,
            )

            fig2 = plt.figure(figsize=(12, 5))
            plt.plot(ft, ff, 'o', label='Detrended light curve')
            for start, stop in zip(flares['tstart'], flares['tstop']):
                plt.axvspan(start, stop, color='green', alpha=0.3)
            plt.title(
                f'Light Curve for {target_name} {planet} visit {idx} ({fltr})'
            )
            plt.xlabel(f'Time - {thisvisit["time"][0]} [days]')
            plt.ylabel('Relative Flux')
            plt.xlim(ft[0], ft[-1])
            plt.legend()
            if verbose:
                print('YEAH ITS SHOWING2!!')
                plt.show()
                fig2.show()
            plt.close(fig2)
            any_figures = True

            flare_results = []
            nflares = len(flares['tstart'])
            if verbose:
                print(f'Detected {nflares} flare(s) for {planet} visit {idx}')

            for index, (start, stop) in enumerate(
                zip(flares['tstart'], flares['tstop'])
            ):
                if verbose:
                    print(
                        f'Fitting flare {index} from {start:.6f} to {stop:.6f}...'
                    )

                buf = 0.01
                start_buf = start - buf
                stop_buf = stop + buf
                obs_duration_days = stop - start
                obs_duration_min = obs_duration_days * MIN_PER_DAY

                fig, axs = plt.subplots(3, 1, figsize=(15, 20))
                thres_ax, model_ax, res_ax = axs

                thres_ax.plot(time, fluxd, zorder=0, label='Raw light curve')
                plot_threshold(
                    ax=thres_ax,
                    bt=ft,
                    bf=ff,
                    med=med,
                    threshold=threshold,
                    isflares=isflares,
                )
                thres_ax.axvspan(
                    start, stop, color='gray', alpha=0.3, label=f'Flare {index}'
                )
                thres_ax.set_title(
                    f'Flare {index} in {target_name} {planet} visit {idx}'
                )
                thres_ax.set_xlabel(f'Time - {thisvisit["time"][0]} [days]')
                thres_ax.set_ylabel('Raw Relative Flux')
                thres_ax.set_xlim(start_buf, stop_buf)
                thres_ax.set_ylim(min(bf) - buf, max(bf) + buf)
                thres_ax.legend()

                mask = (time >= start_buf) & (time <= stop_buf)
                masked_time = time[mask]
                masked_flux = fluxd[mask] - 1.0
                masked_err = fluxd_err[mask]

                fin_params = fit_flare_model(
                    masked_time,
                    masked_flux,
                    masked_err,
                    model,
                    start_buf,
                    stop_buf,
                )
                tpeak = fin_params['tpeak']
                fwhm = fin_params['fwhm']
                ampl = fin_params['ampl']
                fin_flux = model(
                    t=masked_time, tpeak=tpeak, fwhm=fwhm, ampl=ampl
                )

                fwhm_min = fwhm * MIN_PER_DAY
                if verbose:
                    print(
                        f'  Observed duration: {obs_duration_days:.4f} days '
                        f'({obs_duration_min:.1f} min)'
                    )
                    print(
                        f'  FWHM from fit: {fwhm:.4f} days ({fwhm_min:.1f} min)'
                    )

                model_ax.plot(masked_time, masked_flux, label='Data')
                model_ax.plot(masked_time, fin_flux, label='Model')
                plot_params(model_ax, tpeak, fwhm, ampl)
                model_ax.axvspan(
                    start, stop, color='gray', alpha=0.3, label=f'Flare {index}'
                )
                model_ax.set_xlim(start_buf, stop_buf)
                model_ax.set_title('Flare Model Fit')
                model_ax.set_xlabel('Relative Time [days]')
                model_ax.set_ylabel('Raw Relative Flux')
                model_ax.legend()

                all_flares_ax.plot(masked_time, masked_flux)
                all_flares_ax.plot(masked_time, fin_flux)
                all_flares_ax.axvspan(
                    start,
                    stop,
                    color='gray',
                    alpha=0.3,
                    label=f'{planet}{idx}.{index}',
                )

                res_flux = masked_flux - fin_flux
                res_ax.scatter(masked_time, res_flux, label='Residuals')
                plot_params(res_ax, tpeak, fwhm, ampl)
                res_ax.axvspan(
                    start, stop, color='gray', alpha=0.3, label=f'Flare {index}'
                )
                res_ax.plot(
                    masked_time,
                    masked_err * 2,
                    color='green',
                    label='2 sigma',
                )
                res_ax.plot(masked_time, -masked_err * 2, color='green')
                res_ax.set_title('Residuals')
                res_ax.set_xlabel('Relative Time [days]')
                res_ax.set_ylabel('Raw Residuals')
                res_ax.set_xlim(start_buf, stop_buf)

                area, error = get_area_under_lc(
                    model, tpeak, fwhm, ampl, start_buf, stop_buf
                )
                ed_seconds = area * SEC_PER_DAY
                caption_lines = [
                    f'Integrated Area: {area:.2f} +/- {error:.2f}',
                    f'Equivalent Duration: {ed_seconds:.2f} s',
                ]
                if verbose:
                    print(f' Equivalent duration: {ed_seconds:.2f} s')

                flare_data = {
                    'start': start + base_time,
                    'stop': stop + base_time,
                    'tpeak': tpeak,
                    'fwhm': fwhm,
                    'ampl': ampl,
                    'area': area,
                    'err': error,
                    'ED_days': area,
                    'ED_seconds': ed_seconds,
                    'ED_err_days': error,
                }

                if quiescent_luminosity is not None and c_bol is not None:
                    e_band_joules = ed_seconds * quiescent_luminosity
                    e_band_ergs = e_band_joules * 1.0e7
                    e_bol_ergs = e_band_ergs * c_bol
                    peak_flare_luminosity = ampl * quiescent_luminosity
                    flare_data['quiescent_luminosity_w'] = quiescent_luminosity
                    flare_data['peak_flare_luminosity_w'] = (
                        peak_flare_luminosity
                    )
                    flare_data['E_band_ergs'] = e_band_ergs
                    flare_data['E_bol_ergs'] = e_bol_ergs
                    if verbose:
                        print(
                            '  Peak flare luminosity in band: '
                            f'{peak_flare_luminosity:.2e} W'
                        )
                        print(f'  E_band ({fltr}): {e_band_ergs:.2e} ergs')
                        print(
                            f'  E_bol (C_bol={c_bol:.1f}): {e_bol_ergs:.2e} ergs'
                        )
                    caption_lines.append(
                        f'Peak flare L: {peak_flare_luminosity:.2e} W'
                    )
                    caption_lines.append(f'E_band: {e_band_ergs:.2e} ergs')
                    caption_lines.append(
                        f'E_bol: {e_bol_ergs:.2e} ergs (C_bol={c_bol:.1f})'
                    )
                else:
                    print(
                        '  Absolute luminosity/energy not computed: '
                        'stellar_params not provided.'
                    )
                    caption_lines.append(
                        'Absolute luminosity/energy not computed'
                    )

                res_ax.text(
                    start_buf + 0.001,
                    np.min(res_flux) - 0.001,
                    '\n'.join(caption_lines),
                    fontsize=10,
                    verticalalignment='top',
                )
                res_ax.legend()
                fig.subplots_adjust(hspace=0.4, wspace=0.4)
                if verbose:
                    print('YEAH ITS SHOWING4!!')
                    plt.show()
                    fig.show()
                plt.close(fig)

                any_figures = True

                flare_results.append(flare_data)

            results[planet].append(
                {
                    'visit': idx,
                    'n_flares': nflares,
                    'flares': flare_results,
                }
            )
            if verbose:
                print()

            out['data'][planet].append(
                {
                    'visit': idx,
                    'n_flares': nflares,
                    'flares': flare_results,
                }
            )
            flare_plot = save_plot_tosv(all_flares_fig)
            out['data'][planet][-1]['plot_lightcurve'] = flare_plot

    print('any figures?', any_figures)
    print('show plots?', verbose)
    if any_figures:
        all_flares_ax.set_title(f'All flares in {target_name} ({fltr})')
        all_flares_ax.set_xlabel('Relative Time [days]')
        all_flares_ax.set_ylabel('Relative Flux')
        all_flares_ax.legend()
        if verbose:
            print('YEAH ITS SHOWING!!')
            plt.show()
            input('pause until input1')
            all_flares_fig.show()
            input('pause until input2')

    print('results', results)
    return results
