'''phasecurve flare_det ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-branches,too-many-statements,too-many-locals

# -- IMPORTS -- ------------------------------------------------------
import copy
import csv
import json
import os
import sys
import types
import numpy as np
import numpy.linalg as np_linalg
import matplotlib.pyplot as plt

_numpy_linalg_linalg = types.ModuleType('numpy.linalg.linalg')
_numpy_linalg_linalg.LinAlgError = np_linalg.LinAlgError
sys.modules.setdefault('numpy.linalg.linalg', _numpy_linalg_linalg)

from altaipony.flarelc import FlareLightCurve
from altaipony.fakeflares import flare_model_mendoza2022 as model
from altaipony.utils import sigma_clip

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
OBSERVATION_GAP_FACTOR = 10.0
VISIT_COMPLETION_FILENAME = 'visit_status.json'


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


def _as_numpy(values):
    if hasattr(values, 'value'):
        values = values.value
    return np.asarray(values)


def _align_lc_to_cadences(lc, cadenceno):
    if 'cadenceno' not in lc.colnames:
        raise ValueError(
            'Cannot realign detrended light curve: missing cadenceno column.'
        )

    target_cadenceno = np.asarray(cadenceno, dtype=int)
    lc_cadenceno = np.asarray(_as_numpy(lc['cadenceno']), dtype=int)

    sort_idx = np.argsort(lc_cadenceno)
    sorted_cadenceno = lc_cadenceno[sort_idx]
    match_idx = np.searchsorted(sorted_cadenceno, target_cadenceno)
    in_bounds = match_idx < len(sorted_cadenceno)
    valid = np.zeros_like(in_bounds)
    valid[in_bounds] = (
        sorted_cadenceno[match_idx[in_bounds]]
        == target_cadenceno[in_bounds]
    )
    if not np.all(valid):
        missing = target_cadenceno[~valid]
        valid_positions = np.flatnonzero(valid)
        edge_only_missing = False
        first_valid = None
        last_valid = None
        if valid_positions.size:
            first_valid = valid_positions[0]
            last_valid = valid_positions[-1]
            edge_only_missing = np.all(
                valid[first_valid : last_valid + 1]
            )
        if not edge_only_missing:
            raise ValueError(
                'Could not realign detrended light curve to all original '
                f'cadences. Missing {len(missing)} cadence(s): '
                f'{missing[:10].tolist()}'
            )

        aligned_idx = np.empty(len(target_cadenceno), dtype=int)
        aligned_idx[valid] = sort_idx[match_idx[valid]]
        aligned_idx[:first_valid] = aligned_idx[first_valid]
        aligned_idx[last_valid + 1 :] = aligned_idx[last_valid]
        aligned_lc = lc[aligned_idx]
        aligned_lc['cadenceno'] = target_cadenceno
        return aligned_lc

    aligned_idx = sort_idx[match_idx]
    return lc[aligned_idx]


def _detrend_savgol_aligned(
    lc,
    og_flux,
    og_flux_err,
    original_cadenceno,
    original_time=None,
    max_sigma=2.5,
    longdecay=6,
    w=121,
    break_tolerance=10,
    **kwargs,
):
    '''
    Altaipony-compatible Savitzky-Golay detrending with explicit cadence
    realignment back to the original light curve.
    '''
    lcn = lc.normalize()
    sigma_mask = sigma_clip(
        lcn.flux,
        max_sigma=max_sigma,
        longdecay=longdecay,
    )
    mask = ~sigma_mask * 1

    reverse_counts = np.zeros_like(lcn.flux, dtype='int')
    for k in range(1, len(lcn.flux)):
        reverse_counts[-k] = mask[-k] * (
            reverse_counts[-(k - 1)] + mask[-k]
        )

    istart_i = np.where(
        (reverse_counts[1:] >= 1)
        & (reverse_counts[:-1] - reverse_counts[1:] < 0)
    )[0] + 1
    istop_i = istart_i + reverse_counts[istart_i] - 1
    candidates = list(zip(istart_i, istop_i))

    fluxold = lcn.flux.copy()
    lcn.flux[mask] = np.nan

    lcrsf = lcn.flatten(
        window_length=w,
        break_tolerance=break_tolerance,
        **kwargs,
    )
    lcrsf.flux[np.isnan(lcrsf.flux)] = np.nanmedian(lcrsf.flux)

    has_cadence_map = (
        'cadenceno' in lcn.colnames and 'cadenceno' in lcrsf.colnames
    )
    if has_cadence_map:
        lcn_cadenceno = np.asarray(_as_numpy(lcn['cadenceno']), dtype=int)
        lcrsf_cadenceno = np.asarray(
            _as_numpy(lcrsf['cadenceno']),
            dtype=int,
        )
        lcrsf_sort_idx = np.argsort(lcrsf_cadenceno)
        sorted_lcrsf_cadenceno = lcrsf_cadenceno[lcrsf_sort_idx]

    for i, j in candidates:
        mask_ij = np.arange(i, j)
        if has_cadence_map:
            mask_cadenceno = lcn_cadenceno[mask_ij]
            mask_match_idx = np.searchsorted(
                sorted_lcrsf_cadenceno,
                mask_cadenceno,
            )
            mask_in_bounds = mask_match_idx < len(sorted_lcrsf_cadenceno)
            mask_valid = np.zeros_like(mask_in_bounds)
            mask_valid[mask_in_bounds] = (
                sorted_lcrsf_cadenceno[mask_match_idx[mask_in_bounds]]
                == mask_cadenceno[mask_in_bounds]
            )
            if not np.any(mask_valid):
                continue
            lcrsf_mask_ij = lcrsf_sort_idx[mask_match_idx[mask_valid]]
            fluxold_mask_ij = mask_ij[mask_valid]
        else:
            mask_valid = mask_ij < len(lcrsf.flux)
            if not np.any(mask_valid):
                continue
            lcrsf_mask_ij = mask_ij[mask_valid]
            fluxold_mask_ij = mask_ij[mask_valid]

        left_idx = max(0, i - 1)
        right_idx = min(len(lcn.flux) - 1, j)

        left_val = (
            fluxold[left_idx]
            if not np.isnan(fluxold[left_idx])
            else np.nanmedian(fluxold)
        )
        right_val = (
            fluxold[right_idx]
            if not np.isnan(fluxold[right_idx])
            else np.nanmedian(fluxold)
        )

        interpolation_ij = np.interp(
            lcn.time.value[fluxold_mask_ij],
            [lcn.time.value[left_idx], lcn.time.value[right_idx]],
            [left_val, right_val],
        )
        interpolation_ij = np.where(
            interpolation_ij == 0, 1e-10, interpolation_ij
        )
        lcrsf.flux[lcrsf_mask_ij] = fluxold[fluxold_mask_ij] / interpolation_ij

    aligned_lcrsf = _align_lc_to_cadences(lcrsf, original_cadenceno)
    if original_time is not None:
        aligned_lcrsf['time'] = original_time
    og_flux_values = _as_numpy(og_flux)
    og_flux_err_values = _as_numpy(og_flux_err)
    aligned_lcrsf.detrended_flux = (
        _as_numpy(aligned_lcrsf.flux) * np.nanmedian(og_flux_values)
    )
    aligned_lcrsf.detrended_flux_err = og_flux_err_values

    if hasattr(og_flux, 'unit'):
        aligned_lcrsf.flux = og_flux_values * og_flux.unit
        aligned_lcrsf.flux_err = og_flux_err_values * og_flux_err.unit
    else:
        aligned_lcrsf.flux = og_flux_values
        aligned_lcrsf.flux_err = og_flux_err_values

    return aligned_lcrsf


def _safe_savgol_detrend(flc):
    try:
        return flc.detrend('savgol')
    except ValueError as exc:
        if 'Length mismatch after filtering' not in str(exc):
            raise

        log.warning(
            'Altaipony savgol detrend hit a cadence-length mismatch; '
            'retrying with local cadence alignment.'
        )
        clean_lc = copy.deepcopy(flc).remove_nans().find_iterative_median()
        og_flux = clean_lc.flux.copy()
        og_flux_err = clean_lc.flux_err.copy()
        original_cadenceno = np.asarray(
            _as_numpy(clean_lc.cadenceno),
            dtype=int,
        )
        original_time = clean_lc.time.copy()
        interp_lc = clean_lc.interpolate_missing_cadences()
        return _detrend_savgol_aligned(
            interp_lc,
            og_flux,
            og_flux_err,
            original_cadenceno,
            original_time,
        )


def _sanitize_label(label):
    sanitized = ''.join(
        char if char.isalnum() or char in ('-', '_') else '_'
        for char in str(label)
    ).strip('_')
    return sanitized or 'unknown'


def _ensure_directory(path):
    os.makedirs(path, exist_ok=True)
    return path


def _visit_index_label(thisvisit, default_idx):
    try:
        return int(thisvisit.get('visit_index', default_idx))
    except (AttributeError, TypeError, ValueError):
        return int(default_idx)


def _write_text_file(path, lines):
    text = lines if isinstance(lines, str) else '\n'.join(lines)
    if text and not text.endswith('\n'):
        text += '\n'
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write(text)


def _write_binary_file(path, payload):
    with open(path, 'wb') as handle:
        handle.write(payload)


def _read_binary_file(path):
    with open(path, 'rb') as handle:
        return handle.read()


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    return value


def _write_json_file(path, payload):
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(_json_ready(payload), handle, indent=2, sort_keys=True)


def _write_json_file_atomic(path, payload):
    directory = os.path.dirname(path) or '.'
    tmp_path = os.path.join(directory, f'.{os.path.basename(path)}.tmp')
    with open(tmp_path, 'w', encoding='utf-8') as handle:
        json.dump(_json_ready(payload), handle, indent=2, sort_keys=True)
    os.replace(tmp_path, path)


def _write_csv_file(path, fieldnames, rows):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, '') for key in fieldnames})


def _save_figure_file(fig, path):
    fig.savefig(path, dpi=200, bbox_inches='tight')


def _visit_output_dir(results_dir, planet, visit_idx):
    return os.path.join(
        results_dir,
        _sanitize_label(planet),
        f'visit_{visit_idx:02d}',
    )


def _visit_completion_path(visit_dir):
    return os.path.join(visit_dir, VISIT_COMPLETION_FILENAME)


def _load_completed_visit_result(visit_dir, planet, visit_idx):
    completion_path = _visit_completion_path(visit_dir)
    if not os.path.exists(completion_path):
        return None

    try:
        with open(completion_path, 'r', encoding='utf-8') as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        log.warning(
            'Could not read visit completion marker for %s visit %s at %s.',
            planet,
            visit_idx,
            completion_path,
            exc_info=True,
        )
        return None

    if payload.get('status') != 'ok':
        return None

    result = payload.get('result')
    if not isinstance(result, dict):
        return None

    if int(result.get('visit', -1)) != int(visit_idx):
        log.warning(
            'Visit completion marker mismatch for %s visit %s at %s.',
            planet,
            visit_idx,
            completion_path,
        )
        return None

    artifact_names = payload.get('artifacts', [])
    if not isinstance(artifact_names, list) or not artifact_names:
        log.warning(
            'Visit completion marker for %s visit %s is missing artifacts.',
            planet,
            visit_idx,
        )
        return None

    missing_artifacts = [
        name
        for name in artifact_names
        if not os.path.exists(os.path.join(visit_dir, name))
    ]
    if missing_artifacts:
        log.warning(
            'Visit completion marker for %s visit %s is missing %s.',
            planet,
            visit_idx,
            ', '.join(missing_artifacts),
        )
        return None

    return result


def _attach_example_visit_plots(visit_result, visit_dir):
    if visit_result.get('n_flares', 0) <= 0:
        return False

    attached = False
    lightcurve_path = os.path.join(visit_dir, 'detrended_lightcurve.png')
    if os.path.exists(lightcurve_path):
        visit_result['plot_lightcurve'] = _read_binary_file(lightcurve_path)
        attached = True

    flare_list = visit_result.get('flares') or []
    first_flare = flare_list[0] if flare_list else None
    if not isinstance(first_flare, dict):
        return attached

    fit_path = os.path.join(visit_dir, 'flare_00_fit.png')
    if os.path.exists(fit_path):
        first_flare['plot_fit'] = _read_binary_file(fit_path)
        attached = True

    posterior_path = os.path.join(visit_dir, 'flare_00_posterior.png')
    if os.path.exists(posterior_path):
        first_flare['plot_posterior'] = _read_binary_file(posterior_path)
        attached = True

    return attached


def _write_visit_completion_marker(
    visit_dir,
    target_name,
    fltr,
    planet,
    visit_result,
    artifact_names,
):
    if not artifact_names:
        return

    payload = {
        'status': 'ok',
        'target': target_name,
        'filter': fltr,
        'planet': planet,
        'visit': int(visit_result['visit']),
        'artifacts': artifact_names,
        'result': visit_result,
    }
    _write_json_file_atomic(_visit_completion_path(visit_dir), payload)


def _format_metric(value, fmt='.6g', suffix=''):
    if value is None or value == '':
        return 'N/A'
    if isinstance(value, (np.integer, int)):
        text = str(int(value))
    elif isinstance(value, (np.floating, float)):
        text = format(float(value), fmt)
    else:
        text = str(value)
    return f'{text}{suffix}'


def _prefix_lines(lines, prefix):
    return [f'{prefix}{line}' if line else '' for line in lines]


def _flare_summary_lines(flare_index, flare, fltr, c_bol=None):
    c_bol_value = flare.get('c_bol', c_bol)
    return [
        f'Flare {flare_index}',
        f'Start time: {_format_metric(flare.get("start"), ".8f")}',
        f'Stop time: {_format_metric(flare.get("stop"), ".8f")}',
        f'tpeak: {_format_metric(flare.get("tpeak"), ".8f")}',
        (
            'Observed duration: '
            f'{_format_metric(flare.get("observed_duration_days"), ".6f")} days '
            f'({_format_metric(flare.get("observed_duration_minutes"), ".1f")} min)'
        ),
        (
            'FWHM from fit: '
            f'{_format_metric(flare.get("fwhm"), ".6f")} days '
            f'({_format_metric(flare.get("fwhm_minutes"), ".1f")} min)'
        ),
        (
            'Equivalent duration: '
            f'{_format_metric(flare.get("ED_seconds"), ".2f")} s'
        ),
        (
            'Peak flare luminosity in band: '
            f'{_format_metric(flare.get("peak_flare_luminosity_w"), ".2e")} W'
        ),
        (
            f'E_band ({fltr}): '
            f'{_format_metric(flare.get("E_band_ergs"), ".2e")} ergs'
        ),
        (
            f'E_bol (C_bol={_format_metric(c_bol_value, ".1f")}): '
            f'{_format_metric(flare.get("E_bol_ergs"), ".2e")} ergs'
        ),
    ]


def _write_visit_summary_files(
    visit_dir,
    target_name,
    fltr,
    planet,
    visit_idx,
    visit_result,
    time_values,
    c_bol=None,
):
    segments, gaps, summary = _compute_observation_segments(time_values)
    visit_status = 'error' if visit_result.get('error') else 'ok'
    visit_row = {
        'planet': planet,
        'visit': visit_idx,
        'status': visit_status,
        'error': visit_result.get('error', ''),
        'n_points': summary.get('n_points', 0),
        'raw_start': summary.get('raw_start', ''),
        'raw_end': summary.get('raw_end', ''),
        'raw_span_days': summary.get('raw_span_days', ''),
        'median_cadence_days': summary.get('median_cadence_days', ''),
        'gap_threshold_days': summary.get('gap_threshold_days', ''),
        'n_segments': len(segments),
        'n_gaps': len(gaps),
        'n_flares': visit_result.get('n_flares', 0),
    }

    visit_lines = [
        f'Target: {target_name}',
        f'Filter: {fltr}',
        f'Planet: {planet}',
        f'Visit: {visit_idx}',
        f'Status: {visit_status}',
        f'Raw start time: {visit_row["raw_start"]}',
        f'Raw end time: {visit_row["raw_end"]}',
        f'Raw span [days]: {visit_row["raw_span_days"]}',
        f'Estimated cadence [days]: {visit_row["median_cadence_days"]}',
        f'Gap threshold [days]: {visit_row["gap_threshold_days"]}',
        f'Observation segments: {visit_row["n_segments"]}',
        f'Gaps: {visit_row["n_gaps"]}',
        f'Detected flares: {visit_row["n_flares"]}',
    ]

    if visit_result.get('error'):
        visit_lines.extend(['', f'Error: {visit_result["error"]}'])
    elif visit_result.get('flares'):
        for flare_index, flare in enumerate(visit_result['flares']):
            visit_lines.extend([''] + _flare_summary_lines(
                flare_index,
                flare,
                fltr,
                c_bol=c_bol,
            ))
            flare_summary_path = os.path.join(
                visit_dir,
                f'flare_{flare_index:02d}_summary.txt',
            )
            _write_text_file(
                flare_summary_path,
                [
                    f'Target: {target_name}',
                    f'Filter: {fltr}',
                    f'Planet: {planet}',
                    f'Visit: {visit_idx}',
                    '',
                ]
                + _flare_summary_lines(
                    flare_index,
                    flare,
                    fltr,
                    c_bol=c_bol,
                ),
            )
    else:
        visit_lines.extend(['', 'No flares detected in this visit.'])

    if gaps:
        visit_lines.append('')
        visit_lines.append('Gaps:')
        for gap in gaps:
            visit_lines.append(
                (
                    f'  gap {gap["gap_index"]}: '
                    f'{gap["gap_start_raw"]} -> {gap["gap_end_raw"]} '
                    f'({gap["gap_duration_days"]} days)'
                )
            )

    _write_text_file(os.path.join(visit_dir, 'visit_summary.txt'), visit_lines)
    return visit_row, segments, gaps


def _write_planet_summary_file(
    results_dir,
    target_name,
    fltr,
    planet,
    visit_results,
    c_bol=None,
):
    planet_dir = _ensure_directory(os.path.join(results_dir, _sanitize_label(planet)))
    flare_count = sum(visit.get('n_flares', 0) for visit in visit_results)
    error_count = sum(1 for visit in visit_results if visit.get('error'))
    lines = [
        f'Target: {target_name}',
        f'Filter: {fltr}',
        f'Planet: {planet}',
        f'Visits processed: {len(visit_results)}',
        f'Visits with errors: {error_count}',
        f'Total flares: {flare_count}',
    ]

    for visit_result in visit_results:
        visit_idx = int(visit_result.get('visit', -1))
        lines.extend(
            [
                '',
                f'Visit {visit_idx}',
                f'Status: {"error" if visit_result.get("error") else "ok"}',
                f'Detected flares: {visit_result.get("n_flares", 0)}',
            ]
        )
        if visit_result.get('error'):
            lines.append(f'Error: {visit_result["error"]}')
            continue
        if not visit_result.get('flares'):
            lines.append('No flares detected in this visit.')
            continue
        for flare_index, flare in enumerate(visit_result['flares']):
            lines.extend([''] + _prefix_lines(
                _flare_summary_lines(flare_index, flare, fltr, c_bol=c_bol),
                '  ',
            ))

    _write_text_file(os.path.join(planet_dir, 'grand_summary.txt'), lines)


def _compute_observation_segments(time_values, gap_factor=OBSERVATION_GAP_FACTOR):
    times = _as_numpy(time_values).astype(float)
    times = times[np.isfinite(times)]
    if not len(times):
        return [], [], {}

    times = np.unique(np.sort(times))
    if len(times) == 1:
        summary = {
            'n_points': 1,
            'raw_start': float(times[0]),
            'raw_end': float(times[0]),
            'raw_span_days': 0.0,
            'median_cadence_days': 0.0,
            'gap_threshold_days': 0.0,
        }
        segment = {
            'segment_index': 0,
            'raw_start': float(times[0]),
            'raw_end': float(times[0]),
            'start': float(times[0]),
            'end': float(times[0]),
            'duration_days': 0.0,
            'n_points': 1,
            'cadence_days': 0.0,
        }
        return [segment], [], summary

    diffs = np.diff(times)
    positive_diffs = diffs[diffs > 0]
    median_cadence = (
        float(np.median(positive_diffs)) if len(positive_diffs) else 0.0
    )
    gap_threshold = (
        median_cadence * gap_factor if median_cadence > 0 else np.inf
    )
    split_indices = np.where(diffs > gap_threshold)[0]
    seg_starts = np.concatenate(([0], split_indices + 1))
    seg_stops = np.concatenate((split_indices + 1, [len(times)]))

    segments = []
    for seg_index, (start_idx, stop_idx) in enumerate(zip(seg_starts, seg_stops)):
        seg_times = times[start_idx:stop_idx]
        seg_diffs = np.diff(seg_times)
        seg_positive_diffs = seg_diffs[seg_diffs > 0]
        seg_cadence = (
            float(np.median(seg_positive_diffs))
            if len(seg_positive_diffs)
            else median_cadence
        )
        if seg_cadence > 0:
            seg_start = float(seg_times[0] - 0.5 * seg_cadence)
            seg_end = float(seg_times[-1] + 0.5 * seg_cadence)
        else:
            seg_start = float(seg_times[0])
            seg_end = float(seg_times[-1])

        segments.append(
            {
                'segment_index': seg_index,
                'raw_start': float(seg_times[0]),
                'raw_end': float(seg_times[-1]),
                'start': seg_start,
                'end': seg_end,
                'duration_days': max(0.0, seg_end - seg_start),
                'n_points': int(len(seg_times)),
                'cadence_days': float(seg_cadence),
            }
        )

    gaps = []
    for gap_index, split_idx in enumerate(split_indices):
        gaps.append(
            {
                'gap_index': gap_index,
                'gap_start_raw': float(times[split_idx]),
                'gap_end_raw': float(times[split_idx + 1]),
                'gap_duration_days': float(times[split_idx + 1] - times[split_idx]),
            }
        )

    summary = {
        'n_points': int(len(times)),
        'raw_start': float(times[0]),
        'raw_end': float(times[-1]),
        'raw_span_days': float(times[-1] - times[0]),
        'median_cadence_days': median_cadence,
        'gap_threshold_days': (
            float(gap_threshold) if np.isfinite(gap_threshold) else 0.0
        ),
    }
    return segments, gaps, summary


def _merge_intervals(intervals):
    clean_intervals = []
    for interval in intervals:
        start = float(interval['start'])
        end = float(interval['end'])
        if end < start:
            start, end = end, start
        clean_intervals.append(
            {
                'start': start,
                'end': end,
                'payload': interval.get('payload', interval),
            }
        )

    if not clean_intervals:
        return []

    clean_intervals.sort(key=lambda item: (item['start'], item['end']))
    merged = []
    current = {
        'start': clean_intervals[0]['start'],
        'end': clean_intervals[0]['end'],
        'members': [clean_intervals[0]['payload']],
    }

    for interval in clean_intervals[1:]:
        if interval['start'] <= current['end']:
            current['end'] = max(current['end'], interval['end'])
            current['members'].append(interval['payload'])
        else:
            current['duration_days'] = current['end'] - current['start']
            merged.append(current)
            current = {
                'start': interval['start'],
                'end': interval['end'],
                'members': [interval['payload']],
            }

    current['duration_days'] = current['end'] - current['start']
    merged.append(current)
    return merged


def _flatten_flare_rows(results):
    flare_rows = []
    error_rows = []
    for planet, visits in results.items():
        if not isinstance(visits, list):
            continue
        for visit_result in visits:
            visit_idx = int(visit_result.get('visit', -1))
            if visit_result.get('error'):
                error_rows.append(
                    {
                        'planet': planet,
                        'visit': visit_idx,
                        'error': visit_result['error'],
                    }
                )
            for flare_index, flare in enumerate(visit_result.get('flares', [])):
                start = float(flare['start'])
                stop = float(flare['stop'])
                flare_rows.append(
                    {
                        'flare_id': (
                            f'{_sanitize_label(planet)}_visit{visit_idx:02d}'
                            f'_flare{flare_index:02d}'
                        ),
                        'planet': planet,
                        'visit': visit_idx,
                        'flare_index': flare_index,
                        'start_time_raw': min(start, stop),
                        'stop_time_raw': max(start, stop),
                        'duration_days': abs(stop - start),
                        'tpeak': flare.get('tpeak', ''),
                        'fwhm': flare.get('fwhm', ''),
                        'ampl': flare.get('ampl', ''),
                        'area': flare.get('area', ''),
                        'err': flare.get('err', ''),
                        'ED_seconds': flare.get('ED_seconds', ''),
                        'E_band_ergs': flare.get('E_band_ergs', ''),
                        'E_bol_ergs': flare.get('E_bol_ergs', ''),
                    }
                )
    return flare_rows, error_rows


def _build_frequency_products(observation_segments, flare_rows):
    observation_groups = _merge_intervals(
        [
            {
                'start': row['segment_start_raw'],
                'end': row['segment_end_raw'],
                'payload': row,
            }
            for row in observation_segments
        ]
    )
    total_observed_time_days = sum(
        group['duration_days'] for group in observation_groups
    )

    flare_groups = _merge_intervals(
        [
            {
                'start': row['start_time_raw'],
                'end': row['stop_time_raw'],
                'payload': row,
            }
            for row in flare_rows
            if int(row.get('include_in_frequency', 1)) == 1
        ]
    )

    annotated_rows = []
    overlap_rows = []
    for group_index, group in enumerate(flare_groups, start=1):
        members = group['members']
        group_size = len(members)
        for row in members:
            annotated = dict(row)
            annotated['overlap_group_id'] = group_index
            annotated['overlap_group_size'] = group_size
            annotated_rows.append(annotated)

        if group_size > 1:
            overlap_rows.append(
                {
                    'overlap_group_id': group_index,
                    'group_start_raw': group['start'],
                    'group_end_raw': group['end'],
                    'group_duration_days': group['duration_days'],
                    'member_ids': ', '.join(
                        member['flare_id'] for member in members
                    ),
                    'member_labels': ', '.join(
                        f"{member['planet']} visit {member['visit']} "
                        f"flare {member['flare_index']}"
                        for member in members
                    ),
                }
            )

    annotated_rows.sort(
        key=lambda row: (
            row['start_time_raw'],
            row['stop_time_raw'],
            row['planet'],
            row['visit'],
            row['flare_index'],
        )
    )

    flare_frequency_per_day = (
        len(flare_groups) / total_observed_time_days
        if total_observed_time_days > 0
        else 0.0
    )
    flare_frequency_per_hour = flare_frequency_per_day / 24.0

    return {
        'observation_groups': observation_groups,
        'flare_groups': flare_groups,
        'annotated_flare_rows': annotated_rows,
        'overlap_rows': overlap_rows,
        'total_observed_time_days': total_observed_time_days,
        'flare_frequency_per_day': flare_frequency_per_day,
        'flare_frequency_per_hour': flare_frequency_per_hour,
    }


def calculate_flare_frequency(whitelight, flare_results):
    '''
    Calculate overlap-aware flare frequency products from saved flare results.
    '''
    whitelight_data = _normalize_whitelight(whitelight)
    observation_rows = []
    for planet, visits_list in whitelight_data.items():
        for visit_idx, thisvisit in enumerate(visits_list):
            visit_label = _visit_index_label(thisvisit, visit_idx)
            segments, _, _ = _compute_observation_segments(thisvisit['time'])
            for segment in segments:
                observation_rows.append(
                    {
                        'planet': planet,
                        'visit': visit_label,
                        'segment_index': segment['segment_index'],
                        'segment_start_raw': segment['start'],
                        'segment_end_raw': segment['end'],
                        'segment_duration_days': segment['duration_days'],
                    }
                )

    flare_rows, _ = _flatten_flare_rows(flare_results)
    for flare_row in flare_rows:
        flare_row['include_in_frequency'] = 1

    return _build_frequency_products(observation_rows, flare_rows)


def _export_results_bundle(
    results_dir,
    target_name,
    fltr,
    whitelight_data,
    results,
    metadata,
    all_flares_fig=None,
):
    if results_dir is None:
        return None

    results_dir = _ensure_directory(results_dir)
    plots_dir = _ensure_directory(os.path.join(results_dir, 'plots'))

    if all_flares_fig is not None:
        _save_figure_file(
            all_flares_fig,
            os.path.join(plots_dir, 'all_flares.png'),
        )

    result_lookup = {}
    for planet, visit_rows in results.items():
        if isinstance(visit_rows, list):
            for visit_row in visit_rows:
                result_lookup[(planet, int(visit_row['visit']))] = visit_row

    observation_rows = []
    gap_rows = []
    visit_rows = []
    for planet, visits_list in whitelight_data.items():
        for visit_idx, thisvisit in enumerate(visits_list):
            visit_label = _visit_index_label(thisvisit, visit_idx)
            visit_result = result_lookup.get((planet, visit_label), {})
            visit_dir = _ensure_directory(
                os.path.join(
                    results_dir,
                    _sanitize_label(planet),
                    f'visit_{visit_label:02d}',
                )
            )
            visit_row, segments, gaps = _write_visit_summary_files(
                visit_dir=visit_dir,
                target_name=target_name,
                fltr=fltr,
                planet=planet,
                visit_idx=visit_label,
                visit_result=visit_result,
                time_values=thisvisit['time'],
                c_bol=metadata.get('c_bol'),
            )
            visit_rows.append(visit_row)

            for segment in segments:
                observation_rows.append(
                    {
                        'planet': planet,
                        'visit': visit_label,
                        'segment_index': segment['segment_index'],
                        'segment_start_raw': segment['start'],
                        'segment_end_raw': segment['end'],
                        'segment_duration_days': segment['duration_days'],
                        'segment_raw_start': segment['raw_start'],
                        'segment_raw_end': segment['raw_end'],
                        'segment_n_points': segment['n_points'],
                        'segment_cadence_days': segment['cadence_days'],
                    }
                )

            for gap in gaps:
                gap_rows.append(
                    {
                        'planet': planet,
                        'visit': visit_label,
                        **gap,
                    }
                )

    for planet, visit_results in results.items():
        if isinstance(visit_results, list):
            _write_planet_summary_file(
                results_dir=results_dir,
                target_name=target_name,
                fltr=fltr,
                planet=planet,
                visit_results=visit_results,
                c_bol=metadata.get('c_bol'),
            )

    flare_rows, error_rows = _flatten_flare_rows(results)
    for flare_row in flare_rows:
        flare_row['include_in_frequency'] = 1
        flare_row['notes'] = ''

    frequency_products = _build_frequency_products(observation_rows, flare_rows)
    group_lookup = {
        row['flare_id']: row for row in frequency_products['annotated_flare_rows']
    }
    for flare_row in flare_rows:
        group_row = group_lookup.get(flare_row['flare_id'], {})
        flare_row['overlap_group_id'] = group_row.get('overlap_group_id', '')
        flare_row['overlap_group_size'] = group_row.get(
            'overlap_group_size', 1
        )

    run_summary_lines = [
        f'Target: {target_name}',
        f'Filter: {fltr}',
        f'Results directory: {results_dir}',
        f'Visits processed: {len(visit_rows)}',
        f'Visits with errors: {len(error_rows)}',
        f'Detected flare rows: {len(flare_rows)}',
        (
            'Unique flare intervals used in frequency: '
            f'{len(frequency_products["flare_groups"])}'
        ),
        (
            'Total observed time with overlaps removed [days]: '
            f'{frequency_products["total_observed_time_days"]:.8f}'
        ),
        (
            'Flare frequency [1/day]: '
            f'{frequency_products["flare_frequency_per_day"]:.8f}'
        ),
        (
            'Flare frequency [1/hour]: '
            f'{frequency_products["flare_frequency_per_hour"]:.8f}'
        ),
        '',
        'Method notes:',
        (
            'Observation time is built from contiguous raw-time segments in '
            'each visit; gaps larger than 10x the visit median cadence are '
            'excluded before system-wide overlaps are removed.'
        ),
        (
            'Flare frequency uses merged flare intervals across all included '
            'detections, so overlapping flare windows count once.'
        ),
    ]

    if frequency_products['overlap_rows']:
        run_summary_lines.extend(['', 'Overlapping flare groups:'])
        for overlap in frequency_products['overlap_rows']:
            run_summary_lines.append(
                (
                    f'  Group {overlap["overlap_group_id"]}: '
                    f'{overlap["group_start_raw"]} -> '
                    f'{overlap["group_end_raw"]} '
                    f'({overlap["group_duration_days"]} days)'
                )
            )
            run_summary_lines.append(f'    {overlap["member_labels"]}')

    _write_text_file(
        os.path.join(results_dir, 'run_summary.txt'),
        run_summary_lines,
    )
    _write_text_file(
        os.path.join(results_dir, 'flare_frequency_summary.txt'),
        run_summary_lines[:10],
    )

    _write_text_file(
        os.path.join(results_dir, 'manual_review_instructions.txt'),
        [
            'Review the flare rows in flare_frequency_curation.csv.',
            (
                'Set include_in_frequency to 0 for any false positives you '
                'want excluded.'
            ),
            (
                'You can then run '
                f'{os.path.join(os.path.expanduser("~"), "recalculate_flare_frequency.py")} '
                'on this results directory to regenerate the flare-frequency '
                'report without rerunning flare detection.'
            ),
        ],
    )

    _write_csv_file(
        os.path.join(results_dir, 'visit_summary.csv'),
        [
            'planet',
            'visit',
            'status',
            'error',
            'n_points',
            'raw_start',
            'raw_end',
            'raw_span_days',
            'median_cadence_days',
            'gap_threshold_days',
            'n_segments',
            'n_gaps',
            'n_flares',
        ],
        visit_rows,
    )
    _write_csv_file(
        os.path.join(results_dir, 'observation_segments.csv'),
        [
            'planet',
            'visit',
            'segment_index',
            'segment_start_raw',
            'segment_end_raw',
            'segment_duration_days',
            'segment_raw_start',
            'segment_raw_end',
            'segment_n_points',
            'segment_cadence_days',
        ],
        observation_rows,
    )
    _write_csv_file(
        os.path.join(results_dir, 'observation_gaps.csv'),
        [
            'planet',
            'visit',
            'gap_index',
            'gap_start_raw',
            'gap_end_raw',
            'gap_duration_days',
        ],
        gap_rows,
    )
    _write_csv_file(
        os.path.join(results_dir, 'flare_frequency_curation.csv'),
        [
            'flare_id',
            'planet',
            'visit',
            'flare_index',
            'start_time_raw',
            'stop_time_raw',
            'duration_days',
            'overlap_group_id',
            'overlap_group_size',
            'tpeak',
            'fwhm',
            'ampl',
            'area',
            'err',
            'ED_seconds',
            'E_band_ergs',
            'E_bol_ergs',
            'include_in_frequency',
            'notes',
        ],
        flare_rows,
    )
    _write_csv_file(
        os.path.join(results_dir, 'flare_overlaps.csv'),
        [
            'overlap_group_id',
            'group_start_raw',
            'group_end_raw',
            'group_duration_days',
            'member_ids',
            'member_labels',
        ],
        frequency_products['overlap_rows'],
    )
    _write_json_file(
        os.path.join(results_dir, 'detection_results.json'),
        {
            'metadata': metadata,
            'results': results,
            'frequency_summary': {
                'total_observed_time_days': (
                    frequency_products['total_observed_time_days']
                ),
                'unique_flare_intervals': len(
                    frequency_products['flare_groups']
                ),
                'flare_frequency_per_day': (
                    frequency_products['flare_frequency_per_day']
                ),
                'flare_frequency_per_hour': (
                    frequency_products['flare_frequency_per_hour']
                ),
            },
        },
    )

    return {
        'results_dir': results_dir,
        'total_observed_time_days': (
            frequency_products['total_observed_time_days']
        ),
        'flare_frequency_per_day': (
            frequency_products['flare_frequency_per_day']
        ),
        'flare_frequency_per_hour': (
            frequency_products['flare_frequency_per_hour']
        ),
        'unique_flare_intervals': len(frequency_products['flare_groups']),
        'detected_flare_rows': len(flare_rows),
    }


def detect_flares(
    whitelight,
    fin,
    fltr,
    out,
    target=None,
    stellar_params=None,
    results_dir=None,
    resume_completed=True,
    force_rerun=False,
    show_plots=False,
    verbose=False,
):
    '''
    Detect and characterize flares in Spitzer phase-curve white-light data.
    '''
    whitelight_data = _normalize_whitelight(whitelight)
    priors = _normalize_priors(fin)
    target_name = _coerce_target_name(target, whitelight)
    if results_dir is not None:
        results_dir = _ensure_directory(results_dir)

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
    metadata = {
        'target': target_name,
        'filter': fltr,
        'results_dir': results_dir,
        'distance_pc': distance_pc,
        'flux_density_mjy': flux_density_mjy,
        'quiescent_luminosity_w': quiescent_luminosity,
        'c_bol': c_bol,
        'stellar_params': stellar_params,
        'resumed_completed_visits': 0,
    }
    out['data']['target'] = target_name
    out['data']['filter'] = fltr
    out['data']['metadata'] = copy.deepcopy(metadata)
    aggregate_plot_enabled = results_dir is not None or show_plots
    all_flares_fig = None
    all_flares_ax = None
    aggregate_flare_count = 0
    if aggregate_plot_enabled:
        all_flares_fig, all_flares_ax = plt.subplots(figsize=(20, 8))
    resumed_visit_count = 0
    example_plots_attached = False

    for planet, visits_list in whitelight_data.items():
        out['data'][planet] = []
        results[planet] = []

        for idx, thisvisit in enumerate(visits_list):
            visit_label = _visit_index_label(thisvisit, idx)
            if verbose:
                print('-----------------------------------------------------')
                print(
                    f'Planet {planet} visit {visit_label} in {target_name} ({fltr})'
                )
                print('-----------------------------------------------------')

            visit_output_dir = None
            visit_artifacts = []
            if results_dir is not None:
                visit_output_dir = _visit_output_dir(
                    results_dir, planet, visit_label
                )
                if resume_completed and not force_rerun:
                    completed_visit_result = _load_completed_visit_result(
                        visit_output_dir,
                        planet,
                        visit_label,
                    )
                    if completed_visit_result is not None:
                        results[planet].append(completed_visit_result)
                        out_visit_result = copy.deepcopy(completed_visit_result)
                        if not example_plots_attached:
                            example_plots_attached = _attach_example_visit_plots(
                                out_visit_result,
                                visit_output_dir,
                            )
                        out['data'][planet].append(out_visit_result)
                        resumed_visit_count += 1
                        if verbose:
                            print(
                                'Skipping '
                                f'{planet} visit {visit_label}: found completed visit '
                                f'marker in {visit_output_dir}'
                            )
                        continue

                visit_output_dir = _ensure_directory(
                    visit_output_dir
                )

            fig = None
            fig2 = None
            try:
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
                flcd = _safe_savgol_detrend(flc)
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

                fig2, fig2_ax = plt.subplots(figsize=(12, 5))
                fig2_ax.plot(ft, ff, 'o', label='Detrended light curve')
                for start, stop in zip(flares['tstart'], flares['tstop']):
                    fig2_ax.axvspan(start, stop, color='green', alpha=0.3)
                fig2_ax.set_title(
                    f'Light Curve for {target_name} {planet} visit {visit_label} ({fltr})'
                )
                fig2_ax.set_xlabel(f'Time - {thisvisit["time"][0]} [days]')
                fig2_ax.set_ylabel('Relative Flux')
                fig2_ax.set_xlim(ft[0], ft[-1])
                fig2_ax.legend()

                flare_results = []
                flare_output_results = []
                nflares = len(flares['tstart'])
                capture_example_plots = (
                    not example_plots_attached and nflares > 0
                )
                visit_plot = (
                    save_plot_tosv(fig2) if capture_example_plots else None
                )
                if visit_output_dir is not None:
                    _save_figure_file(
                        fig2,
                        os.path.join(
                            visit_output_dir,
                            'detrended_lightcurve.png',
                        ),
                    )
                    visit_artifacts.append('detrended_lightcurve.png')
                plt.close(fig2)
                if verbose:
                    print(
                        f'Detected {nflares} flare(s) for {planet} visit {visit_label}'
                    )

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

                    thres_ax.plot(
                        time, fluxd, zorder=0, label='Raw light curve'
                    )
                    plot_threshold(
                        ax=thres_ax,
                        bt=ft,
                        bf=ff,
                        med=med,
                        threshold=threshold,
                        isflares=isflares,
                    )
                    thres_ax.axvspan(
                        start,
                        stop,
                        color='gray',
                        alpha=0.3,
                        label=f'Flare {index}',
                    )
                    thres_ax.set_title(
                        f'Flare {index} in {target_name} {planet} visit {visit_label}'
                    )
                    thres_ax.set_xlabel(
                        f'Time - {thisvisit["time"][0]} [days]'
                    )
                    thres_ax.set_ylabel('Raw Relative Flux')
                    thres_ax.set_xlim(start_buf, stop_buf)
                    thres_ax.set_ylim(min(bf) - buf, max(bf) + buf)
                    thres_ax.legend()

                    mask = (time >= start_buf) & (time <= stop_buf)
                    masked_time = time[mask]
                    masked_flux = fluxd[mask] - 1.0
                    masked_err = fluxd_err[mask]

                    fit_result = fit_flare_model(
                        masked_time,
                        masked_flux,
                        masked_err,
                        model,
                        start_buf,
                        stop_buf,
                        include_posterior_plot=(
                            capture_example_plots and index == 0
                        ),
                    )
                    if isinstance(fit_result, tuple):
                        fin_params, posterior_plot = fit_result
                    else:
                        fin_params = fit_result
                        posterior_plot = None
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
                            f'  FWHM from fit: {fwhm:.4f} days '
                            f'({fwhm_min:.1f} min)'
                        )

                    model_ax.plot(masked_time, masked_flux, label='Data')
                    model_ax.plot(masked_time, fin_flux, label='Model')
                    plot_params(model_ax, tpeak, fwhm, ampl)
                    model_ax.axvspan(
                        start,
                        stop,
                        color='gray',
                        alpha=0.3,
                        label=f'Flare {index}',
                    )
                    model_ax.set_xlim(start_buf, stop_buf)
                    model_ax.set_title('Flare Model Fit')
                    model_ax.set_xlabel('Relative Time [days]')
                    model_ax.set_ylabel('Raw Relative Flux')
                    model_ax.legend()

                    if all_flares_ax is not None:
                        all_flares_ax.plot(masked_time, masked_flux)
                        all_flares_ax.plot(masked_time, fin_flux)
                        all_flares_ax.axvspan(
                            start,
                            stop,
                            color='gray',
                            alpha=0.3,
                            label=f'{planet}{visit_label}.{index}',
                        )
                        aggregate_flare_count += 1

                    res_flux = masked_flux - fin_flux
                    res_ax.scatter(masked_time, res_flux, label='Residuals')
                    plot_params(res_ax, tpeak, fwhm, ampl)
                    res_ax.axvspan(
                        start,
                        stop,
                        color='gray',
                        alpha=0.3,
                        label=f'Flare {index}',
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
                        'observed_duration_days': obs_duration_days,
                        'observed_duration_minutes': obs_duration_min,
                        'fwhm': fwhm,
                        'fwhm_minutes': fwhm_min,
                        'ampl': ampl,
                        'area': area,
                        'err': error,
                        'ED_days': area,
                        'ED_seconds': ed_seconds,
                        'ED_err_days': error,
                    }

                    if (
                        quiescent_luminosity is not None
                        and c_bol is not None
                    ):
                        e_band_joules = ed_seconds * quiescent_luminosity
                        e_band_ergs = e_band_joules * 1.0e7
                        e_bol_ergs = e_band_ergs * c_bol
                        peak_flare_luminosity = ampl * quiescent_luminosity
                        flare_data['quiescent_luminosity_w'] = (
                            quiescent_luminosity
                        )
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
                            print(
                                f'  E_band ({fltr}): {e_band_ergs:.2e} ergs'
                            )
                            print(
                                f'  E_bol (C_bol={c_bol:.1f}): '
                                f'{e_bol_ergs:.2e} ergs'
                            )
                        caption_lines.append(
                            f'Peak flare L: {peak_flare_luminosity:.2e} W'
                        )
                        caption_lines.append(
                            f'E_band: {e_band_ergs:.2e} ergs'
                        )
                        caption_lines.append(
                            f'E_bol: {e_bol_ergs:.2e} ergs '
                            f'(C_bol={c_bol:.1f})'
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
                    fit_plot = None
                    if capture_example_plots and index == 0:
                        fit_plot = save_plot_tosv(fig)
                    if visit_output_dir is not None:
                        _save_figure_file(
                            fig,
                            os.path.join(
                                visit_output_dir,
                                f'flare_{index:02d}_fit.png',
                            ),
                        )
                        visit_artifacts.append(f'flare_{index:02d}_fit.png')
                    plt.close(fig)

                    flare_results.append(flare_data)
                    flare_output_data = copy.deepcopy(flare_data)
                    if fit_plot is not None:
                        flare_output_data['plot_fit'] = fit_plot
                    if posterior_plot is not None:
                        flare_output_data['plot_posterior'] = posterior_plot
                        if visit_output_dir is not None:
                            posterior_filename = (
                                f'flare_{index:02d}_posterior.png'
                            )
                            _write_binary_file(
                                os.path.join(
                                    visit_output_dir,
                                    posterior_filename,
                                ),
                                posterior_plot,
                            )
                            visit_artifacts.append(posterior_filename)
                    flare_output_results.append(flare_output_data)

                visit_result = {
                    'visit': visit_label,
                    'n_flares': nflares,
                    'flares': flare_results,
                }
                results[planet].append(visit_result)
                if verbose:
                    print()

                out_visit_result = copy.deepcopy(visit_result)
                out_visit_result['flares'] = flare_output_results
                if visit_plot is not None:
                    out_visit_result['plot_lightcurve'] = visit_plot
                out['data'][planet].append(out_visit_result)
                if capture_example_plots:
                    example_plots_attached = True
                if visit_output_dir is not None:
                    _write_visit_summary_files(
                        visit_dir=visit_output_dir,
                        target_name=target_name,
                        fltr=fltr,
                        planet=planet,
                        visit_idx=visit_label,
                        visit_result=visit_result,
                        time_values=thisvisit['time'],
                        c_bol=c_bol,
                    )
                    visit_artifacts.append('visit_summary.txt')
                    for flare_index in range(nflares):
                        visit_artifacts.append(
                            f'flare_{flare_index:02d}_summary.txt'
                        )
                    _write_visit_completion_marker(
                        visit_output_dir,
                        target_name,
                        fltr,
                        planet,
                        visit_result,
                        visit_artifacts,
                    )
            except Exception as exc:  # pragma: no cover
                if fig is not None:
                    plt.close(fig)
                if fig2 is not None:
                    plt.close(fig2)
                log.exception(
                    'Visit processing failed for %s %s visit %s; partial outputs may exist.',
                    target_name,
                    planet,
                    idx,
                )
                error_result = {
                    'visit': visit_label,
                    'n_flares': 0,
                    'flares': [],
                    'error': str(exc),
                }
                results[planet].append(error_result)
                out['data'][planet].append(error_result)
                if verbose:
                    print(
                        f'Visit {planet} {visit_label} failed after partial processing: '
                        f'{exc}'
                    )
                    print()
                continue

    if all_flares_ax is not None and aggregate_flare_count:
        all_flares_ax.set_title(f'All flares in {target_name} ({fltr})')
        all_flares_ax.set_xlabel('Relative Time [days]')
        all_flares_ax.set_ylabel('Relative Flux')
        if aggregate_flare_count <= 20:
            all_flares_ax.legend()
        if show_plots:
            plt.show()

    export_all_flares_fig = all_flares_fig
    if resumed_visit_count:
        export_all_flares_fig = None
        log.info(
            'Skipped %s completed visit(s); leaving all_flares.png unchanged '
            'because the aggregate figure cannot be reconstructed from cached '
            'per-visit metadata alone.',
            resumed_visit_count,
        )
        if verbose:
            print(
                'Skipped '
                f'{resumed_visit_count} completed visit(s); not overwriting '
                'plots/all_flares.png during this resumed run.'
            )
    if not aggregate_flare_count:
        export_all_flares_fig = None

    export_summary = _export_results_bundle(
        results_dir=results_dir,
        target_name=target_name,
        fltr=fltr,
        whitelight_data=whitelight_data,
        results=results,
        metadata={
            **metadata,
            'resumed_completed_visits': resumed_visit_count,
        },
        all_flares_fig=export_all_flares_fig,
    )

    metadata['resumed_completed_visits'] = resumed_visit_count
    out['data']['metadata'] = copy.deepcopy(metadata)

    if export_summary is not None:
        print('saved flare results to', export_summary['results_dir'])
        print(
            'flare frequency [1/day]:',
            export_summary['flare_frequency_per_day'],
        )
        print(
            'unique flare intervals used:',
            export_summary['unique_flare_intervals'],
        )
        results['summary'] = export_summary
        out['data']['frequency_summary'] = copy.deepcopy(export_summary)
    else:
        frequency_products = calculate_flare_frequency(whitelight_data, results)
        flare_rows, _ = _flatten_flare_rows(results)
        out['data']['frequency_summary'] = {
            'total_observed_time_days': (
                frequency_products['total_observed_time_days']
            ),
            'flare_frequency_per_day': (
                frequency_products['flare_frequency_per_day']
            ),
            'flare_frequency_per_hour': (
                frequency_products['flare_frequency_per_hour']
            ),
            'unique_flare_intervals': len(
                frequency_products['flare_groups']
            ),
            'detected_flare_rows': len(flare_rows),
        }
        results['summary'] = copy.deepcopy(out['data']['frequency_summary'])

    if all_flares_fig is not None:
        plt.close(all_flares_fig)
    return results
