import numpy as np
from typing import Tuple
from ca_imaging.models.peak import Peak


def calculate_fwhm(peak: Peak, zero_baseline=True) -> Tuple[float | None, float | None, float | None]:
    half_max = peak.value / 2

    # Find the left base
    for n in range(peak.peak_index_local, -1, -1):
        if peak.y[n] < half_max:
            dy_dt = (peak.y[n + 1] - peak.y[n]) / (peak.times[n + 1] - peak.times[n])
            m = peak.y[n] - dy_dt * peak.times[n]
            left_it = (half_max - m) / dy_dt
            break
    else:
        left_it = None

    # Find the right base
    for n in range(peak.peak_index_local, peak.y.shape[0]):
        if peak.y[n] < half_max:
            dy_dt = (peak.y[n] - peak.y[n - 1]) / (peak.times[n] - peak.times[n - 1])
            m = peak.y[n] - dy_dt * peak.times[n]
            right_it = (half_max - m) / dy_dt
            break
    else:
        right_it = None

    if None in [left_it, right_it]:
        return None, left_it, right_it

    return right_it - left_it, left_it, right_it


def calculate_rise_time(peak: Peak, zero_baseline: bool = True) -> Tuple[float, float, float, float]:
    MAX_RISE_TIME_IN_S = 20
    s = (peak.y - peak.y.min()) / (peak.value - peak.y.min())

    # Find indexes to use for baseline estimation
    baseline_indexes = np.where(peak.times < peak.time - MAX_RISE_TIME_IN_S)[0]
    if len(baseline_indexes) < 2:
        raise ValueError("Not enough data points t calculate baseline")
    baseline_indexes = baseline_indexes[-10:] if len(baseline_indexes) > 10 else baseline_indexes

    # Calculate baseline
    baseline = s[baseline_indexes].mean()

    # Calculate thresholds
    LOWER_THRESHOLD = baseline * 1.1
    UPPER_THRESHOLD = 0.9

    # Calculate indexes inbetween which the rising of the peak starts
    rise_start_index = np.where((s[: peak.peak_index_local + 1] - LOWER_THRESHOLD * 1.1) < 0)[0][-1]
    rise_end_index = (
        rise_start_index + np.where((s[rise_start_index : peak.peak_index_local + 1] - UPPER_THRESHOLD) > 0)[0][0]
    )

    rise_start_time = np.interp(
        LOWER_THRESHOLD,
        s[rise_start_index - 1 : rise_start_index + 1],
        peak.times[rise_start_index - 1 : rise_start_index + 1],
    )
    rise_end_time = np.interp(
        UPPER_THRESHOLD, s[rise_end_index - 1 : rise_end_index + 1], peak.times[rise_end_index - 1 : rise_end_index + 1]
    )

    if rise_end_time > peak.times[peak.peak_index_local]:
        raise ValueError("Weird peak form")
    if rise_end_time < rise_start_time:
        raise ValueError("Weird peak form")

    return rise_end_time - rise_start_time, rise_start_time, rise_end_time, baseline


def calculate_decay_time(peak: Peak, zero_baseline: bool = False) -> Tuple[float, float, float, float]:
    reversed_peak = Peak(
        y=peak.y[::-1],
        times=-peak.times[::-1],
        peak_index_local=peak.y.shape[0] - peak.peak_index_local - 1,
        peak_index=np.nan,
        peak_indexes=peak.peak_indexes[::-1],
    )
    rise_time, rise_start_time, rise_end_time, baseline = calculate_rise_time(reversed_peak)
    return rise_time, -rise_end_time, -rise_start_time, baseline


def calculate_auc(peak: Peak, time_range: Tuple[int, int] = None) -> Tuple[int, int, int]:
    integration_indexes = np.where((peak.times >= time_range[0]) & (peak.times <= time_range[1]))[0]
    integration_times = peak.times[integration_indexes]
    return (
        ((peak.y[integration_indexes][:-1] + peak.y[integration_indexes][1:]) * np.diff(integration_times) / 2).sum(),
        peak.times[integration_indexes[0]],
        peak.times[integration_indexes[-1]],
    )
