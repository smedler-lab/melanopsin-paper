import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks as _find_peaks
from scipy.stats import median_abs_deviation

from ca_imaging.models.peak import Peak
from ca_imaging.models.peak_properties import PeakProperties
from ca_imaging.peak_properties import calculate_auc, calculate_decay_time, calculate_fwhm, calculate_rise_time
from ca_imaging.utils import get_time_step


def _estimate_min_prominence(y) -> float:
    return 6 * median_abs_deviation(y)


def find_peaks(
    y: np.ndarray,
    times: np.ndarray,
    distance_t: float,
    window_size_t: float,
    min_prominence: float | None = None,
) -> list[Peak]:
    time_step = get_time_step(times)
    min_prominence = _estimate_min_prominence(y) if min_prominence is None else min_prominence
    window_size = int((window_size_t / time_step))
    distance = int(distance_t / time_step)

    peak_indexes, peak_data = _find_peaks(
        y, prominence=min_prominence, distance=distance, rel_height=0.5, width=1, wlen=window_size
    )

    peaks = []
    for peak_index in peak_indexes:
        start_index = int(max(0, peak_index - (window_size - 1) / 2))
        end_index = int(min(len(y) - 1, peak_index + (window_size - 1) / 2))
        peaks.append(
            Peak(
                y=y[start_index:end_index],
                times=times[start_index:end_index],
                peak_indexes=np.arange(start_index, end_index),
                peak_index=peak_index,
                peak_index_local=peak_index - start_index,
            )
        )
    return peaks


def _measure_peak_properties(df_f0: np.ndarray, time: np.ndarray, peak: Peak) -> PeakProperties:
    fwhm, fwhm_left_it, fwhm_right_it = calculate_fwhm(peak=peak, zero_baseline=True)
    try:
        rise_time, rise_time_left_it, rise_time_right_it, rise_baseline = calculate_rise_time(peak)
    except Exception:
        rise_time, rise_time_left_it, rise_time_right_it, rise_baseline = None, None, None, None
    try:
        decay_time, decay_time_left_it, decay_time_right_it, decay_baseline = calculate_decay_time(peak)
    except Exception:
        decay_time, decay_time_left_it, decay_time_right_it, decay_baseline = None, None, None, None

    if None not in [rise_time_left_it, decay_time_right_it]:
        auc, auc_left_it, auc_right_it = calculate_auc(peak, time_range=(rise_time_left_it, decay_time_right_it))
    else:
        auc, auc_left_it, auc_right_it = None, None, None

    return PeakProperties(
        fwhm=fwhm,
        fwhm_left_it=fwhm_left_it,
        fwhm_right_it=fwhm_right_it,
        rise_time=rise_time,
        rise_time_left_it=rise_time_left_it,
        rise_time_right_it=rise_time_right_it,
        rise_baseline=rise_baseline,
        decay_time=decay_time,
        decay_time_left_it=decay_time_left_it,
        decay_time_right_it=decay_time_right_it,
        decay_baseline=decay_baseline,
        auc=auc,
        auc_left_it=auc_left_it,
        auc_right_it=auc_right_it,
        amplitude=peak.value,
    )
    pass


def measure_peak_properties(y: np.ndarray, times: np.ndarray, peaks: list[Peak]) -> list[PeakProperties]:
    return [_measure_peak_properties(y, times, peak) for peak in peaks]
