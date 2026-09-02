from dataclasses import dataclass


@dataclass
class PeakProperties:
    fwhm: float
    fwhm_left_it: float
    fwhm_right_it: float
    rise_time: float
    rise_time_left_it: float
    rise_time_right_it: float
    rise_baseline: float
    decay_time: float
    decay_time_left_it: float
    decay_time_right_it: float
    decay_baseline: float
    auc: float
    auc_left_it: float
    auc_right_it: float
    amplitude: float
