import numpy as np
import pandas as pd

from ca_imaging.models.calcium_signal import CalciumSignal
from ca_imaging.utils import float_s_to_timedelta_ms


def calculate_f0(f: np.ndarray, time: np.ndarray, rolling_quantile: float, window_size_t: pd.Timedelta) -> np.ndarray:
    time_in_ms = float_s_to_timedelta_ms(time)
    return (
        pd.Series(f, time_in_ms)
        .rolling(window=pd.Timedelta(window_size_t, "s"), center=True)
        .quantile(rolling_quantile)
        .values
    )


def calculate_df_f0(
    f: np.ndarray, time: np.ndarray, rolling_quantile: float, window_size_t: pd.Timedelta
) -> np.ndarray:
    f0 = calculate_f0(f, time, rolling_quantile, window_size_t)
    return (f - f0) / f0


def calculate_f_f0(f: np.ndarray, time: np.ndarray, rolling_quantile: float, window_size_t: pd.Timedelta) -> np.ndarray:
    f0 = calculate_f0(f, time, rolling_quantile, window_size_t)
    return f / f0


def process_flourescent_signal(
    f: np.ndarray, times: np.ndarray, rolling_quantile: float = 0.1, window_size_t: pd.Timedelta = pd.Timedelta(60, "s")
) -> CalciumSignal:
    f0 = calculate_f0(f, times, rolling_quantile, window_size_t)
    return CalciumSignal(F=f, F0=f0, DF_F0=(f - f0) / f0, F_F0=f / f0, times=times)
