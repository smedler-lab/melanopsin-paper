
import numpy as np


def get_time_step(time: np.ndarray) -> float:
    return np.median(np.diff(time))


def float_s_to_timedelta_ms(time_in_s: np.ndarray) -> np.ndarray:
    return (time_in_s * 1000).astype("timedelta64[ms]")
