from dataclasses import dataclass
import numpy as np


@dataclass
class Peak:
    y: np.ndarray
    times: np.ndarray
    peak_indexes: np.ndarray
    peak_index: int
    peak_index_local: int

    @property
    def value(self) -> float:
        return self.y[self.peak_index_local]

    @property
    def time(self) -> float:
        return self.times[self.peak_index_local]

    def get_interpolated_y(self, t: float) -> float:
        return np.interp(t, self.times, self.y)
