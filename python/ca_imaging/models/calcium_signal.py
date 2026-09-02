from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class CalciumSignal:
    F: np.ndarray
    F0: np.ndarray
    F_F0: np.ndarray
    DF_F0: np.ndarray
    times: np.ndarray

    def __getitem__(self, index):
        """Slice the calcium signal data by time index.

        Args:
            index: Integer, slice, or array of indices to select

        Returns:
            CalciumSignal: A new CalciumSignal object with sliced data
        """
        return CalciumSignal(
            F=self.F[index] if isinstance(self.F, np.ndarray) else self.F,
            F0=self.F0[index] if isinstance(self.F0, np.ndarray) else self.F0,
            F_F0=self.F_F0[index] if isinstance(self.F_F0, np.ndarray) else self.F_F0,
            DF_F0=self.DF_F0[index] if isinstance(self.DF_F0, np.ndarray) else self.DF_F0,
            times=self.times[index] if isinstance(self.times, np.ndarray) else self.times,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Convert the calcium signal data to a pandas DataFrame.

        Returns:
            pd.DataFrame: A DataFrame containing the calcium signal data
        """
        return pd.DataFrame(
            {"F": self.F, "F0": self.F0, "F_F0": self.F_F0, "DF_F0": self.DF_F0, "Times (s)": self.times}
        )

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "CalciumSignal":
        """Create a CalciumSignal object from a pandas DataFrame.

        Args:
            df: A DataFrame containing the calcium signal data with columns "F", "F0", "F_F0", "DF_F0", and "Times (s)"

        Returns:
            CalciumSignal: A new CalciumSignal object created from the DataFrame
        """
        return cls(
            F=df["F"].values,
            F0=df["F0"].values,
            F_F0=df["F_F0"].values,
            DF_F0=df["DF_F0"].values,
            times=df.index.values,
        )

    def to_csv(self, file_path: str):
        """Save the calcium signal data to a CSV file.

        Args:
            file_path: The path to the CSV file where the data will be saved
        """
        self.to_dataframe().to_csv(file_path, index=False)

    @classmethod
    def from_csv(cls, file_path: str) -> "CalciumSignal":
        """Create a CalciumSignal object from a CSV file.

        Args:
            file_path: The path to the CSV file containing the calcium signal data
        Returns:
            CalciumSignal: A new CalciumSignal object created from the CSV file
        """
        df = pd.read_csv(file_path)
        return cls.from_dataframe(df)


def calcium_signals_from_csv(
    folder: str,
    F: str | None = None,
    F0: str | None = None,
    F_F0: str | None = None,
    DF_F0: str | None = None,
) -> list[CalciumSignal]:
    """Load a list of csv files with calcium signals into a list of CalciumSignal objects.

    Args:
        F: Path to the CSV file containing the F values
        F0: Path to the CSV file containing the F0 values
        F_F0: Path to the CSV file containing the F/F0 values
        DF_F0: Path to the CSV file containing the dF/F0 values

    Returns:
        list[CalciumSignal]: A list of CalciumSignal objects created from the CSV files
    """
    F = F or Path(folder) / "f.csv"
    F0 = F0 or Path(folder) / "f0.csv"
    F_F0 = F_F0 or Path(folder) / "f_f0.csv"
    DF_F0 = DF_F0 or Path(folder) / "df_f0.csv"

    df_F = pd.read_csv(F, index_col="Time (s)")
    df_F0 = pd.read_csv(F0, index_col="Time (s)")
    df_F_F0 = pd.read_csv(F_F0, index_col="Time (s)")
    df_DF_F0 = pd.read_csv(DF_F0, index_col="Time (s)")
    if np.any(np.vstack([df_F0.index.values, df_F_F0.index.values, df_DF_F0.index.values]) != df_F.index.values):
        raise ValueError("Time columns in the provided DataFrames do not match.")

    if np.any(
        np.vstack(
            [df_F0.columns.values.astype(int), df_F_F0.columns.values.astype(int), df_DF_F0.columns.values.astype(int)]
        )
        != df_F.columns.values.astype(int)
    ):
        raise ValueError("Columns differ between dataframes")

    return [
        CalciumSignal.from_dataframe(
            pd.DataFrame(
                {
                    "F": df_F[col].values,
                    "F0": df_F0[col].values,
                    "F_F0": df_F_F0[col].values,
                    "DF_F0": df_DF_F0[col].values,
                },
                index=df_F.index,
            )
        )
        for col in df_F.columns
    ]


def calcium_signals_to_csv(calcium_signals: list[CalciumSignal], output_folder: str):
    """Save a list of CalciumSignal objects to CSV files.

    Args:
        calcium_signals: A list of CalciumSignal objects to be saved
        output_folder: The directory where the CSV files will be saved
    """
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    times = calcium_signals[0].times
    f = np.array([cs.F for cs in calcium_signals])
    f_f0 = np.array([cs.F_F0 for cs in calcium_signals])
    df_f0 = np.array([cs.DF_F0 for cs in calcium_signals])

    for signal, name in zip([f, f_f0, df_f0], ["f.csv", "f_f0.csv", "df_f0.csv"]):
        df = pd.DataFrame(signal.T)
        df["Time (s)"] = times
        df = df.set_index("Time (s)")
        df.to_csv(output_folder / name)
