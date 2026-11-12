"""Read tabular observations without losing the first headerless measurement."""

import csv
from pathlib import Path

import numpy as np
import pandas as pd


def load_fatty_acid_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        first = next(csv.reader(handle), [])
    if not first:
        raise ValueError("Empty CSV")
    try:
        [float(value) for value in first]
        headerless = True
    except ValueError:
        headerless = False
    if headerless:
        if len(first) != 2:
            raise ValueError("Headerless data must contain time and one observation")
        frame = pd.read_csv(path, header=None, names=["time", "value"])
    else:
        frame = pd.read_csv(path)
    values = frame.to_numpy(dtype=float)
    if not values.size or not np.all(np.isfinite(values)):
        raise ValueError("Observations must be nonempty and finite")
    if np.any(np.diff(values[:, 0]) <= 0):
        raise ValueError("Observation times must be strictly increasing")
    frame.attrs["source"] = path.name
    frame.attrs["time_unit"] = "s" if frame.columns[0] == "time(s)" else "unverified"
    frame.attrs["measurement_unit"] = (
        "unverified; consult data metadata and source publication"
    )
    return frame
