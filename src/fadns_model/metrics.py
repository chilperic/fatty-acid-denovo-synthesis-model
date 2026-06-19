"""Model-fit metrics."""

from __future__ import annotations

import numpy as np


def ssd(observed, predicted) -> float:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    return float(np.sum((observed - predicted) ** 2))


def rmse(observed, predicted) -> float:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    return float(np.sqrt(np.mean((observed - predicted) ** 2)))
