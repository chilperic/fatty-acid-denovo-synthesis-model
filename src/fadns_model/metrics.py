"""Strict numerical metrics: missing data must be handled explicitly upstream."""

import numpy as np


def _residuals(observed, predicted):
    observed, predicted = (
        np.asarray(observed, dtype=float),
        np.asarray(predicted, dtype=float),
    )
    if not observed.size or observed.shape != predicted.shape:
        raise ValueError(
            "observed and predicted must be nonempty and have identical shapes"
        )
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(predicted)):
        raise ValueError("observed and predicted must be finite")
    return observed - predicted


def ssd(observed, predicted):
    return float(np.sum(_residuals(observed, predicted) ** 2))


def rmse(observed, predicted):
    return float(np.sqrt(np.mean(_residuals(observed, predicted) ** 2)))
