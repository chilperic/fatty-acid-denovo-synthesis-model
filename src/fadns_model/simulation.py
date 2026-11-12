"""Validated integration boundary; equations are supplied by the caller."""

from collections.abc import Callable, Mapping

import numpy as np
from scipy.integrate import solve_ivp


def simulate_time_course(
    rhs: Callable,
    y0,
    params: Mapping | None = None,
    t_span=(0.0, 10.0),
    n_points=200,
    method="LSODA",
    rtol=1e-8,
    atol=1e-10,
    *,
    t_eval=None,
    max_step=np.inf,
):
    """Integrate an ODE, rejecting invalid inputs and incomplete/nonfinite results.

    ``atol`` may be a positive scalar or a vector with one entry per state.
    No clipping or rescaling of scientific states is performed.
    """
    y0 = np.asarray(y0, dtype=float)
    span = np.asarray(t_span, dtype=float)
    if y0.ndim != 1 or not y0.size or not np.all(np.isfinite(y0)):
        raise ValueError("y0 must be a nonempty finite vector")
    if span.shape != (2,) or not np.all(np.isfinite(span)) or span[1] <= span[0]:
        raise ValueError("t_span must contain two finite increasing times")
    tolerance = np.asarray(atol, dtype=float)
    if (
        tolerance.shape not in ((), y0.shape)
        or np.any(~np.isfinite(tolerance))
        or np.any(tolerance <= 0)
    ):
        raise ValueError("atol must be positive and scalar or match y0")
    if not np.isscalar(rtol) or not np.isfinite(rtol) or rtol <= 0:
        raise ValueError("rtol must be a positive finite scalar")
    if not np.isscalar(max_step) or np.isnan(max_step) or max_step <= 0:
        raise ValueError("max_step must be positive")
    if t_eval is None:
        if (
            isinstance(n_points, bool)
            or not isinstance(n_points, (int, np.integer))
            or n_points < 2
        ):
            raise ValueError("n_points must be an integer >= 2")
        times = np.linspace(*span, n_points)
    else:
        times = np.asarray(t_eval, dtype=float)
        if (
            times.ndim != 1
            or not times.size
            or not np.all(np.isfinite(times))
            or np.any(np.diff(times) <= 0)
            or times[0] < span[0]
            or times[-1] > span[1]
        ):
            raise ValueError("t_eval must be finite, increasing, and within t_span")
    parameters = dict(params or {})

    def checked_rhs(t, y):
        value = np.asarray(rhs(t, y, parameters), dtype=float)
        if value.shape != y0.shape or not np.all(np.isfinite(value)):
            raise ValueError("rhs must return a finite vector matching y0")
        return value

    checked_rhs(span[0], y0)
    result = solve_ivp(
        checked_rhs,
        tuple(span),
        y0,
        t_eval=times,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    if (
        not result.success
        or result.y.shape != (y0.size, times.size)
        or not np.all(np.isfinite(result.y))
    ):
        raise RuntimeError(f"Integration failed or incomplete: {result.message}")
    result.solver_config = {
        "method": method,
        "rtol": float(rtol),
        "atol": tolerance.tolist(),
        "max_step": float(max_step),
    }
    return result
