"""ODE simulation utilities."""

from __future__ import annotations

from collections.abc import Callable, Mapping

import numpy as np
from scipy.integrate import solve_ivp


def simulate_time_course(rhs: Callable, y0, params: Mapping[str, float] | None = None, t_span=(0.0, 10.0), n_points: int = 200, method: str = "LSODA"):
    params = {} if params is None else params
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    return solve_ivp(lambda t, y: rhs(t, y, params), t_span=t_span, y0=np.asarray(y0, dtype=float), t_eval=t_eval, method=method)
