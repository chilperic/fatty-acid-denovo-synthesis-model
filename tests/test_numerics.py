import numpy as np
import pytest
from fadns_model.simulation import simulate_time_course


def test_exponential_decay_and_solver_agreement():
    f = lambda t, y, p: -2 * y
    a = simulate_time_course(f, [3.0], {}, t_span=(0, 4), method="LSODA")
    b = simulate_time_course(f, [3.0], {}, t_span=(0, 4), method="Radau")
    np.testing.assert_allclose(a.y[0], 3 * np.exp(-2 * a.t), rtol=2e-7, atol=2e-9)
    np.testing.assert_allclose(a.y, b.y, rtol=3e-7, atol=2e-9)


@pytest.mark.parametrize(
    "kw",
    [
        {"y0": []},
        {"y0": [np.nan]},
        {"n_points": 1},
        {"rtol": 0},
        {"atol": [-1]},
        {"t_span": (2, 1)},
        {"t_eval": [0, 0, 1]},
    ],
)
def test_invalid_solver_inputs(kw):
    args = dict(y0=[1.0])
    args.update(kw)
    with pytest.raises(ValueError):
        simulate_time_course(lambda t, y, p: -y, params={}, **args)


def test_rhs_shape_and_nonfinite_are_errors():
    for f in [lambda t, y, p: [1.0, 2.0], lambda t, y, p: [np.nan]]:
        with pytest.raises(ValueError):
            simulate_time_course(f, [1.0], {})
