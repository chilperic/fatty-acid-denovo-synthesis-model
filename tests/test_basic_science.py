import numpy as np
import pytest
from fadns_model import basic
from fadns_model.legacy_variants import basic_model
from fadns_model.metrics import ssd
from fadns_model.simulation import simulate_time_course


def test_stoichiometric_refactor_matches_original_equations():
    rng = np.random.default_rng(42)
    for _ in range(20):
        y = rng.uniform(0.01, 3, 21)
        p = dict(zip(basic.PARAMETERS, rng.uniform(0.01, 2, 11)))
        np.testing.assert_allclose(
            basic.rhs(0, y, p), basic_model.system(0, y, p), rtol=1e-14, atol=1e-14
        )


def test_enzyme_and_acetyl_unit_conservation():
    enzyme = np.zeros(21)
    enzyme[10:20] = 1
    carbon = np.array([1, 1, *range(2, 10), 0, *range(1, 10), 0])
    np.testing.assert_array_equal(enzyme @ basic.STOICHIOMETRY, 0)
    np.testing.assert_array_equal(carbon @ basic.STOICHIOMETRY, 0)
    y0 = np.zeros(21)
    y0[[0, 1, 10, 20]] = [2, 10, 0.1, 20]
    p = dict.fromkeys(basic.PARAMETERS, 0.1)
    out = simulate_time_course(basic.rhs, y0, p, t_span=(0, 10))
    np.testing.assert_allclose(enzyme @ out.y, enzyme @ y0, atol=1e-10)
    np.testing.assert_allclose(carbon @ out.y, carbon @ y0, atol=1e-9)
    assert np.min(out.y) >= -1e-9


def test_metrics_reject_broadcasting_and_missing_observations():
    with pytest.raises(ValueError):
        ssd([1, 2], [[1, 2]])
    with pytest.raises(ValueError):
        ssd([np.nan], [1])


def test_headerless_observation_retains_first_row(tmp_path):
    from fadns_model.data_io import load_fatty_acid_csv

    path = tmp_path / "observed.csv"
    path.write_text("0,-0.001\n20,0.038\n")
    frame = load_fatty_acid_csv(path)
    assert len(frame) == 2
    assert frame.iloc[0, 0] == 0
    assert frame.iloc[0, 1] == -0.001
