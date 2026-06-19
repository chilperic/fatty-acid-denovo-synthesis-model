import numpy as np
from fadns_model.kinetics import mass_action, conservation
from fadns_model.metrics import ssd, rmse
from fadns_model.variants import list_variants


def test_mass_action():
    assert mass_action(2, 3, 4) == 24


def test_conservation():
    assert conservation(10, 1, 2, 3) == 4


def test_metrics_zero():
    x = np.array([1, 2, 3])
    assert ssd(x, x) == 0
    assert rmse(x, x) == 0


def test_variants_present():
    variants = list_variants()
    assert "basic_model" in variants
    assert "coa_inhibition" in variants
