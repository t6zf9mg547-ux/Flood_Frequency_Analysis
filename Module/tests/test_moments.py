import numpy as np
import pytest
from floodfreq.moments import ordinary_moments, probability_weighted_moments, pwm_to_lmoments, summarize


def test_ordinary_moments_normal_data(reference_data):
    x, _ = reference_data
    m = ordinary_moments(x)
    # cross-checked against the reference workbook's Stat sheet
    assert m["mean"] == pytest.approx(188.63974358974357, abs=1e-6)
    assert m["n"] == 78


def test_pwm_b0_equals_mean(reference_data):
    # b0 (the 0th PWM) is, by definition, just the sample mean
    x, _ = reference_data
    b = probability_weighted_moments(x, formula="weibull")
    assert b["b0"] == pytest.approx(x.mean(), abs=1e-9)


def test_lmoments_l1_equals_mean(reference_data):
    # l1 (the 1st L-moment) is also, by definition, the sample mean
    x, _ = reference_data
    b = probability_weighted_moments(x, formula="weibull")
    lm = pwm_to_lmoments(b)
    assert lm["l1"] == pytest.approx(x.mean(), abs=1e-9)


def test_lmoments_l2_is_positive(reference_data):
    # l2 (L-scale) must be positive for any non-degenerate sample
    x, _ = reference_data
    b = probability_weighted_moments(x, formula="weibull")
    lm = pwm_to_lmoments(b)
    assert lm["l2"] > 0


def test_summarize_bundles_everything(reference_data):
    x, _ = reference_data
    s = summarize(x)
    for key in ("mean", "std", "CV", "CS", "CK", "b0", "b1", "b2", "b3", "l1", "l2", "l3", "l4", "t3", "t4"):
        assert key in s

    # cross-checked against the reference workbook's PWM section (b0-b3)
    assert s["b0"] == pytest.approx(188.639744, abs=1e-3)


def test_symmetric_data_has_near_zero_skew():
    rng = np.random.default_rng(0)
    x = rng.normal(100, 10, 5000)  # large symmetric sample
    m = ordinary_moments(x)
    assert abs(m["CS"]) < 0.15  # should be close to 0 for a large symmetric sample


def test_degenerate_series_raises_or_nans():
    # all-identical values: std = 0, CV/CS/CK become undefined (nan or inf)
    x = np.full(20, 100.0)
    m = ordinary_moments(x)
    assert m["std"] == 0
