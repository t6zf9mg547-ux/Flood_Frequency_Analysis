"""
Tests for stationarity/trend, outlier, and input-validation checks.

The Mann-Kendall direction tests here specifically guard against the sign
bug found and fixed during development (an early version reported
"decreasing" for a clearly increasing series and vice versa) -- these
tests exist so that particular bug class can never silently return.
"""
import numpy as np
import pytest
from floodfreq.data_quality import (
    mann_kendall_test, sens_slope, grubbs_outlier_test, validate_series, run_all,
)


# -- Mann-Kendall: direction must be correct, not just "some" statistic -- #

def test_mann_kendall_detects_increasing_trend():
    rng = np.random.default_rng(0)
    x = np.arange(50) * 2.0 + rng.normal(0, 5, 50)  # strong, unambiguous upward trend
    result = mann_kendall_test(x)
    assert result["S"] > 0, "S must be positive for an increasing series"
    assert result["Z"] > 0
    assert result["significant"] is True
    assert "increasing" in result["trend"]


def test_mann_kendall_detects_decreasing_trend():
    rng = np.random.default_rng(0)
    x = -np.arange(50) * 2.0 + rng.normal(0, 5, 50)  # strong, unambiguous downward trend
    result = mann_kendall_test(x)
    assert result["S"] < 0, "S must be negative for a decreasing series"
    assert result["Z"] < 0
    assert result["significant"] is True
    assert "decreasing" in result["trend"]


def test_mann_kendall_no_trend_on_pure_noise():
    rng = np.random.default_rng(1)
    x = rng.normal(100, 10, 50)
    result = mann_kendall_test(x)
    assert result["significant"] is False
    assert result["trend"] == "no significant trend"


def test_mann_kendall_on_reference_data_is_not_significant(reference_data):
    # this specific station's record shows no significant trend -- if this
    # ever flips to significant, something upstream changed unexpectedly
    x, _ = reference_data
    result = mann_kendall_test(x)
    assert result["significant"] is False


# -- Sen's slope: sign and rough magnitude should match a known trend -- #

def test_sens_slope_recovers_known_slope():
    rng = np.random.default_rng(0)
    years = np.arange(1950, 2000)
    true_slope = 3.5
    x = true_slope * (years - 1950) + 100 + rng.normal(0, 5, len(years))
    result = sens_slope(x, t=years)
    assert result["slope"] == pytest.approx(true_slope, abs=0.5)


def test_sens_slope_sign_matches_direction():
    x_up = np.arange(30) + np.random.default_rng(0).normal(0, 0.1, 30)
    x_down = -np.arange(30) + np.random.default_rng(0).normal(0, 0.1, 30)
    assert sens_slope(x_up)["slope"] > 0
    assert sens_slope(x_down)["slope"] < 0


# -- Grubbs' outlier test -- #

def test_grubbs_detects_injected_high_outlier():
    rng = np.random.default_rng(2)
    x = rng.normal(100, 10, 30)
    x[5] = 500  # obvious outlier
    result = grubbs_outlier_test(x, log_space=False)
    assert result["high_outlier_flagged"] is True
    assert result["high_outlier_value"] == 500


def test_grubbs_no_false_positive_on_clean_data():
    rng = np.random.default_rng(3)
    x = rng.normal(100, 10, 50)
    result = grubbs_outlier_test(x, log_space=False)
    assert result["high_outlier_flagged"] is False
    assert result["low_outlier_flagged"] is False


def test_grubbs_on_reference_data_flags_nothing(reference_data):
    x, _ = reference_data
    result = grubbs_outlier_test(x)
    assert result["high_outlier_flagged"] is False
    assert result["low_outlier_flagged"] is False


# -- Input validation -- #

def test_validate_series_flags_short_record():
    x = np.array([100.0, 110.0, 105.0])
    warnings = validate_series(x)
    assert any("short record" in w.lower() or "record" in w.lower() for w in warnings)


def test_validate_series_flags_nan():
    x = np.array([100.0, np.nan, 105.0, 110.0, 95.0] * 5)
    warnings = validate_series(x)
    assert any("missing" in w.lower() for w in warnings)


def test_validate_series_flags_negative_values():
    x = np.array([100.0, -5.0, 105.0, 110.0, 95.0] * 5)
    warnings = validate_series(x)
    assert any("negative" in w.lower() for w in warnings)


def test_validate_series_flags_duplicate_years():
    x = np.arange(20, dtype=float) + 100
    years = np.array([2000] * 20) + np.arange(20)
    years[5] = years[4]  # duplicate
    warnings = validate_series(x, years=years)
    assert any("duplicate" in w.lower() for w in warnings)


def test_validate_series_clean_data_has_no_warnings(reference_data):
    x, years = reference_data
    warnings = validate_series(x, years=years)
    # this dataset has 2 known missing years (1949, 1951) -- that's the one
    # expected warning, nothing else
    assert len(warnings) == 1
    assert "missing" in warnings[0].lower()


def test_run_all_bundles_every_check(reference_data):
    x, years = reference_data
    result = run_all(x, years=years)
    assert set(result.keys()) == {"validation_warnings", "mann_kendall", "sens_slope", "grubbs"}
