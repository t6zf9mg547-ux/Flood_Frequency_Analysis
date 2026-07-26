import numpy as np
import pytest
from floodfreq.analysis import FloodFrequencyAnalysis
from floodfreq.distributions import DISTRIBUTIONS


def test_fit_all_fits_every_distribution(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    assert len(ffa.fits) == len(DISTRIBUTIONS)


def test_goodness_of_fit_table_sorted_by_aic(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    table = ffa.goodness_of_fit_table()
    assert list(table["AIC"]) == sorted(table["AIC"])


def test_best_fit_matches_table_top_row(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    best_key, best_method = ffa.best_fit()
    table = ffa.goodness_of_fit_table()
    assert best_key == table.iloc[0]["key"]
    assert best_method == table.iloc[0]["method"]


def test_akaike_weights_sum_to_one(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    weights = ffa.akaike_weights()
    assert weights["akaike_weight"].sum() == pytest.approx(1.0, abs=1e-9)
    assert (weights["akaike_weight"] >= 0).all()


def test_akaike_weights_best_model_has_zero_delta(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    weights = ffa.akaike_weights()
    assert weights.iloc[0]["delta"] == pytest.approx(0.0, abs=1e-9)
    assert weights.iloc[0]["akaike_weight"] == weights["akaike_weight"].max()


def test_model_averaged_quantile_increases_with_return_period(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    ma = ffa.model_averaged_quantile_table(return_periods=[2, 10, 100, 1000, 10000])
    assert np.all(np.diff(ma["Q_model_averaged"]) > 0)


def test_model_averaged_spread_grows_with_extrapolation(reference_data):
    # candidate distributions should agree more near the data and disagree
    # more when extrapolated far beyond the record -- this is the core
    # justification for reporting between_model_sd at all
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    ma = ffa.model_averaged_quantile_table(return_periods=[10, 10000])
    sd_near = ma[ma["T"] == 10]["between_model_sd"].iloc[0]
    sd_far = ma[ma["T"] == 10000]["between_model_sd"].iloc[0]
    assert sd_far > sd_near


def test_bootstrap_ci_bounds_are_ordered(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    best_key, best_method = ffa.best_fit()
    ci = ffa.bootstrap_ci(best_key, best_method, [10, 100], n_boot=200, random_state=0)
    assert (ci["lower"] <= ci["median"]).all()
    assert (ci["median"] <= ci["upper"]).all()


def test_data_quality_is_cached(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    first = ffa.data_quality()
    second = ffa.data_quality()
    assert first is second  # same object -> confirms caching, not recomputation


def test_regional_skew_adds_two_extra_candidates(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    n_before = len(DISTRIBUTIONS)
    ffa.fit_all(regional_skew=0.0, regional_mse=0.15)
    assert len(ffa.fits) == n_before + 2
    assert ("pearson3", "mom_weighted_skew") in ffa.fits
    assert ("logpearson3", "mom_weighted_skew") in ffa.fits


def test_generate_recommendation_runs_without_error(reference_data):
    # smoke test: the full summary-generation path should run end to end
    # without raising, across every section (data quality, regional skew,
    # recommendation, model averaging, extrapolation warning)
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all(regional_skew=0.0, regional_mse=0.15)
    text = ffa.generate_recommendation(ci_n_boot=50)  # low n_boot for test speed
    assert "RECOMMENDED DISTRIBUTION" in text
    assert "DATA QUALITY CHECKS" in text
    assert "MODEL-AVERAGED DESIGN FLOOD" in text
    assert "REGIONAL SKEW" in text
    assert "EXTRAPOLATION WARNING" in text  # default return periods go to 10000


# -- Hard input validation: clear errors, not confusing tracebacks -- #

def test_rejects_empty_data():
    with pytest.raises(ValueError, match="empty"):
        FloodFrequencyAnalysis(np.array([]))


def test_rejects_nan():
    x = np.array([100.0, np.nan, 105.0, 110.0, 95.0, 102.0])
    with pytest.raises(ValueError, match="[Mm]issing"):
        FloodFrequencyAnalysis(x)


def test_rejects_infinite_values():
    x = np.array([100.0, np.inf, 105.0, 110.0, 95.0, 102.0])
    with pytest.raises(ValueError, match="[Ii]nfinite"):
        FloodFrequencyAnalysis(x)


def test_rejects_non_positive_values():
    x = np.array([100.0, -5.0, 105.0, 110.0, 95.0, 102.0])
    with pytest.raises(ValueError, match="[Nn]egative"):
        FloodFrequencyAnalysis(x)

    x_zero = np.array([100.0, 0.0, 105.0, 110.0, 95.0, 102.0])
    with pytest.raises(ValueError):
        FloodFrequencyAnalysis(x_zero)


def test_rejects_too_few_points():
    x = np.array([100.0, 110.0, 105.0])  # only 3 points
    with pytest.raises(ValueError, match="5"):
        FloodFrequencyAnalysis(x)


def test_rejects_mismatched_years_length():
    x = np.array([100.0, 110.0, 105.0, 108.0, 112.0, 98.0])
    years = np.array([2000, 2001, 2002])  # wrong length
    with pytest.raises(ValueError, match="years"):
        FloodFrequencyAnalysis(x, years=years)


def test_accepts_minimum_valid_size():
    # exactly 5 points (the minimum) should NOT raise
    x = np.array([100.0, 110.0, 105.0, 108.0, 95.0])
    ffa = FloodFrequencyAnalysis(x)  # should not raise
    assert ffa.n == 5


def test_valid_reference_data_does_not_raise(reference_data):
    x, years = reference_data
    FloodFrequencyAnalysis(x, station_id="test", years=years)  # should not raise
