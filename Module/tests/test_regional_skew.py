import pytest
from floodfreq.regional_skew import station_skew_mse, weighted_skew


def test_mse_decreases_with_more_data():
    # more data should always give a more precise (lower MSE) skew estimate
    assert station_skew_mse(skew=0.5, n=100) < station_skew_mse(skew=0.5, n=25)


def test_mse_increases_with_skew_magnitude():
    # a more extreme skew is inherently harder to estimate precisely
    assert station_skew_mse(skew=2.0, n=25) > station_skew_mse(skew=0.5, n=25)
    assert station_skew_mse(skew=1.0, n=25) > station_skew_mse(skew=0.0, n=25)


def test_mse_symmetric_in_skew_sign():
    # the formula uses |skew|, so +G and -G should give the same MSE
    assert station_skew_mse(skew=0.8, n=40) == pytest.approx(station_skew_mse(skew=-0.8, n=40))


def test_weighted_skew_between_station_and_regional():
    # the weighted result must lie between the two inputs (it's a weighted average)
    result = weighted_skew(station_skew=0.8, n=30, regional_skew=0.0, regional_mse=0.302)
    assert 0.0 < result["weighted_skew"] < 0.8


def test_weighted_skew_favors_more_precise_estimate():
    # if the station estimate has a much smaller MSE (e.g. long record), the
    # weighted skew should sit closer to the station value than the regional one
    result = weighted_skew(station_skew=0.8, n=200, regional_skew=0.0, regional_mse=0.302)
    assert result["weighted_skew"] > 0.4  # closer to station (0.8) than regional (0.0)


def test_weighted_skew_equal_inputs_returns_same_value():
    result = weighted_skew(station_skew=0.5, n=30, regional_skew=0.5, regional_mse=0.302)
    assert result["weighted_skew"] == pytest.approx(0.5, abs=1e-9)


def test_review_flag_triggers_on_large_disagreement():
    result = weighted_skew(station_skew=0.9, n=30, regional_skew=0.0, regional_mse=0.302)
    assert result["review_flag"] is not None


def test_review_flag_absent_on_close_agreement():
    result = weighted_skew(station_skew=0.1, n=30, regional_skew=0.0, regional_mse=0.302)
    assert result["review_flag"] is None
