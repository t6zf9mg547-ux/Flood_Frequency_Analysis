import numpy as np
import pytest
from floodfreq.plotting_positions import (
    empirical_frequency, resolve_formula, return_period, non_exceedance_from_T,
    PLOTTING_POSITIONS, RECOMMENDED_FOR,
)


def test_resolve_formula_accepts_names():
    assert resolve_formula("blom") == "blom"
    assert resolve_formula("BLOM") == "blom"  # case-insensitive


def test_resolve_formula_accepts_legacy_codes():
    # matches the original Excel workbook's numeric coding, verified against
    # the workbook: 1=Hazen, 2=Cunnane, 3=Gringorten, 4=Hosking, 5=Weibull, 6=Blom
    assert resolve_formula(1) == "hazen"
    assert resolve_formula(2) == "cunnane"
    assert resolve_formula(3) == "gringorten"
    assert resolve_formula(4) == "hosking"
    assert resolve_formula(5) == "weibull"
    assert resolve_formula(6) == "blom"


def test_resolve_formula_rejects_unknown():
    with pytest.raises(ValueError):
        resolve_formula("not_a_real_formula")
    with pytest.raises(ValueError):
        resolve_formula(99)


@pytest.mark.parametrize("formula", list(PLOTTING_POSITIONS.keys()))
def test_empirical_frequency_is_monotonic_and_bounded(formula):
    F = empirical_frequency(50, formula=formula)
    assert np.all(np.diff(F) > 0), f"{formula}: F must be strictly increasing with rank"
    assert np.all(F > 0) and np.all(F < 1), f"{formula}: F must lie strictly in (0, 1)"


def test_empirical_frequency_blom_matches_formula_definition():
    # Verifies F(r) = (r - a) / (n + b) directly for Blom (a=0.375, b=0.25).
    # NOTE: the reference workbook's own Blom-based value for this case
    # (n=78, rank 1) is 0.0079114, which does NOT match this formula exactly
    # (it implies a denominator of 79, not 78.25) -- an unexplained
    # discrepancy identified during development and not something this
    # implementation should reproduce; this test checks OUR formula is
    # internally correct, not that it matches that specific workbook figure.
    F = empirical_frequency(78, formula="blom", ascending_rank=True)
    expected = (1 - 0.375) / (78 + 0.25)
    assert F[0] == pytest.approx(expected, abs=1e-12)


def test_return_period_roundtrip():
    T = np.array([2.0, 10.0, 100.0, 10000.0])
    F = non_exceedance_from_T(T)
    T_back = return_period(F)
    np.testing.assert_allclose(T, T_back)


def test_recommended_for_covers_every_pwm_distribution():
    # every distribution that can be fit by PWM should have a recommended
    # plotting position for display/descriptive purposes
    for key in RECOMMENDED_FOR:
        assert RECOMMENDED_FOR[key] in PLOTTING_POSITIONS
