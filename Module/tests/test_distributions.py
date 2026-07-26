"""
Regression tests for distribution fitting. Expected values here were
cross-checked by hand against the original Analyse_Frequentielle_V03.xls
workbook's own computed parameters during development (see conversation
history / summary.txt outputs) -- if these start failing, either a real
bug was introduced, or a fitting-method change was intentional and these
values need to be updated (and re-validated against the workbook again,
not just accepted because the test now passes with new numbers).
"""
import numpy as np
import pytest
from scipy import stats as ss
from floodfreq.distributions import fit, DISTRIBUTIONS, available_distributions


def test_available_distributions_has_nine():
    assert len(available_distributions()) == 9


# -- Parameters cross-checked against the workbook -- #

def test_normal_mle(reference_data):
    x, _ = reference_data
    r = fit("normal", x, method="mle")
    # workbook mean matches exactly; std differs by ddof convention (workbook
    # uses n-1, scipy MLE uses n) -- both are standard, just different
    assert r.params[0] == pytest.approx(188.6397435897436, abs=1e-6)
    assert r.params[1] == pytest.approx(47.356521062191014, abs=1e-6)


def test_lognormal2_mle(reference_data):
    x, _ = reference_data
    r = fit("lognormal2", x, method="mle")
    assert r.params[0] == pytest.approx(5.209076598072831, abs=1e-6)
    assert r.params[1] == pytest.approx(0.2490019956317524, abs=1e-6)


def test_gumbel_pwm_matches_workbook_closely(reference_data):
    x, _ = reference_data
    r = fit("gumbel", x, method="pwm")
    # workbook: a=37.4567, b=167.0191 -- matches to 4 decimal places
    assert r.params[0] == pytest.approx(167.0191, abs=1e-3)
    assert r.params[1] == pytest.approx(37.4567, abs=1e-3)


def test_exponential_pwm_matches_workbook_closely(reference_data):
    x, _ = reference_data
    r = fit("exponential", x, method="pwm")
    # workbook: e=136.7137, a=51.9260
    assert r.params[0] == pytest.approx(136.7137, abs=1e-3)
    assert r.params[1] == pytest.approx(51.9260, abs=1e-3)


def test_gev_pwm_close_to_workbook(reference_data):
    x, _ = reference_data
    r = fit("gev", x, method="pwm")
    # workbook: u=168.1977, alpha=39.6929, k=0.0667 -- small differences vs.
    # our PWM implementation are expected here (different GEV shape
    # approximation order), so tolerance is looser than Gumbel/Exponential
    k, u, alpha = r.params
    assert u == pytest.approx(168.1977, abs=0.5)
    assert alpha == pytest.approx(39.6929, abs=0.5)
    assert k == pytest.approx(0.0667, abs=0.01)


# -- Sanity/invariant checks (not tied to one dataset) -- #

def test_mom_unavailable_raises_for_gev(reference_data):
    x, _ = reference_data
    with pytest.raises(ValueError):
        fit("gev", x, method="mom")


def test_unknown_method_raises(reference_data):
    x, _ = reference_data
    with pytest.raises(ValueError):
        fit("normal", x, method="not_a_method")


def test_quantile_increases_with_return_period(reference_data):
    x, _ = reference_data
    r = fit("gev", x, method="pwm")
    q = r.quantile([2, 10, 100, 1000, 10000])
    assert np.all(np.diff(q) > 0), "design flood must increase monotonically with T"


def test_cdf_ppf_are_inverses(reference_data):
    x, _ = reference_data
    r = fit("lognormal2", x, method="mle")
    F = np.array([0.1, 0.5, 0.9, 0.99])
    q = r.ppf(F)
    F_back = r.cdf(q)
    np.testing.assert_allclose(F, F_back, atol=1e-8)


def test_anderson_darling_matches_scipy_for_normal():
    # cross-checked against scipy.stats.anderson's independent implementation
    # (see conversation history) -- confirms the general A^2 formula is
    # correctly implemented, not just plausible-looking
    rng = np.random.default_rng(0)
    x = rng.normal(100, 15, 60)
    r = fit("normal", x, method="mle")  # MLE => ddof=0, matches scipy's convention here
    our_A2 = r.anderson_darling_statistic(x)

    mean, std = x.mean(), x.std(ddof=0)
    F = np.clip(ss.norm.cdf(np.sort(x), loc=mean, scale=std), 1e-12, 1 - 1e-12)
    n = x.size
    i = np.arange(1, n + 1)
    expected_A2 = -n - np.sum((2 * i - 1) * (np.log(F) + np.log(1 - F[::-1]))) / n

    assert our_A2 == pytest.approx(expected_A2, rel=1e-9)


def test_loglik_is_finite_for_all_default_fits(reference_data):
    x, _ = reference_data
    # Exponential is deliberately excluded: for this dataset, its PWM-fitted
    # location parameter exceeds the minimum observed value, correctly
    # giving -inf log-likelihood (some data points fall outside the fitted
    # distribution's support). That's a real, meaningful result -- see
    # summary.txt / conversation history -- not something to paper over.
    for key in DISTRIBUTIONS:
        if key == "exponential":
            continue
        method = "mom" if key == "gamma2" else ("mle" if key in ("normal", "lognormal2") else "pwm")
        r = fit(key, x, method=method)
        ll = r.loglik(x)
        assert np.isfinite(ll), f"{key} ({method}): log-likelihood should be finite for this well-behaved dataset"


def test_exponential_pwm_can_be_invalid_for_this_dataset(reference_data):
    """
    Documents a known, expected result rather than treating it as a bug:
    the 2-parameter Exponential's PWM-fitted location can exceed the
    minimum observed value, putting some data points outside its support.
    """
    x, _ = reference_data
    r = fit("exponential", x, method="pwm")
    loc = r.params[0]
    assert loc > x.min(), "this assertion documents why loglik is -inf for this case"
    assert not np.isfinite(r.loglik(x))
