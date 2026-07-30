"""
Tests for the CIFAM climate-adjustment module (floodfreq.climate_adjustment).

Coverage, in the spirit of the rest of this project's tests:
  (a) hand-computable toy cases that pin the closed-form algebra for each of
      the three distributions with a published closed form;
  (b) an integration / sanity test reproducing the paper's Lom Pangar worked
      example (Grijsen & Lino 2026, Fig. 4). The baseline Gumbel parameters
      are NOT given numerically in the paper -- they were digitized from
      Fig. 4 (see the module and the session notes), so this test asserts the
      paper's own *rounded, qualitative* headline claims (percent increases,
      PMF agreement, design range) with deliberately loose tolerances, not
      exact figures;
  (c) cross-validation that the Monte Carlo path agrees with the closed form
      for all three distributions within simulation noise.

The math in the module was re-derived from first principles and cross-checked
against the PDF's Table 1 / eqs. 5-8 images (not merely transcribed); these
tests guard that derivation.
"""
import numpy as np
import pytest

from floodfreq.climate_adjustment import (
    climate_adjusted_quantiles,
    mc_climate_adjusted_quantiles,
    reduced_variate,
    ClimateInputs,
    EULER,
    GUMBEL_SIGMA_OVER_BETA,
    CLOSED_FORM_DISTRIBUTIONS,
    _gumbel_closed_form,
    _lognormal_closed_form,
    _pearson3_closed_form,
)


# --------------------------------------------------------------------------- #
# Basic contract                                                              #
# --------------------------------------------------------------------------- #
def test_reduced_variate_known_values():
    # t_p = -ln(-ln(1 - 1/T)) for standard return periods
    assert reduced_variate(10000) == pytest.approx(9.2102, abs=1e-3)
    assert reduced_variate(100000) == pytest.approx(11.5129, abs=1e-3)
    # T = 1/(1-exp(-1)) gives t_p = 0
    T0 = 1.0 / (1.0 - np.exp(-1.0))
    assert reduced_variate(T0) == pytest.approx(0.0, abs=1e-9)


def test_climate_inputs_rejects_negative_sd():
    with pytest.raises(ValueError):
        ClimateInputs(0.3, -0.1, 0.1, 0.18)
    with pytest.raises(ValueError):
        ClimateInputs(0.3, 0.1, 0.1, -0.18)


def test_unknown_distribution_closed_form_raises():
    data = np.array([100.0, 120, 90, 110, 130, 95, 105])
    with pytest.raises(ValueError):
        climate_adjusted_quantiles(data, "gev", [100], 0.3, 0.18, 0.1, 0.18)


def test_zero_climate_change_reduces_to_baseline():
    # With delta1=delta2=tau1=tau2=0 the climate-adjusted central estimate and
    # CI must collapse exactly onto the baseline central estimate and CI.
    rng = np.random.default_rng(0)
    data = rng.gumbel(900, 500, size=50)
    for dist in CLOSED_FORM_DISTRIBUTIONS:
        r = climate_adjusted_quantiles(data, dist, [10, 100, 1000],
                                       0.0, 0.0, 0.0, 0.0)
        np.testing.assert_allclose(r.climate_point, r.baseline_point, rtol=1e-9)
        np.testing.assert_allclose(r.climate_lower, r.baseline_lower, rtol=1e-9)
        np.testing.assert_allclose(r.climate_upper, r.baseline_upper, rtol=1e-9)


# --------------------------------------------------------------------------- #
# (a) Hand-computable toy cases -- GUMBEL                                      #
# --------------------------------------------------------------------------- #
def test_gumbel_closed_form_matches_hand_computation():
    # Toy baseline chosen for round numbers.
    mu0, sigma0, N = 1000.0, 500.0, 50
    d1, d2, t1, t2 = 0.30, 0.18, 0.10, 0.18
    T = 1000.0
    tp = reduced_variate(T)
    c = EULER
    beta0 = sigma0 / GUMBEL_SIGMA_OVER_BETA

    # Baseline point (eq. 4): xp0 = mu0 + beta0*(tp - c)
    xp0 = mu0 + beta0 * (tp - c)
    # Climate point (deterministic part of eq. 7)
    xp_cc = xp0 + mu0 * d1 + beta0 * t1 * (tp - c)
    # Variance (eq. 8)
    var_s = (beta0 ** 2 / N) * (1.1086 + 0.514 * tp + 0.6079 * tp ** 2)
    var_c = d2 ** 2 * mu0 ** 2 + beta0 ** 2 * t2 ** 2 * (tp - c) ** 2
    from scipy.stats import norm
    z = norm.ppf(0.975)

    out = _gumbel_closed_form(mu0, sigma0, N, np.array([T]),
                              ClimateInputs(d1, d2, t1, t2), 95.0)
    assert out["baseline_point"][0] == pytest.approx(xp0, rel=1e-12)
    assert out["climate_point"][0] == pytest.approx(xp_cc, rel=1e-12)
    assert out["baseline_upper"][0] == pytest.approx(xp0 + z * np.sqrt(var_s), rel=1e-12)
    assert out["climate_upper"][0] == pytest.approx(xp_cc + z * np.sqrt(var_s + var_c), rel=1e-12)


def test_gumbel_climate_variance_uses_euler_c_not_a_separate_c1():
    # Guards the eq.(8) "(tp - c1)^2" ambiguity: the last term must use Euler's
    # c (0.57722). If someone "fixed" it to c1=0 (or 1), the climate CI width at
    # a given tp would change; assert against the c-based value explicitly.
    mu0, sigma0, N = 1000.0, 500.0, 50
    d2, t2 = 0.18, 0.18
    T = 100000.0
    tp = reduced_variate(T)
    beta0 = sigma0 / GUMBEL_SIGMA_OVER_BETA
    out = _gumbel_closed_form(mu0, sigma0, N, np.array([T]),
                              ClimateInputs(0.0, d2, 0.0, t2), 95.0)
    from scipy.stats import norm
    z = norm.ppf(0.975)
    var_s = (beta0 ** 2 / N) * (1.1086 + 0.514 * tp + 0.6079 * tp ** 2)
    var_c_correct = d2 ** 2 * mu0 ** 2 + beta0 ** 2 * t2 ** 2 * (tp - EULER) ** 2
    half_width = out["climate_upper"][0] - out["climate_point"][0]
    assert half_width == pytest.approx(z * np.sqrt(var_s + var_c_correct), rel=1e-12)
    # And it must NOT equal the c1=0 variant.
    var_c_wrong = d2 ** 2 * mu0 ** 2 + beta0 ** 2 * t2 ** 2 * tp ** 2
    assert half_width != pytest.approx(z * np.sqrt(var_s + var_c_wrong), rel=1e-6)


# --------------------------------------------------------------------------- #
# (a) Hand-computable toy cases -- LOG-NORMAL                                  #
# --------------------------------------------------------------------------- #
def test_lognormal_closed_form_matches_hand_computation():
    mu0, sigma0, N = 1000.0, 300.0, 60
    d1, d2, t1, t2 = 0.20, 0.10, 0.05, 0.12
    T = 500.0
    tp = reduced_variate(T, kind="normal")

    cv2 = (sigma0 / mu0) ** 2
    sigy2 = np.log(1 + cv2)
    sigy = np.sqrt(sigy2)
    muy = np.log(mu0) - 0.5 * sigy2
    xp0 = np.exp(muy + sigy * tp)

    var_yp_s = (1 + 0.5 * tp ** 2) * sigy2 / N
    # central estimate via EXACT moment mapping of the shifted physical moments
    mu_cc = mu0 * (1 + d1)
    sig_cc = sigma0 * (1 + t1)
    sigy2_cc = np.log(1 + (sig_cc / mu_cc) ** 2)
    sigy_cc = np.sqrt(sigy2_cc)
    muy_cc = np.log(mu_cc) - 0.5 * sigy2_cc
    yp_cc = muy_cc + sigy_cc * tp
    # exact delta-method climate variance (see module): g1,g2 sensitivities
    r_cc = sig_cc / mu_cc
    rho = r_cc ** 2 / (1 + r_cc ** 2)
    g1 = (1 + rho) - rho * tp / sigy_cc
    g2 = rho * (tp / sigy_cc - 1)
    var_yp_c = g1 ** 2 * d2 ** 2 + g2 ** 2 * t2 ** 2
    from scipy.stats import norm
    z = norm.ppf(0.975)

    out = _lognormal_closed_form(mu0, sigma0, N, np.array([T]),
                                 ClimateInputs(d1, d2, t1, t2), 95.0)
    assert out["baseline_point"][0] == pytest.approx(xp0, rel=1e-12)
    assert out["climate_point"][0] == pytest.approx(np.exp(yp_cc), rel=1e-12)
    assert out["baseline_upper"][0] == pytest.approx(np.exp((muy + sigy * tp) + z * np.sqrt(var_yp_s)), rel=1e-12)
    assert out["climate_upper"][0] == pytest.approx(np.exp(yp_cc + z * np.sqrt(var_yp_s + var_yp_c)), rel=1e-12)


# --------------------------------------------------------------------------- #
# (a) Hand-computable toy cases -- PEARSON III                                 #
# --------------------------------------------------------------------------- #
def test_pearson3_partial_derivative_matches_finite_difference():
    # dxp/dgamma from the module vs. a central finite difference of xp=beta*zp.
    from floodfreq.climate_adjustment import _pe3_zp, _pe3_dxp_dgamma
    for gamma in (0.5, 1.0, 2.0, 5.0, 20.0, 100.0):
        for beta in (10.0, 300.0):
            for tp in (-1.0, 0.0, 1.0, 3.0, 5.0, 9.0):
                h = 1e-6 * gamma
                num = (beta * _pe3_zp(gamma + h, tp) - beta * _pe3_zp(gamma - h, tp)) / (2 * h)
                ana = _pe3_dxp_dgamma(beta, gamma, tp)
                assert ana == pytest.approx(num, rel=1e-4, abs=1e-4)


def test_pearson3_parameter_shifts_follow_moment_equations():
    # The climate-shifted shape/scale must equal the EXACT moment mapping
    #   gamma_cc = gamma0*(1+delta1)^2/(1+tau1)^2
    #   beta_cc  = beta0 *(1+tau1)^2/(1+delta1)
    # (the paper's d(ln gamma)=2A-2B, d(ln beta)=2B-A are the first-order form).
    mu0, sigma0, N = 1000.0, 300.0, 60
    d1, t1 = 0.20, 0.05
    gamma0 = (mu0 / sigma0) ** 2
    beta0 = sigma0 ** 2 / mu0
    gamma_expected = gamma0 * (1 + d1) ** 2 / (1 + t1) ** 2
    beta_expected = beta0 * (1 + t1) ** 2 / (1 + d1)

    from floodfreq.climate_adjustment import _pe3_zp
    T = np.array([100.0, 1000.0])
    tp = reduced_variate(T, kind="normal")
    out = _pearson3_closed_form(mu0, sigma0, N, T, ClimateInputs(d1, 0, t1, 0), 95.0)
    # climate_point = beta_cc * zp(gamma_cc, tp); recover gamma_cc from the
    # ratio at two tp (monotone in gamma), then beta_cc.
    from scipy.optimize import brentq
    r = out["climate_point"][1] / out["climate_point"][0]

    def f(g):
        return _pe3_zp(g, tp[1]) / _pe3_zp(g, tp[0]) - r
    # bracket around the expected value
    lo, hi = gamma_expected * 0.1, gamma_expected * 10
    g_rec = brentq(f, lo, hi)
    b_rec = out["climate_point"][0] / _pe3_zp(g_rec, tp[0])
    assert g_rec == pytest.approx(gamma_expected, rel=1e-4)
    assert b_rec == pytest.approx(beta_expected, rel=1e-4)


# --------------------------------------------------------------------------- #
# (b) Lom Pangar integration / sanity test (approximate -- digitized baseline) #
# --------------------------------------------------------------------------- #
# Baseline Gumbel parameters digitized from Grijsen & Lino (2026) Fig. 4:
#   the red "Daily max. with CC" central line reads xp_cc ~ 1566 + 181*tp;
#   inverting the CIFAM shift (delta1=0.30, tau1=0.10) gives u0~1190, beta0~165
#   => mu0~1285, sigma0~211. N~80 best-fits the digitized 97.5%-with-CC line.
# These are read off a figure and are therefore approximate; the assertions
# below use loose tolerances matching the paper's own qualitative wording.
LP_MU0, LP_SIGMA0, LP_N = 1285.0, 211.0, 80
LP_CLIMATE = dict(delta1=0.30, delta2=0.18, tau1=0.10, tau2=0.18)
LP_PMF = 4140.0


def _lp_data():
    """Synthetic Gumbel sample whose moments hit the digitized (mu0, sigma0)
    at exactly N=LP_N, so climate_adjusted_quantiles (which recomputes moments
    from data) sees the intended baseline."""
    beta0 = LP_SIGMA0 / GUMBEL_SIGMA_OVER_BETA
    u0 = LP_MU0 - EULER * beta0
    # deterministic quantiles of a Gumbel at the plotting positions, then
    # rescale to force mean/std to (mu0, sigma0) exactly.
    p = (np.arange(1, LP_N + 1) - 0.44) / (LP_N + 0.12)  # Gringorten
    x = u0 - beta0 * np.log(-np.log(p))
    x = (x - x.mean()) / x.std() * LP_SIGMA0 + LP_MU0
    return x


def test_lompangar_reproduces_paper_headline_claims():
    data = _lp_data()
    m = data.mean(); s = data.std()
    assert m == pytest.approx(LP_MU0, rel=1e-6)
    assert s == pytest.approx(LP_SIGMA0, rel=1e-6)

    T = np.array([50, 100, 1000, 10000, 100000.0])
    r = climate_adjusted_quantiles(data, "gumbel", T, confidence_level=95.0, **LP_CLIMATE)

    # Claim 1: climate 97.5% upper is "broadly ~50%" above the baseline POINT
    # estimate at high return periods. Accept 40-60%.
    ratio_pt = r.climate_upper / r.baseline_point - 1.0
    assert np.all(ratio_pt[-3:] > 0.40)
    assert np.all(ratio_pt[-3:] < 0.60)

    # Claim 2: climate 97.5% upper is roughly +30% to +35% above the BASELINE
    # 97.5% upper. The paper's magnitude is ~30-35%; digitization + N make the
    # exact figure and its T-trend uncertain, so assert the magnitude band
    # (25-45%) rather than the precise slope direction.
    ratio_ci = r.climate_upper / r.baseline_upper - 1.0
    assert np.all(ratio_ci > 0.25)
    assert np.all(ratio_ci < 0.45)

    # PMF agreement: climate 97.5% upper at 100,000 yr sits in the neighborhood
    # of the stated PMF (4140). Accept within 15%.
    assert r.climate_upper[-1] == pytest.approx(LP_PMF, rel=0.15)

    # Design/safety range 3500-4500: the climate 97.5% upper for long T should
    # land inside/near this band.
    assert 3300 < r.climate_upper[-1] < 4700


# --------------------------------------------------------------------------- #
# (c) Monte Carlo vs closed form                                              #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("dist", ["gumbel", "lognormal2", "pearson3"])
def test_monte_carlo_agrees_with_closed_form(dist):
    rng = np.random.default_rng(42)
    # A moderately-skewed baseline so Pearson III has a real shape.
    data = rng.gumbel(900, 400, size=60)
    T = np.array([10, 100, 1000, 10000.0])
    d1, d2, t1, t2 = 0.30, 0.18, 0.10, 0.18

    cf = climate_adjusted_quantiles(data, dist, T, d1, d2, t1, t2, confidence_level=90.0)
    mc = mc_climate_adjusted_quantiles(data, dist, T, d1, d2, t1, t2,
                                       confidence_level=90.0, n_sim=20000, random_state=7)

    # Central estimates: closed form and MC use the same exact moment-based
    # parameter mapping. Tight for Gumbel/Log-Normal; Pearson III uses
    # Wilson-Hilferty (closed form) vs scipy's exact quantile (MC), ~1-3% apart.
    np.testing.assert_allclose(mc.climate_point, cf.climate_point, rtol=0.03)
    np.testing.assert_allclose(mc.baseline_point, cf.baseline_point, rtol=0.03)

    # CI bounds: the closed form is a first-order delta method (Gumbel's is
    # exact; Log-Normal and Pearson III are first-order and, for Log-Normal,
    # exponentiate a symmetric log-space interval into an asymmetric physical
    # one). Agreement with the exact Monte Carlo is within simulation noise
    # plus that first-order residual -- 10% covers all three across T here.
    np.testing.assert_allclose(mc.climate_lower, cf.climate_lower, rtol=0.10)
    np.testing.assert_allclose(mc.climate_upper, cf.climate_upper, rtol=0.10)
    np.testing.assert_allclose(mc.baseline_lower, cf.baseline_lower, rtol=0.10)
    np.testing.assert_allclose(mc.baseline_upper, cf.baseline_upper, rtol=0.10)


def test_monte_carlo_climate_ci_wider_than_baseline():
    # Sanity: adding climate uncertainty must widen the interval.
    rng = np.random.default_rng(1)
    data = rng.gumbel(900, 400, size=60)
    T = np.array([100, 1000.0])
    mc = mc_climate_adjusted_quantiles(data, "gumbel", T, 0.30, 0.18, 0.10, 0.18,
                                       n_sim=8000, random_state=3)
    base_width = mc.baseline_upper - mc.baseline_lower
    clim_width = mc.climate_upper - mc.climate_lower
    assert np.all(clim_width > base_width)


def test_to_frame_roundtrip():
    data = np.random.default_rng(0).gumbel(900, 400, size=40)
    r = climate_adjusted_quantiles(data, "gumbel", [100, 1000], 0.3, 0.18, 0.1, 0.18)
    df = r.to_frame()
    assert list(df.columns) == [
        "T", "baseline_point", "baseline_lower", "baseline_upper",
        "climate_point", "climate_lower", "climate_upper",
    ]
    assert len(df) == 2


# --------------------------------------------------------------------------- #
# (d) End-to-end via the io_utils template loader                             #
# --------------------------------------------------------------------------- #
def test_end_to_end_from_climate_inputs_file(tmp_path):
    """Baseline series + a Data/<CaseName>/climate_adjustment/<scenario>.csv
    inputs file flow through resolve_climate_case + load_climate_inputs into
    climate_adjusted_quantiles, matching a direct call with the same numbers."""
    import pandas as pd
    from floodfreq.io_utils import (
        resolve_climate_case, load_climate_inputs, read_series,
    )
    # fake project skeleton
    (tmp_path / "Module").mkdir()
    module_file = tmp_path / "Module" / "run_analysis.py"

    # case-first baseline series at Data/<CaseName>/<CaseName>.csv
    case = tmp_path / "Data" / "DemoCase"
    (case / "climate_adjustment").mkdir(parents=True)
    rng = np.random.default_rng(0)
    q = rng.gumbel(900, 400, size=50)
    pd.DataFrame({"year": np.arange(1970, 1970 + q.size), "Q": q}).to_csv(
        case / "DemoCase.csv", index=False)

    # climate inputs for one scenario
    pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1", "tau2",
                      "confidence_level", "return_periods", "distribution"],
        "value": ["0.30", "0.18", "0.10", "0.18", "95", "100,1000,10000", "gumbel"],
    }).to_csv(case / "climate_adjustment" / "rcp85.csv", index=False)

    paths = resolve_climate_case("DemoCase", "rcp85", module_file)
    Q, years = read_series(paths.baseline_csv)
    # read_series returns (values, years) in THAT order -- guard against a
    # swap: the flood series must be the ~900-scale gumbel draws, not the
    # 1970s year labels. (A value/year unpacking swap silently produces
    # year-valued "floods" and absurd quantiles.)
    assert 300 < float(np.mean(Q)) < 3000, "read_series returned years, not the flood series"
    assert years is not None and int(years[0]) == 1970

    ci = load_climate_inputs(paths.climate_csv)
    r = climate_adjusted_quantiles(
        Q, ci["distribution"], ci["return_periods"],
        ci["delta1"], ci["delta2"], ci["tau1"], ci["tau2"],
        confidence_level=ci["confidence_level"],
    )
    # same as a direct call with identical numbers
    direct = climate_adjusted_quantiles(
        Q, "gumbel", (100, 1000, 10000), 0.30, 0.18, 0.10, 0.18, confidence_level=95.0)
    np.testing.assert_allclose(r.climate_point, direct.climate_point)
    np.testing.assert_allclose(r.climate_upper, direct.climate_upper)
    assert list(r.T) == [100, 1000, 10000]
