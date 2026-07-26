"""
Descriptive statistics used in flood frequency analysis:
ordinary (product) moments and probability-weighted moments (PWM).

Mirrors the "Stat" sheet of the reference workbook (m1..m4, CV, CS, CK,
b0..b3), plus L-moments (a linear transform of the PWMs) which are the
preferred basis for parameter estimation in modern practice
(Hosking & Wallis, 1997).
"""
from __future__ import annotations
import numpy as np
from .plotting_positions import empirical_frequency


def ordinary_moments(x: np.ndarray) -> dict:
    """Sample mean, variance-type moments (population convention, /n),
    coefficient of variation, skewness and kurtosis."""
    x = np.asarray(x, dtype=float)
    n = x.size
    m1 = x.mean()
    m2 = np.mean((x - m1) ** 2)
    m3 = np.mean((x - m1) ** 3)
    m4 = np.mean((x - m1) ** 4)
    std = np.sqrt(m2)
    cv = std / m1
    # unbiased-ish skew/kurtosis estimators (common hydrology convention)
    g1 = (n ** 2 / ((n - 1) * (n - 2))) * np.sum((x - m1) ** 3) / n / std ** 3 if n > 2 else np.nan
    g2 = m4 / std ** 4
    return {
        "n": n, "mean": m1, "m2": m2, "m3": m3, "m4": m4,
        "std": std, "CV": cv, "CS": g1, "CK": g2,
    }


def probability_weighted_moments(x: np.ndarray, formula="weibull") -> dict:
    """
    Sample PWMs b0..b3, using the given plotting-position formula for the
    empirical non-exceedance probability F used as the weight.

        b_r = (1/n) * sum_i  x_(i) * F_i^r          i = 1..n, x sorted ascending
    """
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    F = empirical_frequency(n, formula=formula, ascending_rank=True)
    b0 = np.mean(x)
    b1 = np.mean(x * F)
    b2 = np.mean(x * F ** 2)
    b3 = np.mean(x * F ** 3)
    return {"b0": b0, "b1": b1, "b2": b2, "b3": b3}


def pwm_to_lmoments(b: dict) -> dict:
    """Convert PWMs (b0..b3) to L-moments (l1..l4) and L-moment ratios (t3, t4)."""
    b0, b1, b2, b3 = b["b0"], b["b1"], b["b2"], b["b3"]
    l1 = b0
    l2 = 2 * b1 - b0
    l3 = 6 * b2 - 6 * b1 + b0
    l4 = 20 * b3 - 30 * b2 + 12 * b1 - b0
    t3 = l3 / l2 if l2 else np.nan
    t4 = l4 / l2 if l2 else np.nan
    return {"l1": l1, "l2": l2, "l3": l3, "l4": l4, "t3": t3, "t4": t4}


def summarize(x: np.ndarray, formula="weibull") -> dict:
    """Full descriptive-statistics bundle for a sample: moments + PWM + L-moments."""
    out = ordinary_moments(x)
    b = probability_weighted_moments(x, formula=formula)
    out.update(b)
    out.update(pwm_to_lmoments(b))
    return out
