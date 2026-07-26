"""
Empirical (plotting-position) frequency formulas.

General form (Rao & Hamed, 2000):

    F(r) = (r - a) / (n + b)

where r is the rank of the observation (ascending order, r = 1 .. n)
and n is the sample size. Six standard formulas are provided, matching
the options historically used in the reference Excel workbook.
"""
from __future__ import annotations
import numpy as np

# name -> (a, b)
PLOTTING_POSITIONS = {
    "weibull":    (0.0,   1.0),
    "hazen":      (0.5,   0.0),
    "cunnane":    (0.4,   0.2),
    "gringorten": (0.44,  0.12),
    "hosking":    (0.35,  0.0),
    "blom":       (0.375, 0.25),
}

# Numeric codes, kept for continuity with the original workbook
PLOTTING_POSITION_CODES = {
    1: "hazen",
    2: "cunnane",
    3: "gringorten",
    4: "hosking",
    5: "weibull",
    6: "blom",
}

RECOMMENDED_FOR = {
    "normal": "blom",
    "lognormal2": "blom",
    "lognormal3": "blom",
    "gumbel": "gringorten",
    "gev": "gringorten",
    "exponential": "gringorten",
    "gamma2": "gringorten",
    "pearson3": "gringorten",
    "logpearson3": "gringorten",
}


def resolve_formula(formula) -> str:
    """Accept either a formula name or a legacy numeric code (1-6)."""
    if isinstance(formula, (int, np.integer)):
        try:
            return PLOTTING_POSITION_CODES[int(formula)]
        except KeyError:
            raise ValueError(
                f"Unknown plotting-position code {formula}. "
                f"Valid codes: {PLOTTING_POSITION_CODES}"
            )
    name = str(formula).lower()
    if name not in PLOTTING_POSITIONS:
        raise ValueError(
            f"Unknown plotting-position formula '{formula}'. "
            f"Valid names: {list(PLOTTING_POSITIONS)}"
        )
    return name


def empirical_frequency(n: int, formula="weibull", ascending_rank: bool = True) -> np.ndarray:
    """
    Return the vector of non-exceedance probabilities F for ranks 1..n.

    Parameters
    ----------
    n : sample size
    formula : plotting-position name (e.g. "blom") or legacy code (1-6)
    ascending_rank : if True, rank 1 = smallest value (non-exceedance F increases
        with the data value, as is standard). If False, rank 1 = largest value
        (i.e. r is a "descending" rank, common in some hydrology texts when
        working with exceedance probability / return period directly).
    """
    name = resolve_formula(formula)
    a, b = PLOTTING_POSITIONS[name]
    r = np.arange(1, n + 1)
    if not ascending_rank:
        r = r[::-1]
    F = (r - a) / (n + b)
    return F


def return_period(F: np.ndarray) -> np.ndarray:
    """Return period T (years) from a non-exceedance probability F."""
    F = np.asarray(F, dtype=float)
    return 1.0 / (1.0 - F)


def non_exceedance_from_T(T) -> np.ndarray:
    """Non-exceedance probability F from a return period T (years)."""
    T = np.asarray(T, dtype=float)
    return 1.0 - 1.0 / T
