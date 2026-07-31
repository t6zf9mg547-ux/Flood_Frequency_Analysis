"""
Data quality checks for an annual-maximum series, run before fitting.

Covers three concerns:
  - stationarity (Mann-Kendall trend test) — the whole flood-frequency
    framework assumes the data are i.i.d. over time (Rao & Hamed, Ch. 1);
    a significant trend means that assumption is questionable.
  - outliers (Grubbs' test, ASTM E178) — a single unusually high/low value
    can distort the sample skewness and drag the whole fit with it.
  - basic input validation — short records, missing/negative values,
    duplicate or non-monotonic years, zero variance.
"""
from __future__ import annotations
import numpy as np
from scipy import stats


def _years_to_numeric(years):
    """Coerce a year column to a numeric array for arithmetic (diffs, spans,
    Sen's slope), while the caller keeps the original labels for display.

    Handles the common formats:
      - already numeric (int/float years)               -> unchanged
      - hydrological-year labels 'YYYY/YYYY' or 'YYYY-YYYY'
                                                        -> the START year (int)
      - plain numeric strings '1950'                    -> 1950
    Returns None if the column can't be interpreted numerically (so callers can
    skip year-based checks rather than crash on an unsupported format)."""
    if years is None:
        return None
    arr = np.asarray(years)
    if np.issubdtype(arr.dtype, np.number):
        return arr.astype(float)
    out = []
    for label in arr.tolist():
        s = str(label).strip()
        # leading run of digits = start year of 'YYYY/YYYY', 'YYYY-YYYY', or a plain year
        i = 0
        while i < len(s) and s[i].isdigit():
            i += 1
        if i == 0:
            return None  # no leading number -> can't interpret
        out.append(float(s[:i]))
    return np.asarray(out, dtype=float)


def sens_slope(x: np.ndarray, t=None) -> dict:
    """
    Theil-Sen slope estimator: the median of all pairwise slopes
    (x_j - x_i) / (t_j - t_i) for i < j. The standard nonparametric
    companion to the Mann-Kendall test — used here purely to draw a
    representative trend line, not as a significance test in itself
    (significance comes from mann_kendall_test).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    t = np.arange(n, dtype=float) if t is None else np.asarray(t, dtype=float)

    slopes = []
    for i in range(n - 1):
        dt = t[i + 1:] - t[i]
        dx = x[i + 1:] - x[i]
        valid = dt != 0
        slopes.extend((dx[valid] / dt[valid]).tolist())

    slope = float(np.median(slopes))
    intercept = float(np.median(x) - slope * np.median(t))
    return {"slope": slope, "intercept": intercept}


def mann_kendall_test(x: np.ndarray, alpha: float = 0.05) -> dict:
    """
    Nonparametric Mann-Kendall trend test (Mann, 1945; Kendall, 1975).

    Returns S (the raw statistic), Z (normal-approximation test statistic,
    tie-corrected), the two-sided p-value, and a plain-language trend label.
    Does not assume any particular distribution for the data, and is the
    standard choice for checking stationarity of hydrologic time series.
    """
    x = np.asarray(x, dtype=float)
    n = x.size

    # S = number of concordant pairs minus discordant pairs
    # S = sum_{i<j} sign(x_j - x_i)  (NOT sign(x_i - x_j) -- direction matters)
    diffs = np.sign(x[None, :] - x[:, None])  # diffs[i, j] = sign(x[j] - x[i])
    iu = np.triu_indices(n, k=1)
    S = diffs[iu].sum()

    # tie correction for the variance
    _, counts = np.unique(x, return_counts=True)
    tie_term = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var_S = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if S > 0:
        Z = (S - 1) / np.sqrt(var_S)
    elif S < 0:
        Z = (S + 1) / np.sqrt(var_S)
    else:
        Z = 0.0

    p_value = 2 * (1 - stats.norm.cdf(abs(Z)))
    significant = p_value < alpha
    if not significant:
        trend = "no significant trend"
    else:
        trend = "significant increasing trend" if S > 0 else "significant decreasing trend"

    return {"S": int(S), "Z": float(Z), "p_value": float(p_value),
            "significant": bool(significant), "trend": trend}


def grubbs_outlier_test(x: np.ndarray, alpha: float = 0.05, log_space: bool = True) -> dict:
    """
    Grubbs' test (ASTM E178) for a single high and a single low outlier,
    applied on the log-transformed series by default (standard practice
    for positively-skewed hydrologic data, matching the space Log-Pearson
    III/Bulletin-17-style outlier testing operates in).

    Tests the single most extreme high value and the single most extreme
    low value independently (not iteratively), so it can miss masked
    outliers (e.g. two adjacent extreme values hiding each other) — a
    reasonable first check, not a substitute for visual inspection of the
    probability plot.
    """
    x = np.asarray(x, dtype=float)
    y = np.log(x) if log_space else x
    n = y.size
    mean, std = y.mean(), y.std(ddof=1)

    # Grubbs' critical value (two-sided, per-tail alpha/(2n))
    t_crit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    G_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit**2 / (n - 2 + t_crit**2))

    i_high = np.argmax(y)
    i_low = np.argmin(y)
    G_high = (y[i_high] - mean) / std
    G_low = (mean - y[i_low]) / std

    return {
        "G_critical": float(G_crit),
        "high_outlier_value": float(x[i_high]), "high_outlier_G": float(G_high),
        "high_outlier_flagged": bool(G_high > G_crit),
        "low_outlier_value": float(x[i_low]), "low_outlier_G": float(G_low),
        "low_outlier_flagged": bool(G_low > G_crit),
        "space": "log" if log_space else "physical",
    }


def validate_series(x: np.ndarray, years=None) -> list:
    """
    Basic input-quality checks. Returns a list of human-readable warning
    strings (empty list if nothing is flagged).
    """
    warnings = []
    x = np.asarray(x, dtype=float)
    n = x.size

    if n < 15:
        warnings.append(f"Very short record (n={n}): parameter estimates, especially for "
                         f"3-parameter distributions, will be highly uncertain.")
    elif n < 30:
        warnings.append(f"Short record (n={n}): treat higher-return-period estimates "
                         f"with extra caution.")

    if np.isnan(x).any():
        warnings.append(f"{int(np.isnan(x).sum())} missing (NaN) value(s) in the series.")
    if np.any(x <= 0):
        warnings.append(f"{int(np.sum(x <= 0))} zero or negative value(s) in the series — "
                         f"invalid for a flood series and for the log-based distributions.")
    if np.nanstd(x) == 0:
        warnings.append("Zero variance: all values are identical.")

    if years is not None:
        years = np.asarray(years)
        if len(years) != n:
            warnings.append("Year column length doesn't match the data column length.")
        else:
            # duplicate check works on the original labels (numeric or text)
            if len(set(years.tolist())) != len(years):
                warnings.append("Duplicate year(s) found in the series.")
            # arithmetic checks need numeric years; hydrological-year labels
            # like 'YYYY/YYYY' are coerced to their start year. If the column
            # can't be interpreted numerically, skip these checks rather than
            # crash.
            ynum = _years_to_numeric(years)
            if ynum is not None:
                if np.any(np.diff(ynum) < 0):
                    warnings.append("Years are not sorted in increasing order.")
                expected_span = ynum.max() - ynum.min() + 1
                if expected_span != n:
                    warnings.append(f"Year range spans {int(expected_span)} years but there are "
                                     f"only {n} values — {int(expected_span - n)} year(s) may be missing.")

    return warnings


def run_all(x: np.ndarray, years=None, alpha: float = 0.05) -> dict:
    """Bundle all checks into one dict, ready for reporting."""
    # Sen's slope needs numeric time; coerce year labels (e.g. 'YYYY/YYYY').
    # Falls back to an integer index if the labels aren't interpretable.
    t_numeric = _years_to_numeric(years)
    return {
        "validation_warnings": validate_series(x, years=years),
        "mann_kendall": mann_kendall_test(x, alpha=alpha),
        "sens_slope": sens_slope(x, t=t_numeric),
        "grubbs": grubbs_outlier_test(x, alpha=alpha),
    }
