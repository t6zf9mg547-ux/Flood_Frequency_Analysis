"""
Regional (pooled) flood-frequency analysis: the index-flood method
(Dalrymple, 1960) with L-moment-based regional statistics
(Hosking & Wallis, 1997; see also Rao & Hamed, 2000, Ch. 9, for the
notation this module follows most closely).

Index-flood idea: pool N hydrologically-similar stations so that the
*shape* of the flood distribution -- a dimensionless regional growth
curve g(F), normalized to a mean of 1 -- is estimated from many
station-years combined, while each station keeps its own scale, the
"index flood" mu_i (here: the station's own sample mean, l1). The design
flood at station i for return period T is then

    Q_i(T) = mu_i * g(F),      F = 1 - 1/T

This module implements the steps in the order Hosking & Wallis recommend
applying them:

1. `station_lmoments`  -- unbiased sample L-moments/ratios per station
2. `discordancy`        -- D_i: flag stations statistically inconsistent
                            with the rest of the group
3. `heterogeneity`      -- H: is the group homogeneous enough to pool,
                            via Monte-Carlo simulation of synthetic
                            homogeneous regions?
4. `zstatistics`        -- Z-statistic: which candidate regional
                            distribution family fits the pooled data best?
5. `fit_growth_curve`   -- fit the chosen family to the pooled,
                            record-length-weighted regional L-moments
6. `RegionalGrowthCurve.station_quantile(T, mu_i)` -- station design flood

`run_regional_analysis()` ties these six steps into one call.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import integrate

import lmoments3.distr as ld
from lmoments3 import lmom_ratios as _sample_lmom_ratios

from . import data_quality as _data_quality


# ---------------------------------------------------------------------------
# Candidate regional distribution families (Hosking & Wallis, 1997, Ch. 4).
# All five are 3-parameter families in lmoments3's `distr` module and are
# fittable directly from the first three L-moments (l1, l2, t3).
# ---------------------------------------------------------------------------
CANDIDATE_FAMILIES = {
    "glo": ld.glo,   # Generalized Logistic
    "gev": ld.gev,   # Generalized Extreme Value
    "gno": ld.gno,   # Generalized Normal (~ LogNormal, 3-parameter)
    "pe3": ld.pe3,   # Pearson Type III
    "gpa": ld.gpa,   # Generalized Pareto
}

FAMILY_LABELS = {
    "glo": "Generalized Logistic (GLO)",
    "gev": "Generalized Extreme Value (GEV)",
    "gno": "Generalized Normal / LogNormal-3p (GNO)",
    "pe3": "Pearson Type III (PE3)",
    "gpa": "Generalized Pareto (GPA)",
}

# Hosking & Wallis (1997), Table 3.1: critical values of the discordancy
# measure D_i for a region of N sites, below which a site is not flagged
# as discordant.
DISCORDANCY_CRITICAL = {
    5: 1.333, 6: 1.648, 7: 1.917, 8: 2.140, 9: 2.329, 10: 2.491,
    11: 2.632, 12: 2.757, 13: 2.869, 14: 2.971, 15: 3.061,
}
DISCORDANCY_CRITICAL_LARGE_N = 3.000  # N > 15 (asymptotic chi-sq(3) value)


def discordancy_critical_value(n_sites: int) -> float:
    """Critical D_i value for screening a region of n_sites stations."""
    if n_sites in DISCORDANCY_CRITICAL:
        return DISCORDANCY_CRITICAL[n_sites]
    if n_sites > 15:
        return DISCORDANCY_CRITICAL_LARGE_N
    # N < 5: S (the 3x3 sum-of-squares matrix) is poorly conditioned with
    # so few points and D_i is not a very meaningful screen; fall back to
    # the smallest tabulated critical value as a soft floor.
    return DISCORDANCY_CRITICAL[5]


# ---------------------------------------------------------------------------
# 1. At-site L-moments
# ---------------------------------------------------------------------------
@dataclass
class StationLMoments:
    """Unbiased sample L-moments/ratios for one station in a pooling group."""
    name: str
    n: int
    mean: float   # index flood: the station's own sample mean (== l1)
    l1: float
    l2: float
    t: float      # L-CV  = l2 / l1
    t3: float     # L-skewness
    t4: float     # L-kurtosis
    t5: float = float("nan")


def station_lmoments(name: str, x) -> StationLMoments:
    """
    Unbiased sample L-moments (Hosking & Wallis, 1997) for one station's
    annual-maximum series.

    Uses the standard *unbiased* PWM/L-moment estimator -- the same one
    `floodfreq.distributions.fit(..., method="pwm")` uses for single-station
    parameter fitting -- NOT the plotting-position-weighted approximation in
    `floodfreq.moments.summarize()`, which is a display/diagnostic
    convention specific to the single-station descriptive-statistics table
    and is documented there as not matching the unbiased estimator.
    """
    x = np.asarray(x, dtype=float)
    if x.size < 5:
        raise ValueError(
            f"Station '{name}' has only {x.size} value(s); at least 5 years "
            f"of record are needed to compute L-moment ratios up to t4.")
    if np.any(~np.isfinite(x)):
        raise ValueError(f"Station '{name}': data contains NaN/Inf values.")
    if np.any(x <= 0):
        raise ValueError(f"Station '{name}': data must be strictly positive.")

    l1, l2, t3, t4, t5 = _sample_lmom_ratios(x, nmom=5)
    t = l2 / l1 if l1 else float("nan")
    return StationLMoments(name=name, n=x.size, mean=float(np.mean(x)),
                            l1=float(l1), l2=float(l2), t=float(t),
                            t3=float(t3), t4=float(t4), t5=float(t5))


def stations_table(stations: Sequence[StationLMoments]) -> pd.DataFrame:
    """Tabulate a list of StationLMoments (for CSV output / inspection)."""
    return pd.DataFrame([{
        "station": s.name, "n": s.n, "mean": s.mean,
        "l1": s.l1, "l2": s.l2, "t_LCV": s.t, "t3_Lskew": s.t3, "t4_Lkurt": s.t4,
    } for s in stations])


def station_data_quality(stations: Sequence[StationLMoments], station_data: dict,
                          station_years: dict | None = None, alpha: float = 0.05) -> pd.DataFrame:
    """
    Per-station data quality checks, run BEFORE pooling: a trending or
    outlier-contaminated station can quietly bias the pooled growth curve,
    the same way it would bias a single-station fit. Reuses
    `floodfreq.data_quality` -- the same Mann-Kendall stationarity test,
    Grubbs' outlier test, and basic input validation the single-station
    tool runs -- applied here to each station's own record independently.

    `station_years`, as returned by `io_utils.load_region_stations()`, maps
    station name -> years array (or None). When available for a station,
    this is passed through to the Mann-Kendall/Sen's-slope/validation
    checks, which then: (a) can detect missing years within a station's
    nominal span (the same "year range spans N years but there are only M
    values" warning the single-station tool gives), and (b) report Sen's
    slope in real per-calendar-year units. Mann-Kendall's trend
    classification and p-value are rank-based and unaffected either way
    (they don't depend on the actual spacing between observations, only
    their order) -- years change what you learn ABOUT the time axis, not
    whether a trend is detected. For any station without years (either
    `station_years` wasn't passed, or that station's own CSV had no
    detectable year column), record order is used as a stand-in for time,
    same as the single-station tool's fallback when no year column is
    present.
    """
    rows = []
    for s in stations:
        x = np.asarray(station_data[s.name], dtype=float)
        years = station_years.get(s.name) if station_years else None
        dq = _data_quality.run_all(x, years=years, alpha=alpha)
        mk, gr = dq["mann_kendall"], dq["grubbs"]
        rows.append({
            "station": s.name, "n": x.size,
            "years_available": years is not None,
            "mann_kendall_trend": mk["trend"],
            "mann_kendall_p_value": round(mk["p_value"], 4),
            "mann_kendall_significant": mk["significant"],
            "sens_slope_per_year": round(dq["sens_slope"]["slope"], 4),
            "grubbs_high_outlier_value": gr["high_outlier_value"],
            "grubbs_high_outlier_flagged": gr["high_outlier_flagged"],
            "grubbs_low_outlier_value": gr["low_outlier_value"],
            "grubbs_low_outlier_flagged": gr["low_outlier_flagged"],
            "validation_warnings": " | ".join(dq["validation_warnings"]) or "none",
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 2. Discordancy measure D_i (Hosking & Wallis, 1997, sec. 3.2)
# ---------------------------------------------------------------------------
def discordancy(stations: Sequence[StationLMoments]) -> pd.DataFrame:
    """
    Discordancy measure D_i for each station.

    Each station is represented by u_i = [t_i, t3_i, t4_i]^T (L-CV,
    L-skewness, L-kurtosis). With ubar the unweighted average across the
    N stations and S = sum_i (u_i - ubar)(u_i - ubar)^T,

        D_i = (N / 3) * (u_i - ubar)^T S^-1 (u_i - ubar)

    A station with D_i above the tabulated critical value (see
    `discordancy_critical_value`) is unusually different from the rest of
    the group in its L-moment ratios and is a candidate for exclusion or
    closer scrutiny before pooling.
    """
    n_sites = len(stations)
    if n_sites < 4:
        raise ValueError(
            f"Discordancy needs at least 4 stations to invert the 3x3 "
            f"covariance matrix S (got {n_sites}). Add more stations to "
            f"the region.")

    U = np.array([[s.t, s.t3, s.t4] for s in stations])  # (N, 3)
    ubar = U.mean(axis=0)
    dev = U - ubar
    S = dev.T @ dev  # 3x3, sum (not averaged) of outer products
    try:
        S_inv = np.linalg.inv(S)
    except np.linalg.LinAlgError as e:
        raise ValueError(
            "Discordancy: the 3x3 L-moment-ratio covariance matrix S is "
            "singular (e.g. two stations with identical t/t3/t4, or fewer "
            "distinct points than dimensions) and cannot be inverted."
        ) from e

    D = np.array([
        (n_sites / 3.0) * (dev[i] @ S_inv @ dev[i]) for i in range(n_sites)
    ])
    crit = discordancy_critical_value(n_sites)
    return pd.DataFrame({
        "station": [s.name for s in stations],
        "n": [s.n for s in stations],
        "t_LCV": U[:, 0], "t3_Lskew": U[:, 1], "t4_Lkurt": U[:, 2],
        "D_i": D,
        "D_critical": crit,
        "discordant": D > crit,
    })


# ---------------------------------------------------------------------------
# Regional-average (record-length-weighted) L-moment ratios
# ---------------------------------------------------------------------------
def regional_average_ratios(stations: Sequence[StationLMoments]) -> dict:
    """
    Record-length-weighted regional-average L-moment ratios (Hosking &
    Wallis, 1997, eq. 3.1-ish weighting convention used throughout their
    heterogeneity/goodness-of-fit statistics):

        t_R  = sum_i n_i t_i  / sum_i n_i
        t3_R = sum_i n_i t3_i / sum_i n_i
        t4_R = sum_i n_i t4_i / sum_i n_i
    """
    n = np.array([s.n for s in stations], dtype=float)
    t = np.array([s.t for s in stations])
    t3 = np.array([s.t3 for s in stations])
    t4 = np.array([s.t4 for s in stations])
    W = n.sum()
    return {
        "t_R": float(np.sum(n * t) / W),
        "t3_R": float(np.sum(n * t3) / W),
        "t4_R": float(np.sum(n * t4) / W),
        "n_total": float(W),
        "n_sites": len(stations),
    }


def _weighted_V1(t_sim: np.ndarray, n: np.ndarray, t_R_sim: float) -> float:
    """Weighted standard deviation of at-site L-CV around the regional
    average -- the heterogeneity measure's H1 "V statistic" (Hosking &
    Wallis, 1997, eq. 4.4)."""
    return float(np.sqrt(np.sum(n * (t_sim - t_R_sim) ** 2) / np.sum(n)))


# ---------------------------------------------------------------------------
# 3 & 4. Monte-Carlo simulation shared by heterogeneity (H) and the
# distribution-selection Z-statistics
# ---------------------------------------------------------------------------
def _fit_kappa_region(t_R: float, t3_R: float, t4_R: float):
    """
    Fit the 4-parameter kappa distribution to the regional-average
    L-moments (l1=1, l2=t_R, t3_R, t4_R), as recommended by Hosking &
    Wallis (1997) for simulating synthetic homogeneous regions -- kappa
    nests GLO/GEV/GNO/GPA/Gumbel as special/limiting cases, so it does not
    bias the simulation toward any one candidate family.

    Falls back to a 3-parameter GEV fit (matching l1, l2, t3_R only, i.e.
    ignoring t4_R) if the (t3_R, t4_R) pair falls outside the kappa
    distribution's feasible region -- this can happen for small, noisy
    regions -- with a warning, since GEV is a reasonable general-purpose
    stand-in and the simulation still reproduces the regional L-CV/L-skew.
    """
    try:
        params = ld.kap.lmom_fit(lmom_ratios=[1.0, t_R, t3_R, t4_R])
        frozen = ld.kap(*params.values())
        # sanity check: the fit must round-trip to a valid distribution
        _ = frozen.ppf(0.5)
        if not np.isfinite(_):
            raise ValueError("kappa fit produced a non-finite median")
        return frozen, "kap"
    except Exception as e:
        warnings.warn(
            f"Could not fit a 4-parameter kappa distribution to the "
            f"regional L-moments (t_R={t_R:.4f}, t3_R={t3_R:.4f}, "
            f"t4_R={t4_R:.4f}): {e}. Falling back to GEV (matching l1, l2, "
            f"t3_R only) for the Monte-Carlo simulation.", stacklevel=2)
        params = ld.gev.lmom_fit(lmom_ratios=[1.0, t_R, t3_R])
        frozen = ld.gev(*params.values())
        return frozen, "gev"


def _simulate_regions(stations: Sequence[StationLMoments], n_sim: int = 500,
                       seed: int | None = None) -> pd.DataFrame:
    """
    Simulate n_sim synthetic "homogeneous" regions: for each replicate,
    draw a series of length n_i for every station from the kappa (or GEV
    fallback) distribution fitted to the observed regional-average
    L-moments, then recompute the same regional summary statistics
    (weighted t_R, t3_R, t4_R, and the V1 heterogeneity statistic) on the
    simulated data. Returns one row per replicate.
    """
    obs = regional_average_ratios(stations)
    kappa_dist, family_used = _fit_kappa_region(obs["t_R"], obs["t3_R"], obs["t4_R"])
    n = np.array([s.n for s in stations], dtype=float)
    rng = np.random.default_rng(seed)

    rows = []
    for _ in range(n_sim):
        t_sim = np.empty(len(stations))
        t3_sim = np.empty(len(stations))
        t4_sim = np.empty(len(stations))
        for i, ni in enumerate(n.astype(int)):
            x_sim = kappa_dist.rvs(size=ni, random_state=rng)
            l1, l2, t3, t4, _t5 = _sample_lmom_ratios(x_sim, nmom=5)
            t_sim[i] = l2 / l1 if l1 else np.nan
            t3_sim[i] = t3
            t4_sim[i] = t4
        W = n.sum()
        t_R_sim = float(np.sum(n * t_sim) / W)
        t3_R_sim = float(np.sum(n * t3_sim) / W)
        t4_R_sim = float(np.sum(n * t4_sim) / W)
        V1_sim = _weighted_V1(t_sim, n, t_R_sim)
        rows.append((t_R_sim, t3_R_sim, t4_R_sim, V1_sim))

    df = pd.DataFrame(rows, columns=["t_R", "t3_R", "t4_R", "V1"])
    df.attrs["simulation_family"] = family_used
    return df


@dataclass
class HeterogeneityResult:
    H1: float
    V1_obs: float
    V1_sim_mean: float
    V1_sim_std: float
    n_sim: int
    simulation_family: str
    interpretation: str


def heterogeneity(stations: Sequence[StationLMoments], n_sim: int = 500,
                   seed: int | None = None) -> HeterogeneityResult:
    """
    Heterogeneity measure H1 (Hosking & Wallis, 1997, sec. 4.3.3): compares
    the observed weighted dispersion of at-site L-CV, V1_obs, against the
    distribution of V1 across n_sim Monte-Carlo-simulated homogeneous
    regions (mu_V, sigma_V):

        H1 = (V1_obs - mu_V) / sigma_V

    Conventional interpretation (Hosking & Wallis, 1997):
        H1 < 1   -- the region is "acceptably homogeneous"
        1 <= H1 < 2 -- "possibly heterogeneous"
        H1 >= 2  -- "definitely heterogeneous"
    """
    if len(stations) < 2:
        raise ValueError("Heterogeneity needs at least 2 stations.")
    n = np.array([s.n for s in stations], dtype=float)
    t = np.array([s.t for s in stations])
    obs = regional_average_ratios(stations)
    V1_obs = _weighted_V1(t, n, obs["t_R"])

    sim = _simulate_regions(stations, n_sim=n_sim, seed=seed)
    mu_V, sigma_V = sim["V1"].mean(), sim["V1"].std(ddof=1)
    H1 = (V1_obs - mu_V) / sigma_V if sigma_V > 0 else float("nan")

    if not np.isfinite(H1):
        interp = "undefined (zero simulated spread; try more replicates)"
    elif H1 < 1:
        interp = "acceptably homogeneous"
    elif H1 < 2:
        interp = "possibly heterogeneous"
    else:
        interp = "definitely heterogeneous"

    return HeterogeneityResult(
        H1=float(H1), V1_obs=float(V1_obs), V1_sim_mean=float(mu_V),
        V1_sim_std=float(sigma_V), n_sim=n_sim,
        simulation_family=sim.attrs["simulation_family"], interpretation=interp)


# ---------------------------------------------------------------------------
# 4. Regional distribution selection: Z-statistic
# ---------------------------------------------------------------------------
def _population_lmoments_from_ppf(ppf, nmom: int = 4, eps: float = 1e-9) -> list[float]:
    """
    Population L-moments 1..nmom of a distribution from its quantile
    function Q(u) = ppf(u), via Hosking's (1990) shifted-Legendre-
    polynomial integral representation:

        lambda_{r+1} = integral_0^1  Q(u) * P*_r(u) du

    with shifted Legendre polynomials P*_0..P*_3 = 1, 2u-1, 6u^2-6u+1,
    20u^3-30u^2+12u-1. Used here to get the *theoretical* t4 of a
    candidate family at a given t3 (see `theoretical_tau4`), by first
    fitting that family's parameters to (l1=0, l2=1, t3=target) -- which
    for these 3-parameter L-moment-fittable families is a direct,
    non-iterative solve -- and then integrating the resulting quantile
    function for its true (population) t4.
    """
    legendre = [
        lambda u: np.ones_like(u),
        lambda u: 2 * u - 1,
        lambda u: 6 * u ** 2 - 6 * u + 1,
        lambda u: 20 * u ** 3 - 30 * u ** 2 + 12 * u - 1,
    ]
    lam = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=integrate.IntegrationWarning)
        for r in range(nmom):
            def integrand(u, _r=r):
                return ppf(u) * legendre[_r](u)
            val, _err = integrate.quad(integrand, eps, 1 - eps, limit=200)
            lam.append(val)
    return lam


def theoretical_tau4(tau3: float, family: str) -> float:
    """
    Theoretical (population) L-kurtosis t4 of `family` at a given
    L-skewness t3 -- the t4(t3) relationship used as the reference curve
    for that family on an L-moment ratio diagram, and as the target value
    in the distribution-selection Z-statistic.

    `family` must be a key of CANDIDATE_FAMILIES. Computed generically by
    fitting the family to (l1=0, l2=1, t3) and numerically integrating its
    quantile function for the true population t4 -- this is exact to
    integration tolerance rather than relying on published polynomial
    approximations, and is validated in the test suite against the
    closed-form GLO/GPA relations (t4 = (1+5t3^2)/6 and
    t4 = t3(1+5t3)/(5+t3) respectively).
    """
    if family not in CANDIDATE_FAMILIES:
        raise ValueError(f"Unknown family '{family}'. Choose from {list(CANDIDATE_FAMILIES)}.")
    dist = CANDIDATE_FAMILIES[family]
    params = dist.lmom_fit(lmom_ratios=[0.0, 1.0, tau3])
    frozen = dist(*params.values())
    l1, l2, l3, l4 = _population_lmoments_from_ppf(frozen.ppf, nmom=4)
    return float(l4 / l2)


def zstatistics(stations: Sequence[StationLMoments],
                 candidates: Sequence[str] = ("glo", "gev", "gno", "pe3", "gpa"),
                 n_sim: int = 500, seed: int | None = None) -> pd.DataFrame:
    """
    Z-statistic for regional distribution selection (Hosking & Wallis,
    1997, sec. 5.2). For each candidate family:

        Z_DIST = (t4_DIST(t3_R) - t4_R + B4) / sigma4

    where t4_DIST(t3_R) is that family's theoretical L-kurtosis at the
    observed regional-average L-skewness t3_R; B4 and sigma4 are the bias
    and standard deviation of the simulated regional-average L-kurtosis
    across n_sim Monte-Carlo-simulated homogeneous regions (reusing the
    same simulation as `heterogeneity`). |Z_DIST| <= 1.64 is conventionally
    taken as an adequate fit at roughly the 90% level.
    """
    obs = regional_average_ratios(stations)
    sim = _simulate_regions(stations, n_sim=n_sim, seed=seed)
    B4 = float(sim["t4_R"].mean() - obs["t4_R"])
    sigma4 = float(sim["t4_R"].std(ddof=1))

    rows = []
    for fam in candidates:
        t4_dist = theoretical_tau4(obs["t3_R"], fam)
        Z = (t4_dist - obs["t4_R"] + B4) / sigma4 if sigma4 > 0 else float("nan")
        rows.append({
            "family": fam, "label": FAMILY_LABELS[fam],
            "t4_theoretical": t4_dist, "t4_regional": obs["t4_R"],
            "Z": Z, "acceptable": bool(np.isfinite(Z) and abs(Z) <= 1.64),
        })
    df = pd.DataFrame(rows).sort_values("Z", key=lambda s: s.abs()).reset_index(drop=True)
    df.attrs["B4"] = B4
    df.attrs["sigma4"] = sigma4
    df.attrs["t3_R"] = obs["t3_R"]
    df.attrs["t4_R"] = obs["t4_R"]
    return df


def recommend_family(zstat_table: pd.DataFrame) -> str:
    """
    Pick a regional distribution family from a `zstatistics()` table: the
    acceptable (|Z| <= 1.64) family with the smallest |Z|, or -- if none
    are acceptable -- the smallest-|Z| family overall with a caveat that
    no candidate fits well.
    """
    if zstat_table.empty:
        raise ValueError("Empty Z-statistic table.")
    ordered = zstat_table.sort_values("Z", key=lambda s: s.abs())
    acceptable = ordered[ordered["acceptable"]]
    chosen = acceptable.iloc[0] if not acceptable.empty else ordered.iloc[0]
    return str(chosen["family"])


# ---------------------------------------------------------------------------
# 5 & 6. Regional growth curve
# ---------------------------------------------------------------------------
@dataclass
class RegionalGrowthCurve:
    """
    A fitted, dimensionless regional growth curve g(F), normalized so that
    its mean (l1) is 1: station design floods are g(F) * (station's own
    index flood).
    """
    family: str
    params: dict
    t_R: float
    t3_R: float
    n_sites: int
    n_total: int

    def _frozen(self):
        return CANDIDATE_FAMILIES[self.family](*self.params.values())

    def ppf(self, F):
        """Growth factor at non-exceedance probability F."""
        return self._frozen().ppf(np.asarray(F, dtype=float))

    def quantile(self, T):
        """Growth factor at return period T (years)."""
        T = np.asarray(T, dtype=float)
        return self.ppf(1.0 - 1.0 / T)

    def station_quantile(self, T, index_flood: float):
        """Design flood at return period T for a station with the given
        index flood (typically that station's own sample mean)."""
        return self.quantile(T) * index_flood

    @property
    def label(self) -> str:
        return FAMILY_LABELS[self.family]


def fit_growth_curve(stations: Sequence[StationLMoments], family: str) -> RegionalGrowthCurve:
    """
    Fit the regional growth curve: `family`'s parameters from the pooled,
    record-length-weighted regional L-moments (l1=1, l2=t_R, t3=t3_R) --
    the index-flood assumption is encoded exactly here by fixing l1=1, so
    the fitted curve is a pure growth factor, not a magnitude.
    """
    if family not in CANDIDATE_FAMILIES:
        raise ValueError(f"Unknown family '{family}'. Choose from {list(CANDIDATE_FAMILIES)}.")
    obs = regional_average_ratios(stations)
    dist = CANDIDATE_FAMILIES[family]
    params = dist.lmom_fit(lmom_ratios=[1.0, obs["t_R"], obs["t3_R"]])
    return RegionalGrowthCurve(family=family, params=dict(params), t_R=obs["t_R"],
                                t3_R=obs["t3_R"], n_sites=len(stations),
                                n_total=int(obs["n_total"]))


def growth_curve_quantile_table(growth_curve: RegionalGrowthCurve,
                                 stations: Sequence[StationLMoments],
                                 return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000)) -> pd.DataFrame:
    """Design-flood table: rows = return periods, one column per station
    (station design flood = growth curve * that station's index flood)
    plus the dimensionless growth factor itself."""
    T = np.asarray(return_periods, dtype=float)
    g = growth_curve.quantile(T)
    data = {"T_years": T, "growth_factor": g}
    for s in stations:
        data[s.name] = g * s.mean
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# Confidence intervals: growth-curve shape uncertainty (Monte-Carlo,
# Hosking & Wallis, 1997, sec. 5.3's accuracy-assessment approach) combined
# with at-site index-flood uncertainty (bootstrap of each station's own mean)
# ---------------------------------------------------------------------------
def _simulate_growth_curve_quantiles(stations: Sequence[StationLMoments], family: str,
                                      return_periods, n_sim: int = 500,
                                      seed: int | None = None) -> np.ndarray:
    """
    Distribution of the regional growth curve's quantile at each requested
    return period, across n_sim Monte-Carlo-simulated homogeneous regions
    (Hosking & Wallis, 1997, sec. 5.3): reuses the same kappa-distribution
    simulation as `heterogeneity()`/`zstatistics()`, refits `family`'s
    growth curve to each replicate's simulated regional L-moments the same
    way `fit_growth_curve()` does for the observed data, and evaluates its
    quantile at each T. This captures uncertainty in the *shape* of the
    growth curve (how well the pooled record pins down t_R, t3_R) -- not
    uncertainty in any individual station's index flood, which is handled
    separately (see `station_quantile_ci`).

    Returns an (n_sim, len(return_periods)) array of growth factors; a
    replicate that fails to fit (rare, e.g. an infeasible simulated
    L-moment combination) is filled with NaN and excluded by the
    percentile step downstream.
    """
    sim = _simulate_regions(stations, n_sim=n_sim, seed=seed)
    dist = CANDIDATE_FAMILIES[family]
    T = np.asarray(return_periods, dtype=float)
    F = 1.0 - 1.0 / T
    out = np.full((len(sim), len(T)), np.nan)
    for i, row in sim.iterrows():
        try:
            params = dist.lmom_fit(lmom_ratios=[1.0, row["t_R"], row["t3_R"]])
            frozen = dist(*params.values())
            out[i] = frozen.ppf(F)
        except Exception:
            continue
    return out


def station_quantile_ci(result: "RegionalAnalysisResult", station_name: str,
                         return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000),
                         n_sim: int = 500, seed: int | None = None,
                         alpha: float = 0.05) -> pd.DataFrame:
    """
    Confidence interval for one station's regional design flood at each
    return period, combining two independent Monte-Carlo-simulated
    sources of uncertainty (the regional counterpart of the single-station
    tool's `bootstrap_ci()`):

    1. Regional growth-curve shape uncertainty -- `n_sim` synthetic
       homogeneous regions simulated from the kappa distribution fitted to
       the observed regional L-moments (same simulation as
       `heterogeneity()`/`zstatistics()`), growth curve refit to each.
    2. At-site index-flood uncertainty -- a nonparametric bootstrap of
       *this* station's own record (`n_sim` resamples with replacement of
       its own n_i values), since the index flood is just this station's
       sample mean and its own record is the only information available
       about it.

    The two are combined by pairing independent Monte-Carlo replicates,
    Q_sim = g_sim(T) * mu_boot_sim -- the practical way to combine two
    uncertainty sources when their exact joint distribution isn't
    tractable analytically, in the same spirit as combining at-site and
    regional error components in Bulletin 17B / FEH-style accuracy
    statements, but by simulation here rather than an analytical RMSE
    formula. Percentiles of Q_sim give the CI; note this treats the two
    sources as independent, which is standard practice but does mean the
    band can be a bit narrow if a station's own record was itself part of
    what pinned down the regional shape (unavoidable with only one method
    of this kind and not corrected for here).
    """
    station = next((s for s in result.stations if s.name == station_name), None)
    if station is None:
        raise ValueError(f"Unknown station '{station_name}'. Available: "
                          f"{[s.name for s in result.stations]}")
    x = np.asarray(result.station_data[station_name], dtype=float)

    g_sim = _simulate_growth_curve_quantiles(result.stations, result.growth_curve.family,
                                              return_periods, n_sim=n_sim, seed=seed)
    boot_seed = None if seed is None else seed + 1
    rng = np.random.default_rng(boot_seed)
    mu_boot = np.array([rng.choice(x, size=x.size, replace=True).mean() for _ in range(g_sim.shape[0])])

    Q_sim = g_sim * mu_boot[:, None]
    T = np.asarray(return_periods, dtype=float)
    Q_point = result.growth_curve.station_quantile(T, station.mean)
    lower = np.nanpercentile(Q_sim, 100 * alpha / 2, axis=0)
    upper = np.nanpercentile(Q_sim, 100 * (1 - alpha / 2), axis=0)
    return pd.DataFrame({"station": station_name, "T_years": T, "Q_design": Q_point,
                          "CI_lower": lower, "CI_upper": upper})


def regional_quantile_ci(result: "RegionalAnalysisResult",
                          return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000),
                          n_sim: int = 500, seed: int | None = None,
                          alpha: float = 0.05) -> pd.DataFrame:
    """`station_quantile_ci()` for every station in the region, stacked
    into one long-format table (station, T_years, Q_design, CI_lower,
    CI_upper). Reuses one shared growth-curve simulation across all
    stations (only the at-site bootstrap differs per station), so this is
    much cheaper than calling `station_quantile_ci()` once per station."""
    g_sim = _simulate_growth_curve_quantiles(result.stations, result.growth_curve.family,
                                              return_periods, n_sim=n_sim, seed=seed)
    T = np.asarray(return_periods, dtype=float)
    rows = []
    for i, s in enumerate(result.stations):
        x = np.asarray(result.station_data[s.name], dtype=float)
        boot_seed = None if seed is None else seed + 1 + i
        rng = np.random.default_rng(boot_seed)
        mu_boot = np.array([rng.choice(x, size=x.size, replace=True).mean() for _ in range(g_sim.shape[0])])
        Q_sim = g_sim * mu_boot[:, None]
        Q_point = result.growth_curve.station_quantile(T, s.mean)
        lower = np.nanpercentile(Q_sim, 100 * alpha / 2, axis=0)
        upper = np.nanpercentile(Q_sim, 100 * (1 - alpha / 2), axis=0)
        rows.append(pd.DataFrame({"station": s.name, "T_years": T, "Q_design": Q_point,
                                   "CI_lower": lower, "CI_upper": upper}))
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------
@dataclass
class RegionalAnalysisResult:
    region_name: str
    stations: list
    stations_df: pd.DataFrame
    data_quality_df: pd.DataFrame
    discordancy_df: pd.DataFrame
    heterogeneity_result: HeterogeneityResult
    zstat_df: pd.DataFrame
    chosen_family: str
    growth_curve: RegionalGrowthCurve
    quantile_table: pd.DataFrame
    station_data: dict = field(default_factory=dict)  # {station_name: raw values array}
    station_years: dict = field(default_factory=dict)  # {station_name: years array or None}


def run_regional_analysis(region_name: str, station_data: dict,
                           candidates: Sequence[str] = ("glo", "gev", "gno", "pe3", "gpa"),
                           n_sim: int = 500, seed: int | None = None,
                           return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000),
                           family: str | None = None,
                           station_years: dict | None = None) -> RegionalAnalysisResult:
    """
    Run the full six-step regional (index-flood, L-moment-based) pipeline
    for one pooling group.

    station_data: {station_name: 1-D array-like of annual maxima, ...}
    station_years: optional {station_name: 1-D array-like of years, ...}
        (or with some/all values None), as returned by
        `io_utils.load_region_stations()`. When given, feeds the
        Mann-Kendall/Sen's-slope/validation checks in
        `station_data_quality()` real calendar time instead of falling
        back to record order -- see that function's docstring for what
        this does and doesn't change.
    family: force a specific regional distribution instead of using
        `recommend_family()` on the Z-statistic table (e.g. to match a
        published study, or to compare families deliberately).
    """
    if len(station_data) < 4:
        raise ValueError(
            f"Regional analysis needs at least 4 stations to compute "
            f"discordancy (got {len(station_data)}). Add more station CSVs "
            f"to this region's Data/Regional/<RegionName>/ folder.")

    stations = [station_lmoments(name, x) for name, x in station_data.items()]
    stations_df = stations_table(stations)
    dq_df = station_data_quality(stations, station_data, station_years=station_years)
    disc_df = discordancy(stations)
    het = heterogeneity(stations, n_sim=n_sim, seed=seed)
    zdf = zstatistics(stations, candidates=candidates, n_sim=n_sim, seed=seed)
    chosen = family if family is not None else recommend_family(zdf)
    growth_curve = fit_growth_curve(stations, chosen)
    q_table = growth_curve_quantile_table(growth_curve, stations, return_periods=return_periods)

    return RegionalAnalysisResult(
        region_name=region_name, stations=stations, stations_df=stations_df,
        data_quality_df=dq_df, discordancy_df=disc_df, heterogeneity_result=het, zstat_df=zdf,
        chosen_family=chosen, growth_curve=growth_curve, quantile_table=q_table,
        station_data={k: np.asarray(v, dtype=float) for k, v in station_data.items()},
        station_years=dict(station_years) if station_years else {})
