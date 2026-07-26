"""
High-level API: FloodFrequencyAnalysis ties together data, plotting
positions, distribution fitting, goodness-of-fit ranking, and confidence
intervals into one convenient object.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .distributions import DISTRIBUTIONS, fit as fit_distribution
from .moments import summarize
from .plotting_positions import empirical_frequency, RECOMMENDED_FOR

# Default fitting method per distribution family — mirrors the logic the
# original Excel workbook used per tab (MLE for Normal/LogNormal, since it
# coincides with MOM there; PWM for the extreme-value/skewed families where
# it is more robust; MOM for Gamma).
DEFAULT_METHOD = {
    "normal": "mle", "lognormal2": "mle",
    "gumbel": "pwm", "gev": "pwm", "exponential": "pwm",
    "gamma2": "mom", "pearson3": "pwm", "logpearson3": "pwm",
    "lognormal3": "pwm",
}


class FloodFrequencyAnalysis:
    def __init__(self, data, station_id: str = "station", years=None):
        self.data = np.asarray(data, dtype=float)
        self.n = self.data.size
        self.station_id = station_id
        self.years = years
        self.fits = {}  # (dist_key, method) -> FitResult

    # ------------------------------------------------------------------ #
    def descriptive_stats(self, plotting_position="weibull") -> dict:
        return summarize(self.data, formula=plotting_position)

    # ------------------------------------------------------------------ #
    def fit(self, dist_key: str, method: str = "pwm", plotting_position="weibull"):
        """Fit one distribution/method combination and cache the result."""
        result = fit_distribution(dist_key, self.data, method=method,
                                   plotting_position=plotting_position)
        self.fits[(dist_key, method)] = result
        return result

    def fit_all(self, methods=None, plotting_position=None, plotting_positions=None):
        """
        Fit every distribution in the registry.

        By default each distribution uses its recommended (method,
        plotting-position) pair:
          - method: MLE for Normal/LogNormal-2p (mathematically ≡ MOM there),
            PWM for Gumbel/GEV/Exponential/Pearson III/LogPearson III/
            LogNormal-3p, MOM for Gamma-2p.
          - plotting position (used only by PWM fits): Blom for Normal/
            LogNormal families, Gringorten for the extreme-value/skewed
            families (see plotting_positions.RECOMMENDED_FOR).

        Overrides:
          methods: {dist_key: method} to override the method for specific
              distributions.
          plotting_positions: {dist_key: formula} to override the plotting
              position for specific distributions.
          plotting_position: a single formula name applied to *every*
              distribution's PWM fit, overriding the per-distribution
              recommendations entirely (useful for sensitivity testing, or
              to reproduce the old "one formula for everything" behaviour).
        """
        methods = methods or {}
        plotting_positions = plotting_positions or {}
        results = {}
        for key in DISTRIBUTIONS:
            method = methods.get(key, DEFAULT_METHOD.get(key, "pwm"))
            if plotting_position is not None:
                pp = plotting_position
            else:
                pp = plotting_positions.get(key, RECOMMENDED_FOR.get(key, "weibull"))
            try:
                results[key] = self.fit(key, method=method, plotting_position=pp)
            except Exception as e:  # pragma: no cover - defensive
                results[key] = e
        return results

    # ------------------------------------------------------------------ #
    def goodness_of_fit_table(self) -> pd.DataFrame:
        """
        AIC / BIC / KS statistic for every distribution fitted so far.

        NOTE: the "plotting_position" column records which formula was
        *requested* for PWM fits, for traceability with the rest of the
        report — it does not affect the fitted parameters themselves. PWM
        parameter estimation always uses the standard unbiased L-moment
        estimator (see distributions.fit() for why). The plotting-position
        choice does affect where empirical points are placed on probability
        plots and the separate descriptive-statistics summary.
        """
        rows = []
        for (key, method), r in self.fits.items():
            if not hasattr(r, "aic"):
                continue
            try:
                ks_d, ks_p = r.ks_statistic(self.data)
            except Exception:
                ks_d, ks_p = np.nan, np.nan
            rows.append({
                "distribution": DISTRIBUTIONS[key]["label"],
                "key": key, "method": method,
                "plotting_position": r.plotting_position or "-",
                "n_params": r.k,
                "loglik": r.loglik(self.data),
                "AIC": r.aic(self.data),
                "BIC": r.bic(self.data),
                "KS_stat": ks_d, "KS_pvalue": ks_p,
            })
        df = pd.DataFrame(rows).sort_values("AIC").reset_index(drop=True)
        return df

    def best_fit(self, criterion="AIC"):
        """Return the (key, method) pair with the lowest AIC/BIC."""
        table = self.goodness_of_fit_table()
        if table.empty:
            raise RuntimeError("No distributions have been fitted yet. Call fit_all() first.")
        best_row = table.sort_values(criterion).iloc[0]
        return best_row["key"], best_row["method"]

    # ------------------------------------------------------------------ #
    def quantile_table(self, return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000),
                        keys=None) -> pd.DataFrame:
        """Flood quantile (design flood) for each fitted distribution at each T."""
        keys = keys or list(self.fits.keys())
        out = {"T": list(return_periods)}
        for (dist_key, method) in keys:
            r = self.fits.get((dist_key, method))
            if r is None:
                continue
            label = f"{DISTRIBUTIONS[dist_key]['label']} ({method.upper()})"
            out[label] = [float(r.quantile(T)) for T in return_periods]
        return pd.DataFrame(out)

    # ------------------------------------------------------------------ #
    def bootstrap_ci(self, dist_key: str, method: str, return_periods,
                      n_boot: int = 1000, alpha: float = 0.05, random_state=None,
                      plotting_position=None):
        """
        Non-parametric bootstrap confidence interval for the T-year flood.
        More robust than the classical normal-approximation CI for skewed
        distributions (GEV, LP3, etc.), at the cost of computation time.

        plotting_position: formula to use when method == "pwm". If not given,
        reuses whatever plotting position was used for the original
        (dist_key, method) fit (so the CI is computed consistently with the
        point estimate), falling back to the per-distribution recommendation
        if that fit hasn't been run yet.
        """
        if plotting_position is None:
            existing = self.fits.get((dist_key, method))
            if existing is not None and existing.plotting_position is not None:
                plotting_position = existing.plotting_position
            else:
                plotting_position = RECOMMENDED_FOR.get(dist_key, "weibull")

        rng = np.random.default_rng(random_state)
        return_periods = np.asarray(return_periods, dtype=float)
        boot_quantiles = np.empty((n_boot, return_periods.size))
        for i in range(n_boot):
            sample = rng.choice(self.data, size=self.n, replace=True)
            try:
                r = fit_distribution(dist_key, sample, method=method,
                                      plotting_position=plotting_position)
                boot_quantiles[i] = r.quantile(return_periods)
            except Exception:
                boot_quantiles[i] = np.nan
        lo = np.nanpercentile(boot_quantiles, 100 * alpha / 2, axis=0)
        hi = np.nanpercentile(boot_quantiles, 100 * (1 - alpha / 2), axis=0)
        med = np.nanpercentile(boot_quantiles, 50, axis=0)
        return pd.DataFrame({"T": return_periods, "lower": lo, "median": med, "upper": hi})

    # ------------------------------------------------------------------ #
    def empirical_points(self, plotting_position=None, dist_key=None):
        """
        (F, T, sorted data) triples for a probability plot. If dist_key is
        given, uses the plotting-position formula recommended for that
        distribution family; otherwise falls back to Weibull.
        """
        formula = plotting_position or RECOMMENDED_FOR.get(dist_key, "weibull")
        x_sorted = np.sort(self.data)
        F = empirical_frequency(self.n, formula=formula, ascending_rank=True)
        T = 1.0 / (1.0 - F)
        return F, T, x_sorted

    def summary(self) -> str:
        stats = self.descriptive_stats()
        lines = [
            f"Station: {self.station_id}",
            f"n = {self.n}",
            f"Mean = {stats['mean']:.2f}, Std = {stats['std']:.2f}, "
            f"CV = {stats['CV']:.3f}, CS = {stats['CS']:.3f}, CK = {stats['CK']:.3f}",
        ]
        if self.fits:
            best_key, best_method = self.best_fit()
            lines.append(f"Best fit by AIC: {DISTRIBUTIONS[best_key]['label']} ({best_method.upper()})")
        return "\n".join(lines)

    # ------------------------------------------------------------------ #
    def generate_recommendation(self, return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000),
                                 criterion="AIC", confidence_level=95.0,
                                 ci_n_boot=1000, random_state=0) -> str:
        """
        Build a plain-text summary that states the recommended distribution,
        the reasoning behind it, and the resulting design-flood quantiles
        with a bootstrap confidence interval. Intended to be both printed
        to the console and written to Output/<CaseName>/summary.txt.

        confidence_level: the confidence level in percent (e.g. 95 for a
            95% CI, 90 for a 90% CI). Must be strictly between 0 and 100.
        """
        if not self.fits:
            raise RuntimeError("No distributions have been fitted yet. Call fit_all() first.")
        if not (0 < confidence_level < 100):
            raise ValueError(f"confidence_level must be between 0 and 100 (got {confidence_level}).")
        alpha = 1.0 - confidence_level / 100.0

        stats = self.descriptive_stats()
        table = self.goodness_of_fit_table().sort_values(criterion).reset_index(drop=True)
        best_row = table.iloc[0]
        best_key, best_method = best_row["key"], best_row["method"]

        # How decisive is the choice? Compare the top candidate to the runner-up.
        runner_up = table.iloc[1] if len(table) > 1 else None
        gap = None
        if runner_up is not None:
            gap = float(runner_up[criterion] - best_row[criterion])

        ks_p = best_row["KS_pvalue"]
        ks_note = ("no evidence against this fit (KS test)" if (pd.notna(ks_p) and ks_p > 0.05)
                   else "the KS test flags some lack of fit — treat design values with extra caution"
                   if pd.notna(ks_p) else "KS test could not be computed")

        q = self.fits[(best_key, best_method)].quantile(np.asarray(return_periods, dtype=float))
        ci = self.bootstrap_ci(best_key, best_method, return_periods,
                                n_boot=ci_n_boot, alpha=alpha, random_state=random_state)

        lines = []
        lines.append("=" * 70)
        lines.append(f"FLOOD FREQUENCY ANALYSIS — SUMMARY & RECOMMENDATION")
        lines.append(f"Station / case: {self.station_id}")
        lines.append("=" * 70)
        lines.append("")
        lines.append(f"Sample size: n = {self.n}")
        lines.append(f"Mean = {stats['mean']:.2f}   Std = {stats['std']:.2f}   "
                      f"CV = {stats['CV']:.3f}   CS (skew) = {stats['CS']:.3f}   CK (kurtosis) = {stats['CK']:.3f}")
        lines.append("")
        lines.append(f"RECOMMENDED DISTRIBUTION: {DISTRIBUTIONS[best_key]['label']}  "
                      f"(fitted by {best_method.upper()}"
                      + (f", requested plotting position: {best_row['plotting_position']}" if best_row['plotting_position'] != "-" else "")
                      + ")")
        lines.append(f"  - Lowest {criterion} among {len(table)} candidates fitted "
                      f"({criterion} = {best_row[criterion]:.2f})")
        if gap is not None:
            verdict = ("a clear margin over the runner-up" if gap > 2
                       else "only a marginal edge over the runner-up — consider both plausible")
            lines.append(f"  - {gap:.2f} {criterion} points ahead of the next-best "
                          f"({runner_up['distribution']}, {runner_up['method'].upper()}): {verdict}")
        lines.append(f"  - KS statistic = {best_row['KS_stat']:.4f}"
                      + (f", p = {ks_p:.3f}" if pd.notna(ks_p) else "")
                      + f" -> {ks_note}")
        lines.append("")
        lines.append("Top candidates (ranked):")
        cols = ["distribution", "method", "plotting_position", "AIC", "BIC", "KS_stat", "KS_pvalue"]
        lines.append(table[cols].head(5).round(3).to_string(index=False))
        lines.append("")
        lines.append(f"Design flood estimates — {DISTRIBUTIONS[best_key]['label']} "
                      f"({best_method.upper()}), with {confidence_level:g}% bootstrap CI "
                      f"({ci_n_boot} resamples):")
        return_periods_arr = np.asarray(return_periods, dtype=float)

        def _flag(T):
            ratio = T / self.n
            if ratio <= 2:
                return ""
            elif ratio <= 10:
                return "extrapolation"
            else:
                return "EXTREME extrapolation"

        design_df = pd.DataFrame({
            "T (years)": ci["T"].astype(int),
            "Q_design": q.round(1),
            "CI_lower": ci["lower"].round(1),
            "CI_upper": ci["upper"].round(1),
            "T / n": (ci["T"] / self.n).round(1),
            "Caution": [_flag(T) for T in return_periods_arr],
        })
        lines.append(design_df.to_string(index=False))
        lines.append("")
        lines.append("Notes:")
        lines.append("  - Ranking is based on AIC/BIC (penalized log-likelihood) plus a")
        lines.append("    Kolmogorov-Smirnov goodness-of-fit check; it does not replace engineering")
        lines.append("    judgement (e.g. regional consistency, physical plausibility of the tail).")
        lines.append(f"  - The {confidence_level:g}% bootstrap CI reflects sampling uncertainty from a")
        lines.append("    finite record only, not model-choice uncertainty or non-stationarity.")
        lines.append("  - PWM-fitted distributions use the standard unbiased L-moment estimator")
        lines.append("    (Hosking & Wallis, 1997) for parameter estimation, regardless of the")
        lines.append("    'plotting position' shown above; that formula affects only where")
        lines.append("    empirical points are drawn on probability plots and the separate")
        lines.append("    descriptive-statistics summary, not the fitted parameters themselves.")
        max_ratio = return_periods_arr.max() / self.n
        if max_ratio > 10:
            lines.append("")
            lines.append(f"  *** EXTRAPOLATION WARNING ***")
            lines.append(f"  The largest return period requested (T = {int(return_periods_arr.max())} years) is "
                         f"{max_ratio:.0f}x the record length (n = {self.n} years).")
            lines.append("  Estimates this far beyond the observed record are dominated by which")
            lines.append("  distribution's tail shape you trust, not by the data itself — different")
            lines.append("  candidate distributions that fit the observed range almost identically well")
            lines.append("  can diverge enormously at T >> n. As a rule of thumb, treat quantile estimates")
            lines.append("  beyond about 2-3x the record length as increasingly uncertain, and beyond")
            lines.append("  ~10x as indicative only, not a number to design against without additional")
            lines.append("  support (regional/pooled analysis, historical/paleoflood information, or a")
            lines.append("  physically-based extreme rainfall-runoff study).")
        elif max_ratio > 2:
            lines.append("")
            lines.append(f"  Note: T = {int(return_periods_arr.max())} years is {max_ratio:.0f}x the record "
                         f"length (n = {self.n}); treat that end of the table with appropriate caution.")
        lines.append("=" * 70)
        return "\n".join(lines)