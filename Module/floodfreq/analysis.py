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
from . import data_quality as dq

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


# -- Plain-language quality assessments, for readers who aren't statisticians -- #

def _fit_assessment(ks_pvalue: float) -> str:
    """KS test has a rigorous p-value; use it (not AD_stat, which has none
    here -- see FitResult.anderson_darling_statistic) for a Good/Warning call."""
    if pd.isna(ks_pvalue):
        return "Not assessed (KS test could not be computed)"
    return "Good: no evidence against this fit" if ks_pvalue > 0.05 else \
           "Warning: KS test rejects this fit at the 5% level"


def _weight_assessment(weight: float) -> str:
    if weight >= 0.5:
        return "Strong relative support"
    elif weight >= 0.2:
        return "Moderate relative support"
    elif weight >= 0.05:
        return "Weak relative support"
    else:
        return "Negligible relative support"


def _confidence_assessment(top_weight: float) -> str:
    if top_weight >= 0.7:
        return "Good: one distribution is clearly favored"
    elif top_weight >= 0.4:
        return "Caution: moderate confidence -- other distributions remain plausible"
    else:
        return "Warning: no distribution stands out -- rely on the model-averaged " \
               "estimate rather than a single 'best' fit"


def _agreement_assessment(relative_spread_pct: float) -> str:
    if relative_spread_pct < 2:
        return "Good: candidate distributions agree closely"
    elif relative_spread_pct < 5:
        return "Caution: moderate disagreement between candidate distributions"
    else:
        return "Warning: candidate distributions disagree substantially -- treat as highly uncertain"


class FloodFrequencyAnalysis:
    def __init__(self, data, station_id: str = "station", years=None):
        self.data = np.asarray(data, dtype=float)
        self._validate_input(self.data, years)
        self.n = self.data.size
        self.station_id = station_id
        self.years = np.asarray(years) if years is not None else None
        self.fits = {}  # (dist_key, method) -> FitResult

    @staticmethod
    def _validate_input(data: np.ndarray, years=None):
        """
        Hard failures only -- structurally invalid input that would
        otherwise crash later with a confusing error deep inside a fit
        routine. Softer statistical concerns (short record, outliers,
        borderline trend) are intentionally NOT raised here; those are
        reported as warnings via data_quality() instead, since a short or
        unusual record is still analyzable, just with more caution.
        """
        if data.size == 0:
            raise ValueError("No data provided: the series is empty.")
        if np.isnan(data).any():
            n_nan = int(np.isnan(data).sum())
            raise ValueError(
                f"{n_nan} missing (NaN) value(s) in the data. Remove or fill them "
                f"before fitting -- flood frequency distributions cannot be fit "
                f"through missing values.")
        if np.isinf(data).any():
            raise ValueError("Infinite value(s) found in the data -- check the input file "
                              "for a parsing error or a sentinel value (e.g. -999) that "
                              "wasn't meant to be read as real data.")
        if np.any(data <= 0):
            n_bad = int(np.sum(data <= 0))
            raise ValueError(
                f"{n_bad} zero or negative value(s) in the data. A flood/flow series should "
                f"be strictly positive -- this usually means a data-entry error or a missing-"
                f"value sentinel (e.g. -999, 0) that needs to be cleaned before fitting. "
                f"(Zero/negative values would otherwise silently produce NaN deep inside the "
                f"log-based distributions -- LogNormal, Log-Pearson III -- rather than a clear error.)")
        if data.size < 5:
            raise ValueError(
                f"Only {data.size} data point(s) provided; at least 5 are needed to "
                f"attempt any distribution fit meaningfully (and even that is very "
                f"little -- results with n < 15-20 should be treated with real caution).")
        if years is not None:
            years_arr = np.asarray(years)
            if years_arr.size != data.size:
                raise ValueError(
                    f"years has {years_arr.size} entries but data has {data.size} -- "
                    f"they must be the same length and correspond index-for-index.")

    # ------------------------------------------------------------------ #
    def data_quality(self, alpha: float = 0.05) -> dict:
        """
        Run stationarity (Mann-Kendall), outlier (Grubbs), and basic input
        validation checks on the series. See floodfreq.data_quality for
        details on each test. Cached after the first call.
        """
        if not hasattr(self, "_dq_cache"):
            self._dq_cache = dq.run_all(self.data, years=self.years, alpha=alpha)
        return self._dq_cache

    def descriptive_stats(self, plotting_position="weibull") -> dict:
        return summarize(self.data, formula=plotting_position)

    # ------------------------------------------------------------------ #
    def fit(self, dist_key: str, method: str = "pwm", plotting_position="weibull",
            regional_skew: float = None, regional_mse: float = 0.302):
        """Fit one distribution/method combination and cache the result."""
        result = fit_distribution(dist_key, self.data, method=method,
                                   plotting_position=plotting_position,
                                   regional_skew=regional_skew, regional_mse=regional_mse)
        self.fits[(dist_key, method)] = result
        return result

    def fit_all(self, methods=None, plotting_position=None, plotting_positions=None,
                regional_skew=None, regional_mse=0.302):
        """
        Fit every distribution in the registry (see class-level docstring
        for the default method/plotting-position choices and overrides).

        regional_skew / regional_mse: if regional_skew is given, ALSO fits
        Pearson III and Log-Pearson III using the Bulletin 17B weighted-skew
        procedure (regional_skew.py), as additional candidates alongside
        their standard PWM fits -- both are then ranked side by side by
        AIC/BIC/KS/AD, so the data (not an assumption) decides whether
        incorporating regional skew actually improves the fit here.
        regional_mse defaults to 0.302 (the Bulletin 17B national skew map's
        documented MSE) but a region-specific study will usually be more
        precise -- supply it explicitly if you have one.
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

        if regional_skew is not None:
            for key in ("pearson3", "logpearson3"):
                try:
                    results[f"{key}_weighted_skew"] = self.fit(
                        key, method="mom_weighted_skew",
                        regional_skew=regional_skew, regional_mse=regional_mse)
                except Exception as e:  # pragma: no cover - defensive
                    results[f"{key}_weighted_skew"] = e

        return results

    # ------------------------------------------------------------------ #
    def goodness_of_fit_table(self) -> pd.DataFrame:
        """
        AIC / BIC / KS statistic / Anderson-Darling statistic for every
        distribution fitted so far.

        NOTE: the "plotting_position" column records which formula was
        *requested* for PWM fits, for traceability with the rest of the
        report — it does not affect the fitted parameters themselves. PWM
        parameter estimation always uses the standard unbiased L-moment
        estimator (see distributions.fit() for why). The plotting-position
        choice does affect where empirical points are placed on probability
        plots and the separate descriptive-statistics summary.

        AD_stat (Anderson-Darling) is tail-weighted, unlike KS which weights
        all parts of the distribution equally -- more relevant for flood
        design values, which live in the tail. It has no p-value here (see
        FitResult.anderson_darling_statistic); use it as a relative ranking
        criterion, lower is better.
        """
        rows = []
        for (key, method), r in self.fits.items():
            if not hasattr(r, "aic"):
                continue
            try:
                ks_d, ks_p = r.ks_statistic(self.data)
            except Exception:
                ks_d, ks_p = np.nan, np.nan
            try:
                ad_stat = r.anderson_darling_statistic(self.data)
            except Exception:
                ad_stat = np.nan
            rows.append({
                "distribution": DISTRIBUTIONS[key]["label"],
                "key": key, "method": method,
                "plotting_position": r.plotting_position or "-",
                "n_params": r.k,
                "loglik": r.loglik(self.data),
                "AIC": r.aic(self.data),
                "BIC": r.bic(self.data),
                "KS_stat": ks_d, "KS_pvalue": ks_p,
                "AD_stat": ad_stat,
                "fit_assessment": _fit_assessment(ks_p),
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
    def akaike_weights(self, criterion: str = "AIC") -> pd.DataFrame:
        """
        Akaike weights (Burnham & Anderson, 2002) for every distribution
        fitted so far: w_i = exp(-Delta_i/2) / sum_j exp(-Delta_j/2), where
        Delta_i = AIC_i - min(AIC). Each w_i is interpretable as the
        (approximate) probability that distribution i is the best model
        among the candidate set, given the data -- the basis for
        model-averaged quantiles rather than betting everything on a
        single AIC "winner".
        """
        table = self.goodness_of_fit_table()
        if table.empty:
            raise RuntimeError("No distributions have been fitted yet. Call fit_all() first.")
        delta = table[criterion] - table[criterion].min()
        w = np.exp(-delta / 2)
        w = w / w.sum()
        out = table[["distribution", "key", "method", criterion]].copy()
        out["delta"] = delta
        out["akaike_weight"] = w
        out["weight_assessment"] = out["akaike_weight"].apply(_weight_assessment)
        return out.sort_values("akaike_weight", ascending=False).reset_index(drop=True)

    def model_averaged_quantile_table(self, return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000),
                                       criterion: str = "AIC") -> pd.DataFrame:
        """
        Akaike-weighted multi-model-averaged design flood: at each return
        period, the weighted average of every fitted distribution's
        quantile, weighted by Akaike weight. Also reports the weighted
        between-model standard deviation (how much the candidates disagree,
        weighted by how plausible each one is) -- a measure of *model-
        selection* uncertainty, distinct from (and complementary to) the
        within-model sampling uncertainty captured by bootstrap_ci().
        """
        weights_table = self.akaike_weights(criterion=criterion)
        return_periods = np.asarray(return_periods, dtype=float)

        rows = []
        for T in return_periods:
            q_vals, w_vals = [], []
            for _, row in weights_table.iterrows():
                r = self.fits.get((row["key"], row["method"]))
                if r is None:
                    continue
                q_vals.append(float(r.quantile(T)))
                w_vals.append(row["akaike_weight"])
            q_vals = np.asarray(q_vals)
            w_vals = np.asarray(w_vals)
            w_vals = w_vals / w_vals.sum()  # renormalize in case any fit failed
            q_avg = float(np.sum(w_vals * q_vals))
            between_model_sd = float(np.sqrt(np.sum(w_vals * (q_vals - q_avg) ** 2)))
            relative_spread_pct = 100 * between_model_sd / q_avg if q_avg else np.nan
            rows.append({"T": T, "Q_model_averaged": q_avg, "between_model_sd": between_model_sd,
                        "relative_spread_pct": relative_spread_pct,
                        "agreement_assessment": _agreement_assessment(relative_spread_pct)})
        return pd.DataFrame(rows)

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

        lines.append("-" * 70)
        lines.append("DATA QUALITY CHECKS")
        lines.append("-" * 70)
        dqr = self.data_quality()
        if dqr["validation_warnings"]:
            for w in dqr["validation_warnings"]:
                lines.append(f"  ! {w}")
        else:
            lines.append("  No input-validation issues found (record length, missing/negative "
                         "values, year sequence).")
        mk = dqr["mann_kendall"]
        lines.append(f"  Stationarity (Mann-Kendall trend test): {mk['trend']} "
                     f"(Z = {mk['Z']:.2f}, p = {mk['p_value']:.3f})")
        if mk["significant"]:
            lines.append("    -> A significant trend means the i.i.d./stationarity assumption "
                         "behind this whole analysis is questionable; consider a non-stationary "
                         "method or investigate the cause (land use, climate, regulation changes).")
        gr = dqr["grubbs"]
        if gr["high_outlier_flagged"]:
            lines.append(f"  Possible HIGH outlier flagged (Grubbs' test, log-space): "
                         f"{gr['high_outlier_value']:.1f} (G = {gr['high_outlier_G']:.2f} > "
                         f"critical {gr['G_critical']:.2f})")
        if gr["low_outlier_flagged"]:
            lines.append(f"  Possible LOW outlier flagged (Grubbs' test, log-space): "
                         f"{gr['low_outlier_value']:.1f} (G = {gr['low_outlier_G']:.2f} > "
                         f"critical {gr['G_critical']:.2f})")
        if not gr["high_outlier_flagged"] and not gr["low_outlier_flagged"]:
            lines.append(f"  No outliers flagged (Grubbs' test, log-space; critical G = "
                         f"{gr['G_critical']:.2f}).")
        lines.append("")

        weighted_skew_fits = {k: r for k, r in self.fits.items()
                              if r.method == "mom_weighted_skew" and r.regional_skew_info}
        if weighted_skew_fits:
            lines.append("-" * 70)
            lines.append("REGIONAL SKEW (Bulletin 17B weighted skew)")
            lines.append("-" * 70)
            for (key, method), r in weighted_skew_fits.items():
                info = r.regional_skew_info
                lines.append(f"  {DISTRIBUTIONS[key]['label']}: station skew = {info['station_skew']:.3f} "
                             f"(MSE = {info['station_mse']:.3f}), regional skew = {info['regional_skew']:.3f} "
                             f"(MSE = {info['regional_mse']:.3f})")
                lines.append(f"    -> weighted skew = {info['weighted_skew']:.3f}")
                if info["review_flag"]:
                    lines.append(f"    ! Warning: {info['review_flag']}")
                else:
                    lines.append(f"    Good: station and regional skew are reasonably consistent "
                                 f"(within 0.5).")
            lines.append("  Note: this is the classical Bulletin 17B weighted-skew procedure, not")
            lines.append("  the full Bulletin 17C Expected Moments Algorithm (EMA). It is included")
            lines.append("  as an additional candidate above/below, ranked alongside the standard")
            lines.append("  PWM fits by AIC/BIC/KS/AD -- check whether it actually improved the fit")
            lines.append("  rather than assuming it automatically does.")
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
                      + f" -> {_fit_assessment(ks_p)}")
        weights_table_preview = self.akaike_weights(criterion=criterion)
        top_weight = float(weights_table_preview["akaike_weight"].iloc[0])
        lines.append(f"  - Confidence in this pick (Akaike weight = {top_weight:.0%} among "
                     f"{len(table)} candidates): {_confidence_assessment(top_weight)}")
        lines.append("")
        lines.append("Top candidates (ranked):")
        cols = ["distribution", "method", "plotting_position", "AIC", "BIC", "KS_stat", "KS_pvalue",
                "AD_stat", "fit_assessment"]
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

        lines.append("-" * 70)
        lines.append("MODEL-AVERAGED DESIGN FLOOD (Akaike-weighted across all fitted distributions)")
        lines.append("-" * 70)
        lines.append("Rather than betting entirely on the single AIC-best distribution above, this")
        lines.append("blends every fitted distribution's quantile, weighted by its Akaike weight")
        lines.append("(~ the probability it is the best model among the candidates). "
                     "'between_model_sd'")
        lines.append("is how much the candidates disagree at that T, weighted by plausibility -- a")
        lines.append("model-selection-uncertainty companion to the bootstrap CI above (which only")
        lines.append("captures sampling uncertainty within a single chosen model).")
        weights_table = self.akaike_weights(criterion=criterion)
        lines.append("")
        lines.append("Akaike weights:")
        lines.append(weights_table[["distribution", "method", "akaike_weight", "weight_assessment"]]
                     .round(3).to_string(index=False))
        ma_table = self.model_averaged_quantile_table(return_periods=return_periods, criterion=criterion)
        lines.append("")
        lines.append(ma_table.rename(columns={"T": "T (years)"}).round(1).to_string(index=False))
        lines.append("")
        lines.append("(Good = models agree closely; Caution = moderate disagreement; "
                     "Warning = substantial disagreement -- read the design value with that in mind.)")
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
