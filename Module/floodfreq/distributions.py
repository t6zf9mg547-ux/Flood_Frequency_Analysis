"""
Unified interface to the candidate distributions used in flood frequency
analysis. Wraps scipy.stats (for MLE fitting and as the numerical engine)
and lmoments3 (for probability-weighted-moment / L-moment fitting).

Each entry in DISTRIBUTIONS describes one distribution family and how to
fit it with each of the three classical methods:
    - "mom"  : method of ordinary moments
    - "mle"  : maximum likelihood
    - "pwm"  : probability-weighted moments / L-moments

Two distributions (LogNormal-2p and Log-Pearson III) are handled as a
log-transform of an existing distribution rather than as a distinct
scipy/lmoments3 object.
"""
from __future__ import annotations
import numpy as np
import scipy.stats as ss
import lmoments3.distr as ld
from scipy import optimize

from .moments import ordinary_moments, probability_weighted_moments, pwm_to_lmoments


class FitResult:
    """A fitted distribution, ready to produce quantiles / cdf / pdf."""

    def __init__(self, key: str, method: str, params: tuple, n: int, transform=None,
                 plotting_position=None):
        self.key = key
        self.method = method
        self.params = tuple(float(p) for p in params)
        self.n = n
        self.transform = transform  # None, "log10", or "ln"
        self.plotting_position = plotting_position  # only meaningful when method == "pwm"
        self.regional_skew_info = None  # populated externally for method == "mom_weighted_skew"
        self.spec = DISTRIBUTIONS[key]
        self.k = len(self.params)  # number of fitted parameters

    # -- core distribution calls, in the *physical* (untransformed) scale --
    def _frozen(self):
        return self.spec["lmom_obj"](*self.params)

    def ppf(self, F):
        F = np.asarray(F, dtype=float)
        q = self._frozen().ppf(F)
        if self.transform == "log10":
            q = 10.0 ** q
        elif self.transform == "ln":
            q = np.exp(q)
        return q

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        if self.transform == "log10":
            x = np.log10(x)
        elif self.transform == "ln":
            x = np.log(x)
        return self._frozen().cdf(x)

    def pdf(self, x):
        x = np.asarray(x, dtype=float)
        xt = x
        if self.transform == "log10":
            xt = np.log10(x)
        elif self.transform == "ln":
            xt = np.log(x)
        p = self._frozen().pdf(xt)
        # Jacobian correction back to the physical scale
        if self.transform == "log10":
            p = p / (x * np.log(10))
        elif self.transform == "ln":
            p = p / x
        return p

    def quantile(self, T):
        """Flood magnitude for return period T (years)."""
        T = np.asarray(T, dtype=float)
        F = 1.0 - 1.0 / T
        return self.ppf(F)

    def loglik(self, data):
        """
        Log-likelihood on the *physical* (untransformed) scale, so AIC/BIC
        are comparable across transformed and untransformed distributions.
        For y = transform(x), log f_X(x) = log f_Y(y) - log|dy/dx|.
        """
        data = np.asarray(data, dtype=float)
        if self.transform == "log10":
            y = np.log10(data)
            jac = np.log(data) + np.log(np.log(10))
        elif self.transform == "ln":
            y = np.log(data)
            jac = np.log(data)
        else:
            y = data
            jac = 0.0
        return float(np.sum(self._frozen().logpdf(y) - jac))

    def aic(self, data):
        return 2 * self.k - 2 * self.loglik(data)

    def bic(self, data):
        return self.k * np.log(self.n) - 2 * self.loglik(data)

    def ks_statistic(self, data):
        data = np.asarray(data, dtype=float)
        if self.transform == "log10":
            data = np.log10(data)
        elif self.transform == "ln":
            data = np.log(data)
        d, p = ss.kstest(data, self._frozen().cdf)
        return d, p

    def anderson_darling_statistic(self, data):
        """
        General Anderson-Darling A^2 statistic (Anderson & Darling, 1954),
        computed via the probability integral transform against this
        distribution's own fitted CDF:

            A^2 = -n - (1/n) * sum_i (2i-1) * [ln F(x_(i)) + ln(1 - F(x_(n+1-i)))]

        More sensitive to tail discrepancies than the Kolmogorov-Smirnov
        statistic, which weights all parts of the distribution equally --
        relevant here since flood design values live in the tail.

        NOTE: unlike ks_statistic(), this does NOT return a p-value. Formal
        Anderson-Darling critical values are distribution-specific and were
        derived assuming known (not estimated) parameters; the standard
        corrections for estimated parameters (Stephens, 1974) only cover a
        handful of distributions (normal, exponential, Weibull, etc.), not
        the full set used here (GEV, Pearson III, ...). A^2 is therefore
        reported as an additional *relative* ranking criterion alongside
        AIC/BIC/KS -- lower is better -- not as a formal hypothesis test.
        """
        data = np.sort(np.asarray(data, dtype=float))
        n = data.size
        F = self.cdf(data)
        eps = 1e-12
        F = np.clip(F, eps, 1 - eps)
        i = np.arange(1, n + 1)
        S = np.sum((2 * i - 1) * (np.log(F) + np.log(1 - F[::-1])))
        A2 = -n - S / n
        return float(A2)

    def __repr__(self):
        return (f"<FitResult {self.spec['label']} method={self.method} "
                f"params={np.round(self.params, 4).tolist()}>")


def _transform_data(x, transform):
    if transform == "log10":
        return np.log10(x)
    elif transform == "ln":
        return np.log(x)
    return x


def _mom_normal(x):
    m = ordinary_moments(x)
    return (m["mean"], m["std"])


def _mom_gumbel(x):
    m = ordinary_moments(x)
    scale = np.sqrt(6) * m["std"] / np.pi
    loc = m["mean"] - 0.5772156649 * scale
    return (loc, scale)


def _mom_exponential(x):
    m = ordinary_moments(x)
    loc = m["mean"] - m["std"]
    scale = m["std"]
    return (loc, scale)


def _mom_gamma(x):
    m = ordinary_moments(x)
    shape = (m["mean"] / m["std"]) ** 2
    scale = m["std"] ** 2 / m["mean"]
    return (shape, 0.0, scale)


def _mom_pearson3(x):
    m = ordinary_moments(x)
    return (m["CS"], m["mean"], m["std"])


DISTRIBUTIONS = {
    "normal": {
        "label": "Normal",
        "lmom_obj": ld.nor,
        "scipy_obj": ss.norm,
        "n_params": 2,
        "fit_mom": _mom_normal,
        "transform": None,
    },
    "lognormal2": {
        "label": "LogNormal (2-parameter)",
        "lmom_obj": ld.nor,          # fitted on ln(x)
        "scipy_obj": ss.norm,
        "n_params": 2,
        "fit_mom": _mom_normal,      # applied to ln(x) by the caller
        "transform": "ln",
    },
    "gumbel": {
        "label": "Gumbel (EV1)",
        "lmom_obj": ld.gum,
        "scipy_obj": ss.gumbel_r,
        "n_params": 2,
        "fit_mom": _mom_gumbel,
        "transform": None,
    },
    "gev": {
        "label": "Generalized Extreme Value (GEV)",
        "lmom_obj": ld.gev,
        "scipy_obj": ss.genextreme,
        "n_params": 3,
        "fit_mom": None,  # no simple closed-form MOM for GEV
        "transform": None,
    },
    "exponential": {
        "label": "Exponential",
        "lmom_obj": ld.exp,
        "scipy_obj": ss.expon,
        "n_params": 2,
        "fit_mom": _mom_exponential,
        "transform": None,
    },
    "gamma2": {
        "label": "Gamma (2-parameter)",
        "lmom_obj": ld.gam,
        "scipy_obj": lambda: None,   # scipy.gamma MLE needs floc=0, handled specially
        "n_params": 2,
        "fit_mom": _mom_gamma,
        "transform": None,
    },
    "pearson3": {
        "label": "Pearson Type III",
        "lmom_obj": ld.pe3,
        "scipy_obj": ss.pearson3,
        "n_params": 3,
        "fit_mom": _mom_pearson3,
        "transform": None,
    },
    "logpearson3": {
        "label": "Log-Pearson Type III",
        "lmom_obj": ld.pe3,           # fitted on log10(x)
        "scipy_obj": ss.pearson3,
        "n_params": 3,
        "fit_mom": _mom_pearson3,     # applied to log10(x) by the caller
        "transform": "log10",
    },
    "lognormal3": {
        "label": "LogNormal (3-parameter) / Generalized Normal",
        "lmom_obj": ld.gno,
        "scipy_obj": None,            # no direct scipy MLE equivalent
        "n_params": 3,
        "fit_mom": None,
        "transform": None,
    },
}


def available_distributions():
    return {k: v["label"] for k, v in DISTRIBUTIONS.items()}


def fit(key: str, data: np.ndarray, method: str = "pwm", plotting_position="weibull",
        regional_skew: float = None, regional_mse: float = 0.302) -> FitResult:
    """
    Fit a distribution to a 1-D array of annual maxima.

    method: "mom", "mle", "pwm", or "mom_weighted_skew" (Pearson III /
        Log-Pearson III only -- see regional_skew.py). regional_skew and
        regional_mse are only used by "mom_weighted_skew"; regional_skew
        MUST be supplied by the caller (there is no sensible default --
        it has to come from a published regional study for your area).
    plotting_position: only used by "pwm" (which empirical-frequency formula
        weights the sample probability-weighted moments).
    """
    spec = DISTRIBUTIONS[key]
    x = np.asarray(data, dtype=float)
    n = x.size
    transform = spec["transform"]
    xt = _transform_data(x, transform)

    method = method.lower()
    regional_skew_info = None

    if method == "mom":
        if spec["fit_mom"] is None:
            raise ValueError(f"Method of moments is not available for '{key}'.")
        params = spec["fit_mom"](xt)

    elif method == "mom_weighted_skew":
        if key not in ("pearson3", "logpearson3"):
            raise ValueError("'mom_weighted_skew' is only defined for pearson3/logpearson3 "
                              "(it operates on the skew parameter of a Pearson-family fit).")
        if regional_skew is None:
            raise ValueError("regional_skew must be supplied for method='mom_weighted_skew' "
                              "-- it comes from a published regional study for your area, "
                              "not something this tool can determine on its own.")
        from .regional_skew import weighted_skew
        station_skew, mean, std = _mom_pearson3(xt)
        regional_skew_info = weighted_skew(station_skew, n, regional_skew, regional_mse)
        params = (regional_skew_info["weighted_skew"], mean, std)

    elif method == "mle":
        scipy_obj = spec["scipy_obj"]
        if key == "gamma2":
            shape, loc, scale = ss.gamma.fit(xt, floc=0)
            params = (shape, loc, scale)
        elif scipy_obj is None:
            # generic numerical MLE via the lmoments3/rv_continuous object,
            # started from the PWM solution for stability
            start = spec["lmom_obj"].lmom_fit(xt)
            x0 = list(start.values())

            def neg_ll(theta):
                try:
                    return -np.sum(spec["lmom_obj"].logpdf(xt, *theta))
                except Exception:
                    return np.inf

            res = optimize.minimize(neg_ll, x0=x0, method="Nelder-Mead")
            params = tuple(res.x)
        else:
            params = scipy_obj.fit(xt)

    elif method == "pwm":
        from .plotting_positions import resolve_formula
        resolve_formula(plotting_position)  # validate only
        # NOTE ON plotting_position AND PWM FITTING:
        # PARAMETER ESTIMATION always uses lmoments3's standard *unbiased*
        # L-moment estimator (Hosking & Wallis, 1997) here, regardless of
        # `plotting_position`. This was verified against the reference
        # workbook: e.g. Gumbel PWM on station P1009 matches the workbook's
        # a=37.4567, b=167.0191 to 4 decimal places using the unbiased
        # estimator, but does NOT match if the sample PWMs are instead
        # computed by weighting each order statistic with a chosen
        # plotting-position formula (e.g. Gringorten) as in moments.py's
        # probability_weighted_moments(). In other words: the workbook's own
        # "Fréquence empirique" selector does not feed into its PWM parameter
        # estimation either — it is a display choice, not a fitting choice.
        # `plotting_position` here therefore only affects:
        #   (a) where empirical data points are placed on probability plots
        #       (see analysis.empirical_points / plots.probability_plot), and
        #   (b) the descriptive PWM/L-moment summary in moments.summarize(),
        #       which mirrors the workbook's "Stat" sheet.
        # It is recorded on the FitResult for that reason (traceability with
        # the rest of the report), not because it changes these parameters.
        d = spec["lmom_obj"].lmom_fit(xt)
        params = tuple(d.values())

    else:
        raise ValueError(f"Unknown fitting method '{method}'. "
                          f"Use 'mom', 'mle', 'pwm', or 'mom_weighted_skew'.")

    result = FitResult(key, method, params, n, transform=transform,
                        plotting_position=plotting_position if method == "pwm" else None)
    result.regional_skew_info = regional_skew_info
    return result