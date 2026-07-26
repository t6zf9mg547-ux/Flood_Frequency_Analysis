"""
Plotting helpers for flood frequency analysis: probability plots (data vs.
fitted distribution on a return-period axis), and an L-moment ratio diagram
for visually comparing candidate distributions.
"""
from __future__ import annotations
import matplotlib

matplotlib.use("Agg")  # non-interactive backend: safe for headless script runs

import numpy as np
import matplotlib.pyplot as plt

from .distributions import DISTRIBUTIONS
from .moments import summarize


def probability_plot(ffa, dist_keys=None, methods=None, ax=None, log_x=True, max_T=10000):
    """
    Plot the empirical data points (return period vs. flood magnitude)
    together with one or more fitted distributions' quantile curves,
    extrapolated out to max_T years.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    keys = dist_keys or [k for (k, m) in ffa.fits.keys()]
    keys = list(dict.fromkeys(keys))  # unique, preserve order

    F, T, x_sorted = ffa.empirical_points(dist_key=keys[0] if keys else None)
    ax.scatter(T, x_sorted, s=25, color="black", zorder=5, label="Observed data")

    T_smooth = np.geomspace(max(1.01, T.min()), max(T.max(), max_T), 300)
    F_smooth = 1 - 1 / T_smooth

    for key in keys:
        for (k2, m2), r in ffa.fits.items():
            if k2 != key:
                continue
            label = f"{DISTRIBUTIONS[key]['label']} ({m2.upper()})"
            ax.plot(T_smooth, r.ppf(F_smooth), lw=1.8, label=label)

    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel("Return period T (years)")
    ax.set_ylabel("Flood magnitude")
    ax.set_title(f"Flood frequency curve — {ffa.station_id}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    return ax


def quantile_plot_with_ci(ffa, dist_key, method, return_periods=None,
                           n_boot=500, alpha=0.05, ax=None, random_state=0):
    """Fitted quantile curve with a bootstrap confidence band."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    return_periods = return_periods if return_periods is not None else \
        np.geomspace(1.05, 10000, 100)

    r = ffa.fits[(dist_key, method)]
    q = r.quantile(return_periods)

    ci = ffa.bootstrap_ci(dist_key, method, return_periods,
                           n_boot=n_boot, alpha=alpha, random_state=random_state)

    ci_pct = 100 * (1 - alpha)
    F, T, x_sorted = ffa.empirical_points(dist_key=dist_key)
    ax.scatter(T, x_sorted, s=25, color="black", zorder=5, label="Observed data")
    ax.plot(return_periods, q, color="crimson", lw=2,
            label=f"{DISTRIBUTIONS[dist_key]['label']} ({method.upper()})")
    ax.fill_between(ci["T"], ci["lower"], ci["upper"], color="crimson", alpha=0.2,
                     label=f"{ci_pct:g}% bootstrap CI")

    ax.set_xscale("log")
    ax.set_xlabel("Return period T (years)")
    ax.set_ylabel("Flood magnitude")
    ax.set_title(f"{DISTRIBUTIONS[dist_key]['label']} fit with bootstrap CI — {ffa.station_id}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    return ax


# Reference (t3, t4) curves for common distributions, used as an overlay on
# the L-moment ratio diagram. Values from Hosking & Wallis (1997).
def _lmr_reference_curves(t3):
    t3 = np.asarray(t3)
    curves = {}
    # Generalized logistic
    curves["GLO"] = t3 ** 2
    # Generalized Pareto
    curves["GPA"] = t3 * (1 - t3) / (1 + t3) if np.all(t3 != -1) else None
    return curves


def moment_ratio_diagram(ffa, ax=None):
    """
    L-moment ratio diagram (t3 = L-skewness on x, t4 = L-kurtosis on y):
    plots the sample point together with theoretical curves for a few
    two-parameter-shape distribution families, to visually suggest which
    family best matches the data's shape.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    stats = summarize(ffa.data)
    t3, t4 = stats["t3"], stats["t4"]

    t3_grid = np.linspace(-0.9, 0.9, 200)
    # Generalized logistic: t4 = (1 + 5 t3^2) / 6
    ax.plot(t3_grid, (1 + 5 * t3_grid ** 2) / 6, label="Gen. Logistic (GLO)")
    # Generalized Pareto: t4 = t3(1+5t3)/(5+t3)
    ax.plot(t3_grid, t3_grid * (1 + 5 * t3_grid) / (5 + t3_grid), label="Gen. Pareto (GPA)")
    # Gumbel (point): t3=0.1699, t4=0.1504
    ax.scatter([0.1699], [0.1504], marker="s", color="orange", label="Gumbel (point)")
    # Normal (point): t3=0, t4=0.1226
    ax.scatter([0.0], [0.1226], marker="^", color="green", label="Normal (point)")

    ax.scatter([t3], [t4], marker="*", s=250, color="crimson", zorder=6,
               label=f"{ffa.station_id} sample")

    ax.set_xlabel("L-skewness (t3)")
    ax.set_ylabel("L-kurtosis (t4)")
    ax.set_title("L-moment ratio diagram")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    return ax


def data_histogram(ffa, dist_key=None, method=None, bins=None, ax=None):
    """
    Histogram of the raw annual-maximum input series. If dist_key/method are
    given (or a fit exists to use as default), overlays that distribution's
    fitted PDF for a quick visual check of fit quality.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    x = ffa.data
    if bins is None:
        # Sturges' rule as a reasonable default for typical flood record lengths
        bins = max(5, int(np.ceil(np.log2(x.size) + 1)))

    ax.hist(x, bins=bins, density=True, color="steelblue", edgecolor="white",
            alpha=0.75, label="Observed data")

    # Overlay a fitted PDF if available (explicit dist_key/method, or the
    # first one available in ffa.fits as a fallback)
    r = None
    if dist_key is not None and method is not None:
        r = ffa.fits.get((dist_key, method))
    elif ffa.fits:
        (dist_key, method), r = next(iter(ffa.fits.items()))

    if r is not None:
        x_grid = np.linspace(x.min() * 0.9, x.max() * 1.1, 300)
        x_grid = x_grid[x_grid > 0]  # guard for log-based distributions
        ax.plot(x_grid, r.pdf(x_grid), color="crimson", lw=2,
                label=f"{DISTRIBUTIONS[dist_key]['label']} ({method.upper()}) fit")

    ax.set_xlabel("Flood magnitude")
    ax.set_ylabel("Density")
    ax.set_title(f"Input data histogram — {ffa.station_id}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    return ax


# --------------------------------------------------------------------- #
# Script-friendly save_* wrappers: build the figure, save a PNG, close it.
# Used by run_analysis.py so batch runs over many cases don't leak figures.
# --------------------------------------------------------------------- #

def save_probability_plot(ffa, path, dist_keys=None, max_T=10000, dpi=150):
    fig, ax = plt.subplots(figsize=(7, 5))
    probability_plot(ffa, dist_keys=dist_keys, ax=ax, max_T=max_T)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def save_data_histogram(ffa, path, dist_key=None, method=None, bins=None, dpi=150):
    fig, ax = plt.subplots(figsize=(7, 5))
    data_histogram(ffa, dist_key=dist_key, method=method, bins=bins, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def save_quantile_ci_plot(ffa, dist_key, method, path, n_boot=500, alpha=0.05, dpi=150, random_state=0):
    fig, ax = plt.subplots(figsize=(7, 5))
    quantile_plot_with_ci(ffa, dist_key, method, n_boot=n_boot, alpha=alpha, ax=ax, random_state=random_state)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def save_moment_ratio_diagram(ffa, path, dpi=150):
    fig, ax = plt.subplots(figsize=(6, 6))
    moment_ratio_diagram(ffa, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path