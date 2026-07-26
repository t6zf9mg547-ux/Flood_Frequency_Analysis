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
    ax.set_ylabel(ffa.axis_label())
    ax.set_title(f"{ffa.short_name} frequency curve — {ffa.station_id}")
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
    ax.set_ylabel(ffa.axis_label())
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

    ax.set_xlabel(ffa.axis_label())
    ax.set_ylabel("Density")
    ax.set_title(f"Input data histogram — {ffa.station_id}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    return ax


def data_quality_plot(ffa, ax=None):
    """
    Time-series plot of the annual-maximum series: raw values against
    year (or index, if no year column), with the Sen's-slope trend line
    and any Grubbs-flagged outliers highlighted, so the Mann-Kendall/
    Grubbs test results in the summary can be sanity-checked visually.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4.5))

    x = ffa.data
    t = ffa.years if ffa.years is not None else np.arange(1, x.size + 1)
    t = np.asarray(t, dtype=float)

    dqr = ffa.data_quality()
    mk = dqr["mann_kendall"]
    sen = dqr["sens_slope"]
    gr = dqr["grubbs"]

    ax.plot(t, x, "o-", color="steelblue", ms=4, lw=1, label="Annual maximum")

    trend_line = sen["intercept"] + sen["slope"] * t
    trend_color = "crimson" if mk["significant"] else "gray"
    trend_style = "-" if mk["significant"] else "--"
    ax.plot(t, trend_line, trend_style, color=trend_color, lw=1.8,
            label=f"Sen's slope trend ({mk['trend']})")

    # highlight Grubbs-flagged outliers, if any
    if gr["high_outlier_flagged"]:
        idx = np.argmax(x)
        ax.scatter(t[idx], x[idx], s=140, facecolors="none", edgecolors="red",
                   linewidths=2, zorder=6, label="Flagged high outlier (Grubbs)")
    if gr["low_outlier_flagged"]:
        idx = np.argmin(x)
        ax.scatter(t[idx], x[idx], s=140, facecolors="none", edgecolors="orange",
                   linewidths=2, zorder=6, label="Flagged low outlier (Grubbs)")

    xlabel = "Year" if ffa.years is not None else "Record index"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ffa.axis_label())
    ax.set_title(f"Data quality: time series & trend — {ffa.station_id}\n"
                 f"Mann-Kendall: {mk['trend']} (p={mk['p_value']:.3f})")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    return ax


# --------------------------------------------------------------------- #
# Script-friendly save_* wrappers: build the figure, save a PNG, close it.
# Used by run_analysis.py so batch runs over many cases don't leak figures.
# --------------------------------------------------------------------- #

def dashboard(ffa, best_key, best_method, n_boot=500, alpha=0.05, fig=None):
    """
    One-page overview combining all 5 individual plots plus a text panel
    with the key numbers, for a quick at-a-glance review instead of
    flipping between 5 separate PNGs.
    """
    if fig is None:
        fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

    ax_ts = fig.add_subplot(gs[0, 0])
    data_quality_plot(ffa, ax=ax_ts)
    ax_ts.set_title(ax_ts.get_title(), fontsize=9)

    ax_hist = fig.add_subplot(gs[0, 1])
    data_histogram(ffa, dist_key=best_key, method=best_method, ax=ax_hist)
    ax_hist.set_title(ax_hist.get_title(), fontsize=9)

    ax_lmr = fig.add_subplot(gs[0, 2])
    moment_ratio_diagram(ffa, ax=ax_lmr)
    ax_lmr.set_title(ax_lmr.get_title(), fontsize=9)

    ax_prob = fig.add_subplot(gs[1, 0])
    probability_plot(ffa, ax=ax_prob)
    ax_prob.set_title(ax_prob.get_title(), fontsize=9)
    ax_prob.legend(fontsize=6)

    ax_ci = fig.add_subplot(gs[1, 1])
    quantile_plot_with_ci(ffa, best_key, best_method, n_boot=n_boot, alpha=alpha, ax=ax_ci)
    ax_ci.set_title(ax_ci.get_title(), fontsize=9)

    # -- Text summary panel --
    ax_txt = fig.add_subplot(gs[1, 2])
    ax_txt.axis("off")
    stats = ffa.descriptive_stats()
    table = ffa.goodness_of_fit_table()
    best_row = table[(table["key"] == best_key) & (table["method"] == best_method)].iloc[0]
    weights = ffa.akaike_weights()
    top_weight = float(weights["akaike_weight"].iloc[0])
    dqr = ffa.data_quality()
    mk = dqr["mann_kendall"]
    gr = dqr["grubbs"]
    outlier_txt = ("none flagged" if not (gr["high_outlier_flagged"] or gr["low_outlier_flagged"])
                   else "SEE WARNING in summary.txt")

    lines = [
        f"Station: {ffa.station_id}",
        f"n = {ffa.n} years",
        "",
        f"Mean = {stats['mean']:.1f}   CV = {stats['CV']:.2f}   CS = {stats['CS']:.2f}",
        "",
        f"Trend: {mk['trend']}",
        f"Outliers: {outlier_txt}",
        "",
        f"Recommended: {DISTRIBUTIONS[best_key]['label']}",
        f"  ({best_method.upper()})",
        f"AIC = {best_row['AIC']:.1f}   KS p = {best_row['KS_pvalue']:.3f}",
        f"Akaike weight = {top_weight:.0%}",
        "",
        "See summary.txt for full",
        "quality assessments and",
        "the extrapolation warning.",
    ]
    ax_txt.text(0.02, 0.98, "\n".join(lines), transform=ax_txt.transAxes,
                fontsize=9, va="top", family="monospace")

    fig.suptitle(f"{ffa.short_name} Frequency Analysis Dashboard — {ffa.station_id}", fontsize=14, fontweight="bold")
    return fig


def save_dashboard(ffa, path, best_key, best_method, n_boot=500, alpha=0.05, dpi=150):
    fig = plt.figure(figsize=(16, 10))
    dashboard(ffa, best_key, best_method, n_boot=n_boot, alpha=alpha, fig=fig)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _text_page(fig, text, fontsize=8):
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=fontsize,
            va="top", family="monospace")


def save_pdf_report(ffa, path, best_key, best_method, recommendation_text,
                     n_boot=500, alpha=0.05, lines_per_page=62, dpi=150):
    """
    Single PDF bundling the text summary (paginated) + the dashboard +
    each individual full-size plot -- the one file to actually hand to a
    colleague or attach to an email, rather than 5 separate PNGs and a
    .txt file.

    recommendation_text: pass the already-computed string from
    ffa.generate_recommendation() -- accepted as a parameter rather than
    recomputed here, since it runs a bootstrap internally and the caller
    (run_analysis.py) has typically already paid that cost once.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        # -- Text summary, paginated -- #
        text_lines = recommendation_text.split("\n")
        for i in range(0, len(text_lines), lines_per_page):
            chunk = "\n".join(text_lines[i:i + lines_per_page])
            fig = plt.figure(figsize=(8.5, 11))
            _text_page(fig, chunk)
            pdf.savefig(fig, dpi=dpi)
            plt.close(fig)

        # -- Dashboard overview -- #
        fig = plt.figure(figsize=(16, 10))
        dashboard(ffa, best_key, best_method, n_boot=n_boot, alpha=alpha, fig=fig)
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)

        # -- Individual full-size plots -- #
        fig, ax = plt.subplots(figsize=(10, 6.5))
        data_quality_plot(ffa, ax=ax)
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 6.5))
        data_histogram(ffa, dist_key=best_key, method=best_method, ax=ax)
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 6.5))
        probability_plot(ffa, ax=ax)
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 8))
        moment_ratio_diagram(ffa, ax=ax)
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 6.5))
        quantile_plot_with_ci(ffa, best_key, best_method, n_boot=n_boot, alpha=alpha, ax=ax)
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)

        d = pdf.infodict()
        d["Title"] = f"{ffa.short_name} Frequency Analysis — {ffa.station_id}"
        d["Author"] = "floodfreq"

    return path


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


def save_data_quality_plot(ffa, path, dpi=150):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    data_quality_plot(ffa, ax=ax)
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
