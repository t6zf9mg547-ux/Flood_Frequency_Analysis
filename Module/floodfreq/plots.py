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


def _label_points_staggered(ax, labels, x, y, fontsize=7, min_px_gap=13):
    """
    Label a scatter of points without the labels overlapping when several
    points sit close together (e.g. stations that pool well cluster
    tightly on the L-moment ratio diagram, which is exactly the case
    where the *default* one-annotation-per-point placement collides
    worst). No external dependency (no adjustText): points are ordered by
    y-value and given a monotonically increasing vertical text offset
    (a small "ladder"), alternating left/right of the point, connected
    back to their actual location with a thin leader line -- enough to
    resolve the common case of a tight cluster without the complexity of
    true collision detection.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(y)
    n = len(order)
    for rank, idx in enumerate(order):
        dy = (rank - (n - 1) / 2.0) * min_px_gap
        left_side = rank % 2 == 1
        dx = -28 if left_side else 14
        ha = "right" if left_side else "left"
        ax.annotate(labels[idx], (x[idx], y[idx]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=fontsize, ha=ha, va="center",
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6, alpha=0.7,
                                     shrinkA=0, shrinkB=2))


def regional_moment_ratio_diagram(stations, region_name="", discordancy_df=None, ax=None):
    """
    Multi-station L-moment ratio diagram: extends `moment_ratio_diagram`
    to plot every station in a pooling group as its own (t3, t4) point,
    plus the record-length-weighted regional average, against the same
    theoretical GLO/GPA reference curves (closed-form) and, additionally,
    GEV/GNO/PE3 reference curves (computed numerically via
    `regional.theoretical_tau4`, since those three don't have a simple
    closed form). Used as a visual homogeneity/distribution-selection
    check: a tight, non-scattered cluster of station points supports
    pooling; the regional average's position relative to the reference
    curves suggests which family fits best (cross-checked numerically by
    the Z-statistic in `regional.zstatistics`).

    `stations`: list of `regional.StationLMoments`.
    `discordancy_df`: optional output of `regional.discordancy(stations)`
        -- if given, discordant stations are marked with a red outline.
    """
    from . import regional as _regional  # local import: avoid a hard
    # plots.py -> regional.py -> ... dependency for callers that only use
    # the single-station plots and don't have lmoments3's Kappa etc. handy.

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    t3_grid = np.linspace(-0.9, 0.9, 200)
    ax.plot(t3_grid, (1 + 5 * t3_grid ** 2) / 6, label="Gen. Logistic (GLO)", lw=1.5)
    ax.plot(t3_grid, t3_grid * (1 + 5 * t3_grid) / (5 + t3_grid), label="Gen. Pareto (GPA)", lw=1.5)

    for fam, style in (("gev", "--"), ("gno", "-."), ("pe3", ":")):
        t3_sub = np.linspace(-0.6, 0.6, 25)
        t4_curve = [_regional.theoretical_tau4(t3, fam) for t3 in t3_sub]
        ax.plot(t3_sub, t4_curve, style, lw=1.5, label=f"{fam.upper()} (numeric)")

    ax.scatter([0.1699], [0.1504], marker="s", color="orange", label="Gumbel (point)")
    ax.scatter([0.0], [0.1226], marker="^", color="green", label="Normal (point)")

    t3_vals = np.array([s.t3 for s in stations])
    t4_vals = np.array([s.t4 for s in stations])
    discordant = np.zeros(len(stations), dtype=bool)
    if discordancy_df is not None:
        d_by_name = dict(zip(discordancy_df["station"], discordancy_df["discordant"]))
        discordant = np.array([d_by_name.get(s.name, False) for s in stations])

    ax.scatter(t3_vals[~discordant], t4_vals[~discordant], marker="o", s=60,
               color="steelblue", zorder=5, label="Stations")
    if discordant.any():
        ax.scatter(t3_vals[discordant], t4_vals[discordant], marker="o", s=90,
                   facecolors="none", edgecolors="crimson", linewidths=2, zorder=6,
                   label="Discordant station(s)")
    _label_points_staggered(ax, [s.name for s in stations], t3_vals, t4_vals, fontsize=7)

    reg = _regional.regional_average_ratios(stations)
    ax.scatter([reg["t3_R"]], [reg["t4_R"]], marker="*", s=300, color="crimson",
               zorder=7, label="Regional average")

    ax.set_xlabel("L-skewness (t3)")
    ax.set_ylabel("L-kurtosis (t4)")
    title = "Regional L-moment ratio diagram"
    if region_name:
        title += f" — {region_name}"
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="upper left")
    return ax


def save_regional_moment_ratio_diagram(stations, path, region_name="", discordancy_df=None, dpi=150):
    fig, ax = plt.subplots(figsize=(7, 7))
    regional_moment_ratio_diagram(stations, region_name=region_name,
                                   discordancy_df=discordancy_df, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


# --------------------------------------------------------------------- #
# Regional (pooled) plot suite -- the index-flood-method equivalents of
# the single-station plots above (probability_plot, quantile_plot_with_ci,
# data_quality_plot, dashboard).
# --------------------------------------------------------------------- #
from .plotting_positions import empirical_frequency, return_period as _return_period_from_F


def regional_growth_curve_plot(result, ax=None, log_x=True, max_T=10000, station_formula="weibull",
                                 show_pooled_rank=True):
    """
    The regional analogue of `probability_plot`: the fitted, dimensionless
    growth curve g(F) (station design flood / that station's own index
    flood), plotted against return period T, together with every pooled
    station's own data rescaled the same way (x_i / station i's mean) so
    all stations land on one common, unit-mean axis -- the classic
    "Dalrymple plot" used to eyeball how well one growth curve serves the
    whole group.

    Each station's points use that station's OWN rank/plotting position
    (n_i years -> ranks 1..n_i), which is what the index-flood method
    actually fits to: the regional L-moments are a record-length-weighted
    AVERAGE of each station's own L-moment ratios, not the L-moments of
    one flat concatenated series.

    `show_pooled_rank` adds a secondary, purely diagnostic overlay: every
    station's dimensionless growth factor pooled into one array and
    RE-ranked together (as if it were a single N = sum(n_i)-year record),
    plotted in light grey behind the per-station points. This is NOT what
    the growth curve is fitted to -- concatenating stations before ranking
    would let whichever station has the most years dominate the shape
    estimate -- but it's a useful sanity check: if the grey pooled-rank
    cloud and the colored per-station points trace visibly different
    shapes, that's a sign the stations don't reduce to one common
    dimensionless distribution as cleanly as the per-station view alone
    might suggest.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    if show_pooled_rank:
        pooled = np.sort(np.concatenate([
            np.asarray(result.station_data[s.name], dtype=float) / s.mean
            for s in result.stations
        ]))
        F_pooled = empirical_frequency(pooled.size, formula=station_formula, ascending_rank=True)
        T_pooled = _return_period_from_F(F_pooled)
        ax.scatter(T_pooled, pooled, s=10, color="0.65", alpha=0.5, marker="x", zorder=2,
                   label=f"Pooled rank (all {pooled.size} station-years combined, diagnostic only)")

    cmap = plt.get_cmap("tab10")
    for i, s in enumerate(result.stations):
        x = np.sort(np.asarray(result.station_data[s.name], dtype=float))
        F = empirical_frequency(x.size, formula=station_formula, ascending_rank=True)
        T = _return_period_from_F(F)
        ax.scatter(T, x / s.mean, s=18, color=cmap(i % 10), alpha=0.75,
                   label=s.name, zorder=4)

    T_smooth = np.geomspace(1.01, max_T, 300)
    F_smooth = 1 - 1 / T_smooth
    ax.plot(T_smooth, result.growth_curve.ppf(F_smooth), color="black", lw=2.2,
            zorder=5, label=f"Regional growth curve ({result.growth_curve.label})")

    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel("Return period T (years)")
    ax.set_ylabel("Growth factor (value / index flood)", fontsize=9)
    ax.set_title(f"Regional growth curve — {result.region_name}")
    ax.grid(True, which="both", alpha=0.3)
    # "lower right" is reliably empty: at high T the curve rises well above
    # the data cloud, leaving the low-growth-factor / high-T corner clear --
    # unlike "best" (matplotlib's default), which tends to land the legend
    # box on top of the y-axis label/ticks in the upper-left.
    ax.legend(fontsize=7, ncol=1, loc="lower right")
    return ax


def save_regional_growth_curve_plot(result, path, max_T=10000, dpi=150, show_pooled_rank=True):
    fig, ax = plt.subplots(figsize=(7, 5))
    regional_growth_curve_plot(result, ax=ax, max_T=max_T, show_pooled_rank=show_pooled_rank)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def regional_pooled_vs_stations_plot(result, log_x=True, max_T=10000, station_formula="weibull",
                                       fig=None):
    """
    Two panels, side by side, at equal visual weight, so the difference
    between the two ways of looking at pooled data is unmistakable rather
    than one being a faint overlay on the other:

    LEFT  -- per-station view: each station's own data at its own
             plotting-position ranks (n_i years -> ranks 1..n_i), colored
             by station. This is what the growth curve is actually fit to
             (a record-length-weighted average of each station's own
             L-moment ratios).
    RIGHT -- "as one station" view: every station's dimensionless growth
             factor pooled into a single array and re-ranked together, as
             if it were one combined N = sum(n_i)-year record. All points
             are the same color (no station identity) because that's the
             point of this view -- once concatenated, station identity is
             gone. This is NOT what the curve is fit to; it's shown here
             purely so you can compare it directly against the left panel.

    Both panels carry the same fitted regional growth curve for reference.
    Notice the right panel reaches further out along the T axis than any
    single station's own record could (up to T ~ N years, N = total
    pooled station-years) -- that reach is an artifact of concatenation
    inflating the apparent sample size, not extra information the data
    actually contain about rare events; it's one reason the index-flood
    method fits the weighted-average-of-ratios way (left panel) instead.
    """
    if fig is None:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    else:
        axes = fig.subplots(1, 2)
    ax_left, ax_right = axes

    T_smooth = np.geomspace(1.01, max_T, 300)
    F_smooth = 1 - 1 / T_smooth
    curve = result.growth_curve.ppf(F_smooth)

    # -- left: per-station --
    cmap = plt.get_cmap("tab10")
    for i, s in enumerate(result.stations):
        x = np.sort(np.asarray(result.station_data[s.name], dtype=float))
        F = empirical_frequency(x.size, formula=station_formula, ascending_rank=True)
        T = _return_period_from_F(F)
        ax_left.scatter(T, x / s.mean, s=20, color=cmap(i % 10), alpha=0.8, label=s.name, zorder=4)
    ax_left.plot(T_smooth, curve, color="black", lw=2.2, zorder=5, label="Regional growth curve")
    ax_left.set_title("Per-station (what the curve is fit to)")
    ax_left.legend(fontsize=7, ncol=2)

    # -- right: pooled as one station --
    pooled = np.sort(np.concatenate([
        np.asarray(result.station_data[s.name], dtype=float) / s.mean
        for s in result.stations
    ]))
    F_pooled = empirical_frequency(pooled.size, formula=station_formula, ascending_rank=True)
    T_pooled = _return_period_from_F(F_pooled)
    ax_right.scatter(T_pooled, pooled, s=20, color="0.25", alpha=0.6, zorder=4,
                     label=f"All {pooled.size} station-years, ranked as one series")
    ax_right.plot(T_smooth, curve, color="crimson", lw=2.2, zorder=5, label="Same regional growth curve")
    ax_right.set_title('"As one station" (diagnostic only -- NOT what is fitted)')
    ax_right.legend(fontsize=7)

    for ax in (ax_left, ax_right):
        if log_x:
            ax.set_xscale("log")
        ax.set_xlabel("Return period T (years)")
        ax.set_ylabel("Growth factor (value / index flood)")
        ax.grid(True, which="both", alpha=0.3)
    ax_right.set_ylim(ax_left.get_ylim())

    fig.suptitle(f"Per-station vs. pooled-as-one-station — {result.region_name}",
                 fontsize=13, fontweight="bold")
    return fig


def save_regional_pooled_vs_stations_plot(result, path, max_T=10000, dpi=150):
    fig = plt.figure(figsize=(13, 5.5))
    regional_pooled_vs_stations_plot(result, max_T=max_T, fig=fig)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def station_design_flood_plot(result, station_name, ax=None, log_x=True, max_T=10000,
                                station_formula="weibull"):
    """
    The regional analogue of `quantile_plot_with_ci` / `probability_plot`
    for ONE station: that station's own observed data (at its own
    plotting-position return periods) against the design-flood curve
    implied by the regional growth curve scaled by this station's index
    flood -- i.e. what an engineer would actually hand to a client for
    this specific site.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    s = next((st for st in result.stations if st.name == station_name), None)
    if s is None:
        raise ValueError(f"Unknown station '{station_name}'. Available: "
                          f"{[st.name for st in result.stations]}")
    x = np.sort(np.asarray(result.station_data[station_name], dtype=float))
    F = empirical_frequency(x.size, formula=station_formula, ascending_rank=True)
    T = _return_period_from_F(F)
    ax.scatter(T, x, s=25, color="black", zorder=5, label="Observed data")

    T_smooth = np.geomspace(max(1.01, T.min()), max(T.max(), max_T), 300)
    ax.plot(T_smooth, result.growth_curve.station_quantile(T_smooth, s.mean),
            color="crimson", lw=2, zorder=4,
            label=f"Regional design flood ({result.growth_curve.label})")

    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel("Return period T (years)")
    ax.set_ylabel("Flood magnitude")
    ax.set_title(f"Regional design flood — {station_name} ({result.region_name})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    return ax


def save_station_design_flood_plot(result, station_name, path, max_T=10000, dpi=150):
    fig, ax = plt.subplots(figsize=(7, 5))
    station_design_flood_plot(result, station_name, ax=ax, max_T=max_T)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def save_regional_station_series_plot(result, path, ncols=3, dpi=150):
    """
    Small multiples: one panel per station, its raw annual-maximum series
    plotted against calendar year when `result.station_years` has years
    for that station (as returned by `io_utils.load_region_stations()`),
    or against plain observation order (1..n) otherwise -- a quick visual
    screen for anything unusual (a station with a wildly different range,
    an obvious jump, a short/noisy record) before trusting the pooling.
    """
    n = len(result.stations)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    for i, s in enumerate(result.stations):
        ax = axes.flatten()[i]
        x = np.asarray(result.station_data[s.name], dtype=float)
        years = result.station_years.get(s.name) if result.station_years else None
        t = np.asarray(years, dtype=float) if years is not None else np.arange(1, x.size + 1)
        ax.plot(t, x, "o-", ms=3, lw=1, color="steelblue")
        ax.axhline(s.mean, color="gray", ls="--", lw=1)
        ax.set_title(s.name, fontsize=9)
        ax.set_xlabel("Year" if years is not None else "Observation #", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
    for j in range(len(result.stations), nrows * ncols):
        axes.flatten()[j].axis("off")
    fig.suptitle(f"Station annual-maximum series — {result.region_name}", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def regional_discordancy_plot(result, ax=None):
    """Bar chart of each station's discordancy measure D_i against the
    group's critical value -- a quick visual companion to discordancy.csv."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))

    df = result.discordancy_df.sort_values("D_i", ascending=False)
    colors = ["crimson" if d else "steelblue" for d in df["discordant"]]
    ax.bar(df["station"], df["D_i"], color=colors)
    crit = float(df["D_critical"].iloc[0])
    ax.axhline(crit, color="black", ls="--", lw=1.3, label=f"Critical D_i = {crit:.3f}")
    ax.set_ylabel("Discordancy D_i")
    ax.set_title(f"Station discordancy — {result.region_name}")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    return ax


def save_regional_discordancy_plot(result, path, dpi=150):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    regional_discordancy_plot(result, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def regional_dashboard(result, fig=None):
    """
    One-page overview for a regional analysis, the pooled-data counterpart
    of the single-station `dashboard()`: growth curve, L-moment ratio
    diagram, discordancy bar chart, and a text summary panel.
    """
    if fig is None:
        fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    ax_growth = fig.add_subplot(gs[0, 0])
    regional_growth_curve_plot(result, ax=ax_growth)
    ax_growth.set_title(ax_growth.get_title(), fontsize=9)
    ax_growth.legend(fontsize=6, ncol=1, loc="lower right")

    ax_lmr = fig.add_subplot(gs[0, 1])
    regional_moment_ratio_diagram(result.stations, region_name=result.region_name,
                                   discordancy_df=result.discordancy_df, ax=ax_lmr)
    ax_lmr.set_title(ax_lmr.get_title(), fontsize=9)
    ax_lmr.legend(fontsize=6, loc="upper left")

    ax_disc = fig.add_subplot(gs[1, 0])
    regional_discordancy_plot(result, ax=ax_disc)
    ax_disc.set_title(ax_disc.get_title(), fontsize=9)

    ax_txt = fig.add_subplot(gs[1, 1])
    ax_txt.axis("off")
    het = result.heterogeneity_result
    n_disc = int(result.discordancy_df["discordant"].sum())
    lines = [
        f"Region: {result.region_name}",
        f"Stations: {len(result.stations)}   "
        f"Total station-years: {int(result.stations_df['n'].sum())}",
        "",
        f"Discordant stations: {n_disc if n_disc else 'none'}",
        f"Heterogeneity H1 = {het.H1:.2f}  ({het.interpretation})",
        "",
        f"Regional distribution: {result.growth_curve.label}",
        f"  t_R (L-CV) = {result.growth_curve.t_R:.3f}",
        f"  t3_R (L-skew) = {result.growth_curve.t3_R:.3f}",
        "",
        "Growth factors:",
    ] + [f"  T={int(T):>6,d}y: {float(result.growth_curve.quantile(T)):.2f}"
         for T in (10, 100, 1000, 10000)] + [
        "",
        "See summary.txt for the full",
        "station-by-station walkthrough.",
    ]
    ax_txt.text(0.02, 0.98, "\n".join(lines), transform=ax_txt.transAxes,
                fontsize=9, va="top", family="monospace")

    fig.suptitle(f"Regional Flood Frequency Analysis Dashboard — {result.region_name}",
                 fontsize=14, fontweight="bold")
    return fig


def save_regional_dashboard(result, path, dpi=150):
    fig = plt.figure(figsize=(15, 10))
    regional_dashboard(result, fig=fig)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path
