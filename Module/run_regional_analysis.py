#!/usr/bin/env python3
"""
Run a regional (pooled, index-flood) flood frequency analysis for one group
of hydrologically-similar stations.

Layout expected (this script lives in <project_root>/Module/):

    <project_root>/
        Data/Regional/<RegionName>/*.csv   <- one file per station
        Module/run_regional_analysis.py    <- this file
        Module/floodfreq/                  <- the package (floodfreq.regional)
        Output/Regional/<RegionName>/      <- CSV results written here
        Plot/Regional/<RegionName>/        <- PNG plots written here

Method: Hosking & Wallis (1997) L-moment-based index-flood procedure
(Dalrymple, 1960). Each station contributes its own L-moment ratios; the
group is screened for discordant stations (D_i) and overall heterogeneity
(H1); a regional growth curve is fit to the pooled, record-length-weighted
L-moments of whichever candidate family best matches the regional-average
L-kurtosis (Z-statistic); station design floods are then that curve's
quantile times the station's own index flood (its sample mean).

Usage:
    python run_regional_analysis.py Piedmont
    python run_regional_analysis.py Piedmont --n-sim 1000 --seed 42
    python run_regional_analysis.py Piedmont --family gev
    python run_regional_analysis.py Piedmont --value-col Q --year-col year
"""
from __future__ import annotations
import argparse
import sys
import datetime
import platform
from pathlib import Path

import pandas as pd

import floodfreq
from floodfreq.io_utils import resolve_region, load_region_stations, project_root_from
from floodfreq.regional import (
    run_regional_analysis,
    CANDIDATE_FAMILIES,
    discordancy_critical_value,
    regional_quantile_ci,
)
from floodfreq.plots import (
    save_regional_moment_ratio_diagram,
    save_regional_growth_curve_plot,
    save_regional_pooled_vs_stations_plot,
    save_station_design_flood_plot,
    save_regional_station_series_plot,
    save_regional_discordancy_plot,
    save_regional_dashboard,
)

DEFAULT_RETURN_PERIODS = [2, 5, 10, 20, 25, 50, 100, 200, 500, 1000]


def build_provenance_header() -> str:
    """Same idea as run_analysis.py's provenance header: stamp when/how/
    with-what-version this run was produced, so a result found later is
    self-documenting."""
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    command = " ".join([Path(sys.argv[0]).name] + sys.argv[1:])
    lines = [
        "=" * 70,
        "RUN PROVENANCE",
        "=" * 70,
        f"Generated:     {now}",
        f"Tool version:  floodfreq {floodfreq.__version__}",
        f"Command:       python {command}",
        f"Python:        {platform.python_version()}",
        "=" * 70,
        "",
    ]
    return "\n".join(lines)


def _csv_hint(path: Path) -> str:
    """Row count + last-modified time, matching run_analysis.py's picker."""
    try:
        with open(path) as f:
            n_rows = max(sum(1 for _ in f) - 1, 0)
    except Exception:
        n_rows = "?"
    mtime = datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    return f"{n_rows} rows, modified {mtime}"


def prompt_region_name(project_root: Path) -> str:
    """
    List every subfolder of Data/Regional/ and let the user pick one
    interactively. Used when region_name isn't passed on the command line.
    """
    regional_dir = project_root / "Data" / "Regional"
    regions = sorted(d for d in regional_dir.glob("*") if d.is_dir()) if regional_dir.is_dir() else []
    if not regions:
        sys.exit(f"ERROR: no region folders found in {regional_dir}. "
                  f"Create one (e.g. {regional_dir}/<RegionName>/) and add "
                  f"at least 4 station CSVs to it, then try again.")

    if len(regions) == 1:
        only = regions[0]
        n_csv = len(list(only.glob("*.csv")))
        raw = input(f"Found one region: {only.name} ({n_csv} station CSV(s)) — use it? [Y/n]: ").strip().lower()
        if raw in ("", "y", "yes"):
            return only.name

    print(f"\nAvailable regions in {regional_dir}:")
    for i, d in enumerate(regions, start=1):
        n_csv = len(list(d.glob("*.csv")))
        print(f"  {i}. {d.name}  ({n_csv} station CSV(s))")

    while True:
        raw = input(f"Choose a region [1-{len(regions)}] or type a region name: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(regions):
            return regions[int(raw) - 1].name
        if (regional_dir / raw).is_dir():
            return raw
        print(f"  Not a valid choice. Enter a number 1-{len(regions)}, or an exact region name.")


def build_summary_text(result, provenance: str, ci_df=None, confidence_level: float = 95.0) -> str:
    """Plain-language regional-analysis summary, matching the tone/structure
    of run_analysis.py's summary.txt (Good/Caution/Warning-style framing,
    but keyed to the pooling-specific diagnostics)."""
    lines = [provenance]
    lines.append(f"REGIONAL (INDEX-FLOOD) FLOOD FREQUENCY ANALYSIS — {result.region_name}")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Stations pooled: {len(result.stations)}")
    lines.append(f"Total station-years: {int(result.stations_df['n'].sum())}")
    lines.append("")

    # -- Per-station data quality (run before pooling) --
    lines.append("-" * 70)
    lines.append("STEP 1: PER-STATION DATA QUALITY (checked before pooling)")
    lines.append("-" * 70)
    dq = result.data_quality_df
    trend_flagged = dq[dq["mann_kendall_significant"]]
    outlier_flagged = dq[dq["grubbs_high_outlier_flagged"] | dq["grubbs_low_outlier_flagged"]]
    if trend_flagged.empty and outlier_flagged.empty:
        lines.append("Good: no station shows a significant Mann-Kendall trend or a Grubbs-flagged")
        lines.append("outlier -- nothing here suggests non-stationarity or a data error before pooling.")
    else:
        if not trend_flagged.empty:
            lines.append("Caution: the following station(s) show a significant trend, which questions")
            lines.append("the stationarity assumption behind pooling them as-is:")
            for _, row in trend_flagged.iterrows():
                lines.append(f"  - {row['station']}: {row['mann_kendall_trend']} "
                              f"(p={row['mann_kendall_p_value']:.3f})")
        if not outlier_flagged.empty:
            lines.append("Caution: the following station(s) have a Grubbs-flagged outlier, worth a")
            lines.append("closer look (data error? a genuinely extreme event?):")
            for _, row in outlier_flagged.iterrows():
                which = []
                if row["grubbs_high_outlier_flagged"]:
                    which.append(f"high={row['grubbs_high_outlier_value']:.1f}")
                if row["grubbs_low_outlier_flagged"]:
                    which.append(f"low={row['grubbs_low_outlier_value']:.1f}")
                lines.append(f"  - {row['station']}: {', '.join(which)}")
    lines.append("")
    lines.append(dq.to_string(index=False))
    lines.append("")

    # -- Discordancy --
    lines.append("-" * 70)
    lines.append("STEP 2: DISCORDANCY (are any stations statistically odd-ones-out?)")
    lines.append("-" * 70)
    crit = discordancy_critical_value(len(result.stations))
    lines.append(f"Critical D_i for a {len(result.stations)}-station region: {crit:.3f}")
    disc = result.discordancy_df
    flagged = disc[disc["discordant"]]
    if flagged.empty:
        lines.append("Good: no station exceeds the critical value — none are flagged discordant.")
    else:
        lines.append("Warning: the following station(s) are flagged discordant and are worth a closer")
        lines.append("look before pooling (data error? a genuinely different flood regime?):")
        for _, row in flagged.iterrows():
            lines.append(f"  - {row['station']}: D_i = {row['D_i']:.3f} (> {crit:.3f})")
    lines.append("")
    lines.append(disc.to_string(index=False))
    lines.append("")

    # -- Heterogeneity --
    lines.append("-" * 70)
    lines.append("STEP 3: HETEROGENEITY (is the group homogeneous enough to pool?)")
    lines.append("-" * 70)
    het = result.heterogeneity_result
    lines.append(f"H1 = {het.H1:.3f}  -> {het.interpretation}")
    lines.append(f"  (observed weighted L-CV dispersion V1 = {het.V1_obs:.4f}; "
                  f"{het.n_sim} simulated homogeneous regions gave "
                  f"V1 = {het.V1_sim_mean:.4f} +/- {het.V1_sim_std:.4f}, "
                  f"simulated via the '{het.simulation_family}' distribution)")
    if het.H1 < 1:
        lines.append("Good: H1 < 1 — conventionally read as an acceptably homogeneous region.")
    elif het.H1 < 2:
        lines.append("Caution: 1 <= H1 < 2 — 'possibly heterogeneous'; growth-curve quantiles for")
        lines.append("individual stations are less reliable than for a clearly homogeneous region.")
    else:
        lines.append("Warning: H1 >= 2 — 'definitely heterogeneous'. Pooling these stations under a")
        lines.append("single growth curve is not recommended as-is; consider removing discordant")
        lines.append("station(s) above, splitting the region, or reviewing station selection.")
    lines.append("")

    # -- Distribution selection --
    lines.append("-" * 70)
    lines.append("STEP 4: REGIONAL DISTRIBUTION SELECTION (Z-statistic)")
    lines.append("-" * 70)
    lines.append("|Z| <= 1.64 is conventionally taken as an adequate fit (~90% level).")
    lines.append("")
    lines.append(result.zstat_df[["family", "label", "t4_theoretical", "Z", "acceptable"]]
                  .to_string(index=False))
    lines.append("")
    lines.append(f"Chosen family: {result.growth_curve.label} ({result.chosen_family})")
    lines.append("")

    # -- Growth curve --
    lines.append("-" * 70)
    lines.append("STEP 5: REGIONAL GROWTH CURVE & STATION DESIGN FLOODS")
    lines.append("-" * 70)
    lines.append(f"Fitted to pooled, record-length-weighted L-moments: "
                  f"t_R (L-CV) = {result.growth_curve.t_R:.4f}, "
                  f"t3_R (L-skew) = {result.growth_curve.t3_R:.4f}")
    lines.append("Station design flood at return period T = growth_factor(T) * station's own")
    lines.append("index flood (its sample mean).")
    lines.append("")
    lines.append(result.quantile_table.round(2).to_string(index=False))
    lines.append("")
    lines.append("Note: as with the single-station tool, treat quantiles far beyond the pooled")
    lines.append("record length as increasingly uncertain — pooling improves the *shape* estimate")
    lines.append("but does not manufacture additional information about how far the tail extends.")
    lines.append("")

    if ci_df is not None:
        lines.append("-" * 70)
        lines.append(f"STEP 6: {confidence_level:g}% CONFIDENCE INTERVALS (Monte-Carlo)")
        lines.append("-" * 70)
        lines.append("Combines two independent, simulated sources of uncertainty: (1) regional")
        lines.append("growth-curve shape uncertainty, via synthetic homogeneous regions simulated")
        lines.append("from the same kappa distribution used for the H1/Z-statistic (Hosking & Wallis,")
        lines.append("1997, sec. 5.3), and (2) each station's own index-flood (sample mean) sampling")
        lines.append("uncertainty, via a nonparametric bootstrap of that station's own record. This")
        lines.append("is the regional counterpart of the single-station tool's bootstrap_ci_*.csv;")
        lines.append("see growth_curve_quantiles_ci.csv for the full table.")
        lines.append("")
        lines.append(ci_df.round(2).to_string(index=False))
        lines.append("")

    return "\n".join(lines)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("region_name", nargs="?", default=None,
                    help="Region name; expects Data/Regional/<RegionName>/ to contain one CSV "
                         "per station. If omitted, you'll be prompted to pick interactively.")
    p.add_argument("--value-col", default=None, help="Column name of the flood series (auto-detected if omitted)")
    p.add_argument("--year-col", default=None, help="Column name of the year index (auto-detected if omitted)")
    p.add_argument("--candidates", nargs="+", default=None,
                    choices=list(CANDIDATE_FAMILIES),
                    help=f"Candidate regional distribution families to test with the Z-statistic "
                         f"(default: all of {list(CANDIDATE_FAMILIES)}).")
    p.add_argument("--family", default=None, choices=list(CANDIDATE_FAMILIES),
                    help="Force a specific regional distribution instead of using the "
                         "Z-statistic recommendation (e.g. to match a published study).")
    p.add_argument("--n-sim", type=int, default=500,
                    help="Monte-Carlo replicates for the heterogeneity measure and Z-statistics "
                         "(default: 500). Higher is more stable but slower.")
    p.add_argument("--seed", type=int, default=None,
                    help="Random seed for the Monte-Carlo simulation (default: unseeded/random). "
                         "Set this for reproducible H1/Z values across runs.")
    p.add_argument("--return-periods", type=float, nargs="+", default=None,
                    help=f"Return periods (years) for the growth-curve/quantile table "
                         f"(default: {DEFAULT_RETURN_PERIODS}).")
    p.add_argument("--confidence-level", type=float, default=95.0,
                    help="Confidence level (percent) for the Monte-Carlo design-flood interval "
                         "(default: 95).")
    p.add_argument("--no-ci", action="store_true",
                    help="Skip the Monte-Carlo confidence-interval step (it re-runs the "
                         "growth-curve simulation plus a per-station bootstrap, so it roughly "
                         "doubles run time; skip it for a quick point-estimate-only run).")
    p.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation")
    return p.parse_args()


def main():
    args = parse_args()

    region_name = args.region_name
    if region_name is None:
        root = project_root_from(__file__)
        region_name = prompt_region_name(root)
        print()

    paths = resolve_region(region_name, __file__)

    print(f"Region:      {paths.region_name}")
    print(f"Project root:{paths.project_root}")
    print(f"Data dir:    {paths.data_dir}")
    print(f"Output dir:  {paths.output_dir}")
    print(f"Plot dir:    {paths.plot_dir}")
    print()

    try:
        station_data = load_region_stations(paths, value_col=args.value_col, year_col=args.year_col)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"ERROR reading regional station data: {e}")
    print(f"Loaded {len(station_data)} station(s): {', '.join(station_data)}")
    for name, values in station_data.items():
        print(f"  {name}: n={len(values)}, mean={values.mean():.1f}")
    print()

    candidates = tuple(args.candidates) if args.candidates else ("glo", "gev", "gno", "pe3", "gpa")
    return_periods = tuple(args.return_periods) if args.return_periods else tuple(DEFAULT_RETURN_PERIODS)

    try:
        result = run_regional_analysis(
            region_name, station_data, candidates=candidates,
            n_sim=args.n_sim, seed=args.seed, return_periods=return_periods,
            family=args.family)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    dq_trend = result.data_quality_df[result.data_quality_df["mann_kendall_significant"]]
    dq_outlier = result.data_quality_df[
        result.data_quality_df["grubbs_high_outlier_flagged"] | result.data_quality_df["grubbs_low_outlier_flagged"]]
    if not dq_trend.empty:
        print(f"! {len(dq_trend)} station(s) show a significant trend: {', '.join(dq_trend['station'])}")
    if not dq_outlier.empty:
        print(f"! {len(dq_outlier)} station(s) have a Grubbs-flagged outlier: {', '.join(dq_outlier['station'])}")
    disc_flagged = result.discordancy_df[result.discordancy_df["discordant"]]
    if not disc_flagged.empty:
        print(f"! {len(disc_flagged)} station(s) flagged discordant: "
              f"{', '.join(disc_flagged['station'])}")
    het = result.heterogeneity_result
    print(f"Heterogeneity: H1 = {het.H1:.3f} -> {het.interpretation}")
    print(f"Chosen regional distribution: {result.growth_curve.label}")
    print()

    # ---- CSV outputs ----
    result.stations_df.to_csv(paths.output_dir / "station_lmoments.csv", index=False)
    result.data_quality_df.to_csv(paths.output_dir / "data_quality.csv", index=False)
    result.discordancy_df.to_csv(paths.output_dir / "discordancy.csv", index=False)
    pd.DataFrame([{
        "H1": het.H1, "V1_obs": het.V1_obs, "V1_sim_mean": het.V1_sim_mean,
        "V1_sim_std": het.V1_sim_std, "n_sim": het.n_sim,
        "simulation_family": het.simulation_family, "interpretation": het.interpretation,
    }]).to_csv(paths.output_dir / "heterogeneity.csv", index=False)
    result.zstat_df.to_csv(paths.output_dir / "zstatistics.csv", index=False)
    result.quantile_table.to_csv(paths.output_dir / "quantile_table.csv", index=False)

    ci_df = None
    if not args.no_ci:
        if not (0 < args.confidence_level < 100):
            sys.exit(f"ERROR: --confidence-level must be between 0 and 100 (got {args.confidence_level}).")
        alpha = 1.0 - args.confidence_level / 100.0
        print(f"Computing {args.confidence_level:g}% Monte-Carlo confidence intervals "
              f"({args.n_sim} replicates x {len(station_data)} stations)...")
        ci_df = regional_quantile_ci(result, return_periods=return_periods,
                                      n_sim=args.n_sim, seed=args.seed, alpha=alpha)
        ci_df.to_csv(paths.output_dir / "growth_curve_quantiles_ci.csv", index=False)

    # ---- Summary ----
    summary = build_summary_text(result, build_provenance_header(), ci_df=ci_df,
                                  confidence_level=args.confidence_level)
    print(summary)
    summary_path = paths.output_dir / "summary.txt"
    summary_path.write_text(summary)
    print(f"Summary written to {summary_path}")

    # ---- Plots ----
    if not args.no_plots:
        save_regional_moment_ratio_diagram(
            result.stations, paths.plot_dir / "regional_moment_ratio_diagram.png",
            region_name=region_name, discordancy_df=result.discordancy_df)
        save_regional_growth_curve_plot(
            result, paths.plot_dir / "regional_growth_curve.png")
        save_regional_pooled_vs_stations_plot(
            result, paths.plot_dir / "regional_pooled_vs_stations.png")
        save_regional_discordancy_plot(
            result, paths.plot_dir / "regional_discordancy.png")
        save_regional_station_series_plot(
            result, paths.plot_dir / "regional_station_series.png")
        for name in station_data:
            save_station_design_flood_plot(
                result, name, paths.plot_dir / f"station_{name}_design_flood.png")
        save_regional_dashboard(
            result, paths.plot_dir / "regional_dashboard.png")
        print(f"PNG plots written to {paths.plot_dir}")

    print("\nDone.")


if __name__ == "__main__":
    main()
