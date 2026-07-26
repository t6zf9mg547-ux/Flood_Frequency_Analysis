#!/usr/bin/env python3
"""
Run a full flood frequency analysis for one case.

Layout expected (this script lives in <project_root>/Module/):

    <project_root>/
        Data/<CaseName>.csv
        Module/run_analysis.py   <- this file
        Module/floodfreq/        <- the package
        Output/<CaseName>/       <- CSV results written here
        Plot/<CaseName>/         <- PNG plots written here

Usage:
    python run_analysis.py P1009
    python run_analysis.py P1009 --value-col Q --year-col year
    python run_analysis.py P1009 --plotting-position gringorten --n-boot 2000
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

from floodfreq.analysis import FloodFrequencyAnalysis
from floodfreq.distributions import DISTRIBUTIONS
from floodfreq.io_utils import read_series, resolve_case, write_report, project_root_from
from floodfreq.plots import (
    save_probability_plot,
    save_quantile_ci_plot,
    save_moment_ratio_diagram,
    save_data_histogram,
)

DEFAULT_RETURN_PERIODS = [2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000]


def _csv_hint(path: Path) -> str:
    """Short descriptive hint for a CSV file: row count + last-modified time,
    so near-identical filenames (or a leftover Template.csv) are easy to
    tell apart in the picker."""
    import datetime
    try:
        with open(path) as f:
            n_rows = max(sum(1 for _ in f) - 1, 0)  # minus header
    except Exception:
        n_rows = "?"
    mtime = datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    return f"{n_rows} rows, modified {mtime}"


def prompt_case_name(project_root: Path) -> str:
    """
    List every CSV in Data/ and let the user pick one interactively.
    Used when case_name isn't passed on the command line (e.g. when running
    the script directly from an editor's "Run" button).
    """
    data_dir = project_root / "Data"
    csvs = sorted(data_dir.glob("*.csv"))
    if not csvs:
        sys.exit(f"ERROR: no .csv files found in {data_dir}. "
                  f"Add a file there (e.g. Data/<CaseName>.csv) and try again.")

    if len(csvs) == 1:
        only = csvs[0]
        raw = input(f"Found one input file: {only.name} ({_csv_hint(only)}) — use it? [Y/n]: ").strip().lower()
        if raw in ("", "y", "yes"):
            return only.stem
        # fall through to the full picker if they said no

    print(f"\nAvailable input files in {data_dir}:")
    for i, f in enumerate(csvs, start=1):
        print(f"  {i}. {f.name}  ({_csv_hint(f)})")

    while True:
        raw = input(f"Choose a file [1-{len(csvs)}] or type a case name: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(csvs):
            return csvs[int(raw) - 1].stem
        # allow typing the case name directly (with or without .csv)
        candidate = raw[:-4] if raw.lower().endswith(".csv") else raw
        if (data_dir / f"{candidate}.csv").exists():
            return candidate
        print(f"  Not a valid choice. Enter a number 1-{len(csvs)}, or an exact case name.")


def prompt_confidence_level(default=95.0) -> float:
    """Interactively ask the user for a confidence level (in percent),
    re-prompting on invalid input. Used when --confidence-level isn't passed."""
    while True:
        raw = input(f"Confidence interval level in % (e.g. 95) [default {default:g}]: ").strip()
        if raw == "":
            return default
        try:
            value = float(raw.rstrip("%"))
        except ValueError:
            print("  Please enter a number, e.g. 95 or 90.")
            continue
        if not (0 < value < 100):
            print("  Please enter a value strictly between 0 and 100.")
            continue
        return value


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("case_name", nargs="?", default=None,
                    help="Case name; expects Data/<CaseName>.csv to exist. "
                         "If omitted, you'll be prompted to pick from Data/ interactively.")
    p.add_argument("--value-col", default=None, help="Column name of the flood series (auto-detected if omitted)")
    p.add_argument("--year-col", default=None, help="Column name of the year index (auto-detected if omitted)")
    p.add_argument("--plotting-position", default=None,
                   choices=["weibull", "hazen", "cunnane", "gringorten", "hosking", "blom"],
                   help="Override the empirical frequency formula recorded/used for probability-"
                        "plot placement for EVERY distribution (default: each distribution uses its "
                        "recommended formula — Blom for Normal/LogNormal, Gringorten for the rest). "
                        "NOTE: PWM parameter estimation always uses the standard unbiased L-moment "
                        "estimator regardless of this setting — see summary.txt for details.")
    p.add_argument("--descriptive-plotting-position", default="weibull",
                   choices=["weibull", "hazen", "cunnane", "gringorten", "hosking", "blom"],
                   help="Plotting-position formula used ONLY for the general descriptive "
                        "statistics (Stat-sheet-style moments/PWM), independent of the "
                        "per-distribution fitting choice above (default: weibull)")
    p.add_argument("--n-boot", type=int, default=1000, help="Bootstrap resamples for the CI plot (default: 1000)")
    p.add_argument("--confidence-level", type=float, default=None,
                   help="Confidence level for the design-flood interval, in percent "
                        "(e.g. 95 for a 95%% CI). If omitted, you'll be prompted interactively.")
    p.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation")
    p.add_argument("--xlsx-report", action="store_true",
                   help="Also write a formatted .xlsx report to Output/<CaseName>/ (in addition to the CSVs)")
    return p.parse_args()


def main():
    args = parse_args()

    case_name = args.case_name
    if case_name is None:
        root = project_root_from(__file__)
        case_name = prompt_case_name(root)
        print()

    paths = resolve_case(case_name, __file__)

    if not paths.data_csv.exists():
        sys.exit(f"ERROR: expected input file not found: {paths.data_csv}")

    print(f"Case:        {paths.case_name}")
    print(f"Project root:{paths.project_root}")
    print(f"Input data:  {paths.data_csv}")
    print(f"Output dir:  {paths.output_dir}")
    print(f"Plot dir:    {paths.plot_dir}")
    print()

    values, years = read_series(paths.data_csv, value_col=args.value_col, year_col=args.year_col)
    print(f"Loaded {len(values)} annual maxima"
          f"{f' ({years.min()}-{years.max()})' if years is not None else ''}.")

    ffa = FloodFrequencyAnalysis(values, station_id=case_name, years=years)
    ffa.fit_all(plotting_position=args.plotting_position)  # None -> per-distribution recommendations

    confidence_level = args.confidence_level
    if confidence_level is None:
        confidence_level = prompt_confidence_level(default=95.0)
    elif not (0 < confidence_level < 100):
        sys.exit(f"ERROR: --confidence-level must be between 0 and 100 (got {confidence_level}).")
    print(f"Using a {confidence_level:g}% confidence interval.\n")

    # ---- CSV outputs ----
    stats = ffa.descriptive_stats(plotting_position=args.descriptive_plotting_position)
    import pandas as pd
    pd.Series(stats, name="value").rename_axis("statistic").reset_index().to_csv(
        paths.output_dir / "descriptive_stats.csv", index=False)

    gof = ffa.goodness_of_fit_table()
    gof.to_csv(paths.output_dir / "goodness_of_fit.csv", index=False)

    q_table = ffa.quantile_table(return_periods=DEFAULT_RETURN_PERIODS)
    q_table.to_csv(paths.output_dir / "quantile_table.csv", index=False)

    best_key, best_method = ffa.best_fit(criterion="AIC")

    # ---- Summary & recommendation (console + text file) ----
    # This also runs the bootstrap CI once; we re-derive the CSV from ffa's cache-friendly call below.
    recommendation = ffa.generate_recommendation(
        return_periods=DEFAULT_RETURN_PERIODS, criterion="AIC",
        confidence_level=confidence_level, ci_n_boot=args.n_boot, random_state=0)
    print("\n" + recommendation)
    summary_path = paths.output_dir / "summary.txt"
    summary_path.write_text(recommendation)
    print(f"\nSummary written to {summary_path}")

    alpha = 1.0 - confidence_level / 100.0
    ci_table = ffa.bootstrap_ci(best_key, best_method, DEFAULT_RETURN_PERIODS,
                                 n_boot=args.n_boot, alpha=alpha, random_state=0)
    ci_table.to_csv(paths.output_dir / f"bootstrap_ci_{best_key}.csv", index=False)

    if args.xlsx_report:
        report_path = paths.output_dir / f"{case_name}_report.xlsx"
        write_report(report_path, ffa, return_periods=DEFAULT_RETURN_PERIODS)
        print(f"Excel report written to {report_path}")

    # ---- Plots ----
    if not args.no_plots:
        save_data_histogram(ffa, paths.plot_dir / "input_data_histogram.png",
                             dist_key=best_key, method=best_method)
        save_probability_plot(ffa, paths.plot_dir / "flood_frequency_curve.png")
        save_moment_ratio_diagram(ffa, paths.plot_dir / "moment_ratio_diagram.png")
        save_quantile_ci_plot(ffa, best_key, best_method,
                               paths.plot_dir / f"bestfit_{best_key}_ci.png",
                               n_boot=args.n_boot, alpha=alpha)
        print(f"PNG plots written to {paths.plot_dir}")

    print("\nDone.")


if __name__ == "__main__":
    main()