#!/usr/bin/env python3
"""
Run a full flood frequency analysis for one case.

Layout expected (this script lives in <project_root>/Module/):

    <project_root>/
        Data/<CaseName>/<CaseName>.csv
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
import datetime
import platform
from pathlib import Path

import pandas as pd

import floodfreq
from floodfreq.analysis import FloodFrequencyAnalysis
from floodfreq.distributions import DISTRIBUTIONS
from floodfreq.io_utils import read_series, resolve_case, write_report, project_root_from, load_case_config
from floodfreq.plots import (
    save_probability_plot,
    save_quantile_ci_plot,
    save_moment_ratio_diagram,
    save_data_histogram,
    save_data_quality_plot,
    save_dashboard,
    save_pdf_report,
)

DEFAULT_RETURN_PERIODS = [2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000]


def build_provenance_header(config_path: Path = None) -> str:
    """
    A short header stamping exactly how this run was produced: when, with
    what command, and which tool version -- so a summary.txt or PDF found
    later (or handed to someone else) is self-documenting instead of an
    orphaned result with no record of the settings that produced it.
    """
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
    ]
    if config_path is not None:
        lines.append(f"Config file:   {config_path} (CLI flags, if any, override its settings)")
    lines.append("=" * 70)
    lines.append("")
    return "\n".join(lines)


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
    List every case folder in Data/ and let the user pick one interactively.
    Used when case_name isn't passed on the command line (e.g. when running
    the script directly from an editor's "Run" button).

    Case-first layout (v0.8.0): a case is a folder Data/<CaseName>/ containing
    a baseline series Data/<CaseName>/<CaseName>.csv. The reserved folders
    Templates/ (shipped examples) and Regional/ (regional analysis) are not
    cases and are skipped.
    """
    data_dir = project_root / "Data"
    reserved = {"Templates", "Regional"}
    cases = sorted(
        d for d in data_dir.iterdir()
        if d.is_dir() and d.name not in reserved and (d / f"{d.name}.csv").exists()
    ) if data_dir.is_dir() else []
    if not cases:
        sys.exit(f"ERROR: no case folders found in {data_dir}. "
                 f"Create Data/<CaseName>/<CaseName>.csv and try again.")

    if len(cases) == 1:
        only = cases[0]
        raw = input(f"Found one case: {only.name} "
                    f"({_csv_hint(only / f'{only.name}.csv')}) — use it? [Y/n]: ").strip().lower()
        if raw in ("", "y", "yes"):
            return only.name
        # fall through to the full picker if they said no

    print(f"\nAvailable cases in {data_dir}:")
    for i, d in enumerate(cases, start=1):
        print(f"  {i}. {d.name}  ({_csv_hint(d / f'{d.name}.csv')})")

    while True:
        raw = input(f"Choose a case [1-{len(cases)}] or type a case name: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(cases):
            return cases[int(raw) - 1].name
        candidate = raw
        if (data_dir / candidate / f"{candidate}.csv").exists():
            return candidate
        print(f"  Not a valid choice. Enter a number 1-{len(cases)}, or an exact case name.")


def resolve_setting(cli_value, config: dict, key: str, hardcoded_default):
    """
    Precedence: an explicit CLI flag wins, then the per-case Data/<CaseName>/<CaseName>.toml
    config file, then the hardcoded default. `cli_value` must be None when the
    flag wasn't passed on the command line (all mergeable CLI flags default to
    None for exactly this reason).
    """
    if cli_value is not None:
        return cli_value
    if key in config:
        return config[key]
    return hardcoded_default


def resolve_bool_setting(cli_flag: bool, config: dict, key: str) -> bool:
    """
    For store_true CLI flags: True if EITHER the CLI flag was passed OR the
    config file sets it. There's no way to force a config-file `true` back to
    False from the CLI with a plain store_true flag -- edit the .toml file
    directly if you need to turn it off for one run.
    """
    return bool(cli_flag) or bool(config.get(key, False))


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
                    help="Case name; expects Data/<CaseName>/<CaseName>.csv to exist. "
                         "If omitted, you'll be prompted to pick from Data/ interactively.")
    p.add_argument("--value-col", default=None, help="Column name of the flood series (auto-detected if omitted)")
    p.add_argument("--year-col", default=None, help="Column name of the year index (auto-detected if omitted)")
    p.add_argument("--variable-name", default=None,
                   help="Label for the analyzed variable, used on plot axes (default: 'Flood magnitude'). "
                        "Lets this tool be used for any annual-maximum series, not just streamflow -- "
                        "e.g. --variable-name 'Rainfall depth' --units mm for extreme rainfall analysis.")
    p.add_argument("--units", default=None,
                   help="Units for --variable-name, appended in parentheses on plot axes (e.g. 'mm', 'm3/s'). "
                        "No default -- omitted from labels if not given.")
    p.add_argument("--short-name", default=None,
                   help="Short label used in titles/headers (default: 'Flood', giving 'Flood frequency "
                        "curve', 'FLOOD FREQUENCY ANALYSIS', etc.). For rainfall, try --short-name Rainfall. "
                        "NOTE: the --regional-skew feature (Bulletin 17B) uses MSE coefficients empirically "
                        "calibrated from streamflow data specifically -- avoid it for non-streamflow variables "
                        "regardless of this setting.")
    p.add_argument("--plotting-position", default=None,
                   choices=["weibull", "hazen", "cunnane", "gringorten", "hosking", "blom"],
                   help="Override the empirical frequency formula recorded/used for probability-"
                        "plot placement for EVERY distribution (default: each distribution uses its "
                        "recommended formula — Blom for Normal/LogNormal, Gringorten for the rest). "
                        "NOTE: PWM parameter estimation always uses the standard unbiased L-moment "
                        "estimator regardless of this setting — see summary.txt for details.")
    p.add_argument("--descriptive-plotting-position", default=None,
                   choices=["weibull", "hazen", "cunnane", "gringorten", "hosking", "blom"],
                   help="Plotting-position formula used ONLY for the general descriptive "
                        "statistics (Stat-sheet-style moments/PWM), independent of the "
                        "per-distribution fitting choice above (default: weibull)")
    p.add_argument("--regional-skew", type=float, default=None,
                   help="Regional skew value for Pearson III / Log-Pearson III (Bulletin 17B "
                        "weighted-skew procedure). Must come from a published regional study "
                        "for your area -- there is no sensible default. If given, adds "
                        "weighted-skew variants of Pearson III/LP3 as additional candidates "
                        "alongside their standard PWM fits.")
    p.add_argument("--regional-skew-mse", type=float, default=None,
                   help="Mean square error of --regional-skew (default: 0.302, the documented "
                        "MSE of the Bulletin 17B national skew map -- prefer a region/state-"
                        "specific study's MSE if you have one, it will usually be smaller/better).")
    p.add_argument("--n-boot", type=int, default=None, help="Bootstrap resamples for the CI plot (default: 1000)")
    p.add_argument("--confidence-level", type=float, default=None,
                   help="Confidence level for the design-flood interval, in percent "
                        "(e.g. 95 for a 95%% CI). If omitted, you'll be prompted interactively.")
    p.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation")
    p.add_argument("--xlsx-report", action="store_true",
                   help="Also write a formatted .xlsx report to Output/<CaseName>/ (in addition to the CSVs)")
    p.add_argument("--pdf-report", action="store_true",
                   help="Also write a single bundled .pdf report (text summary + dashboard + all "
                        "individual plots) to Output/<CaseName>/ -- the one file to hand to a "
                        "colleague instead of the separate CSVs/PNGs.")
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

    try:
        config = load_case_config(case_name, paths.project_root)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")
    config_path = paths.project_root / "Data" / case_name / f"{case_name}.toml"
    if config:
        print(f"Using settings from {config_path} (CLI flags, if any, take precedence).")

    value_col = resolve_setting(args.value_col, config, "value_col", None)
    year_col = resolve_setting(args.year_col, config, "year_col", None)
    variable_name = resolve_setting(args.variable_name, config, "variable_name", "Flood magnitude")
    units = resolve_setting(args.units, config, "units", None)
    short_name = resolve_setting(args.short_name, config, "short_name", "Flood")
    plotting_position = resolve_setting(args.plotting_position, config, "plotting_position", None)
    descriptive_plotting_position = resolve_setting(
        args.descriptive_plotting_position, config, "descriptive_plotting_position", "weibull")
    regional_skew = resolve_setting(args.regional_skew, config, "regional_skew", None)
    regional_skew_mse = resolve_setting(args.regional_skew_mse, config, "regional_skew_mse", 0.302)
    n_boot = resolve_setting(args.n_boot, config, "n_boot", 1000)
    cli_confidence_level = resolve_setting(args.confidence_level, config, "confidence_level", None)
    no_plots = resolve_bool_setting(args.no_plots, config, "no_plots")
    xlsx_report = resolve_bool_setting(args.xlsx_report, config, "xlsx_report")
    pdf_report = resolve_bool_setting(args.pdf_report, config, "pdf_report")

    print(f"Case:        {paths.case_name}")
    print(f"Project root:{paths.project_root}")
    print(f"Input data:  {paths.data_csv}")
    print(f"Output dir:  {paths.output_dir}")
    print(f"Plot dir:    {paths.plot_dir}")
    print()

    try:
        values, years = read_series(paths.data_csv, value_col=value_col, year_col=year_col)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"ERROR reading input data: {e}")
    print(f"Loaded {len(values)} annual maxima"
          f"{f' ({years.min()}-{years.max()})' if years is not None else ''}.")

    try:
        ffa = FloodFrequencyAnalysis(values, station_id=case_name, years=years,
                                      variable_name=variable_name, units=units, short_name=short_name)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    # ---- Data quality checks (before fitting, so issues are visible early) ----
    dqr = ffa.data_quality()
    if dqr["validation_warnings"]:
        print("\nData quality warnings:")
        for w in dqr["validation_warnings"]:
            print(f"  ! {w}")
    mk = dqr["mann_kendall"]
    print(f"Stationarity check: {mk['trend']} (p = {mk['p_value']:.3f})")
    if mk["significant"]:
        print("  -> significant trend detected: the stationarity assumption behind this "
              "analysis is questionable; see summary.txt for details.")
    gr = dqr["grubbs"]
    if gr["high_outlier_flagged"] or gr["low_outlier_flagged"]:
        print("  -> possible outlier(s) flagged (Grubbs' test); see summary.txt for details.")

    dq_rows = [
        {"check": "n_records", "value": ffa.n,
         "meaning": f"You have {ffa.n} years of data in this record."},
        {"check": "validation_warnings",
         "value": " | ".join(dqr["validation_warnings"]) or "none",
         "meaning": "No input-quality issues found." if not dqr["validation_warnings"]
                    else "See the value column for details; worth a quick manual check "
                         "of the raw data before trusting the analysis."},
        {"check": "mann_kendall_trend", "value": mk["trend"],
         "meaning": ("No evidence that flood magnitudes are systematically increasing or "
                     "decreasing over time; the record looks stationary.")
                    if not mk["significant"] else
                    ("A statistically significant trend was detected — the assumption that "
                     "flood risk is constant through time may not hold here, which affects "
                     "how much to trust the fitted distribution and its extrapolation.")},
        {"check": "mann_kendall_p_value", "value": round(mk["p_value"], 4),
         "meaning": f"p={mk['p_value']:.3f}; a value below 0.05 would indicate a statistically "
                    f"significant trend. This is {'above' if mk['p_value'] >= 0.05 else 'below'} "
                    f"that threshold."},
        {"check": "grubbs_high_outlier_value", "value": gr["high_outlier_value"],
         "meaning": "The single largest value in the record."},
        {"check": "grubbs_high_outlier_flagged", "value": gr["high_outlier_flagged"],
         "meaning": "Good: this value is NOT considered a statistical outlier." if not gr["high_outlier_flagged"]
                    else "Warning: this value IS flagged as a possible outlier — worth checking "
                         "whether it's a data error or a genuine extreme event."},
        {"check": "grubbs_low_outlier_value", "value": gr["low_outlier_value"],
         "meaning": "The single smallest value in the record."},
        {"check": "grubbs_low_outlier_flagged", "value": gr["low_outlier_flagged"],
         "meaning": "Good: this value is NOT considered a statistical outlier." if not gr["low_outlier_flagged"]
                    else "Warning: this value IS flagged as a possible outlier — worth checking "
                         "whether it's a data error or a genuine extreme event."},
    ]
    pd.DataFrame(dq_rows).to_csv(paths.output_dir / "data_quality.csv", index=False)
    ffa.fit_all(plotting_position=plotting_position,  # None -> per-distribution recommendations
                regional_skew=regional_skew, regional_mse=regional_skew_mse)

    confidence_level = cli_confidence_level
    if confidence_level is None:
        confidence_level = prompt_confidence_level(default=95.0)
    elif not (0 < confidence_level < 100):
        sys.exit(f"ERROR: --confidence-level must be between 0 and 100 (got {confidence_level}).")
    print(f"Using a {confidence_level:g}% confidence interval.\n")

    # ---- CSV outputs ----
    stats = ffa.descriptive_stats(plotting_position=descriptive_plotting_position)
    pd.Series(stats, name="value").rename_axis("statistic").reset_index().to_csv(
        paths.output_dir / "descriptive_stats.csv", index=False)

    gof = ffa.goodness_of_fit_table()
    gof.to_csv(paths.output_dir / "goodness_of_fit.csv", index=False)

    q_table = ffa.quantile_table(return_periods=DEFAULT_RETURN_PERIODS)
    q_table.to_csv(paths.output_dir / "quantile_table.csv", index=False)

    ffa.akaike_weights().to_csv(paths.output_dir / "akaike_weights.csv", index=False)
    ffa.model_averaged_quantile_table(return_periods=DEFAULT_RETURN_PERIODS).to_csv(
        paths.output_dir / "model_averaged_quantiles.csv", index=False)

    best_key, best_method = ffa.best_fit(criterion="AIC")

    # ---- Summary & recommendation (console + text file) ----
    # This also runs the bootstrap CI once; we re-derive the CSV from ffa's cache-friendly call below.
    recommendation = ffa.generate_recommendation(
        return_periods=DEFAULT_RETURN_PERIODS, criterion="AIC",
        confidence_level=confidence_level, ci_n_boot=n_boot, random_state=0)
    recommendation = build_provenance_header(config_path if config else None) + "\n" + recommendation
    print("\n" + recommendation)
    summary_path = paths.output_dir / "summary.txt"
    summary_path.write_text(recommendation)
    print(f"\nSummary written to {summary_path}")

    alpha = 1.0 - confidence_level / 100.0
    ci_table = ffa.bootstrap_ci(best_key, best_method, DEFAULT_RETURN_PERIODS,
                                 n_boot=n_boot, alpha=alpha, random_state=0)
    ci_table.to_csv(paths.output_dir / f"bootstrap_ci_{best_key}.csv", index=False)

    if xlsx_report:
        report_path = paths.output_dir / f"{case_name}_report.xlsx"
        write_report(report_path, ffa, return_periods=DEFAULT_RETURN_PERIODS)
        print(f"Excel report written to {report_path}")

    if pdf_report:
        pdf_path = paths.output_dir / f"{case_name}_report.pdf"
        save_pdf_report(ffa, pdf_path, best_key, best_method, recommendation,
                         n_boot=n_boot, alpha=alpha)
        print(f"PDF report written to {pdf_path}")

    # ---- Plots ----
    if not no_plots:
        save_data_quality_plot(ffa, paths.plot_dir / "data_quality_timeseries.png")
        save_data_histogram(ffa, paths.plot_dir / "input_data_histogram.png",
                             dist_key=best_key, method=best_method)
        save_probability_plot(ffa, paths.plot_dir / "flood_frequency_curve.png")
        save_moment_ratio_diagram(ffa, paths.plot_dir / "moment_ratio_diagram.png")
        save_quantile_ci_plot(ffa, best_key, best_method,
                               paths.plot_dir / f"bestfit_{best_key}_ci.png",
                               n_boot=n_boot, alpha=alpha)
        save_dashboard(ffa, paths.plot_dir / "dashboard.png",
                        best_key, best_method, n_boot=n_boot, alpha=alpha)
        print(f"PNG plots written to {paths.plot_dir}")

    print("\nDone.")


if __name__ == "__main__":
    main()
