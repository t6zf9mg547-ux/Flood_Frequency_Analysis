#!/usr/bin/env python3
"""
Run a CIFAM climate-informed flood adjustment for one case and scenario.

Implements the rapid Climate-Informed Flood Assessment Methodology (CIFAM) of
Grijsen & Lino (ICOLD 2026): takes a baseline single-station annual-maximum
series plus four projected climate-change numbers, and produces climate-adjusted
flood quantiles whose confidence intervals combine sampling uncertainty with
climate-change uncertainty.

Layout expected (this script lives in <project_root>/Module/):

    <project_root>/
        Data/<CaseName>/<CaseName>.csv                        <- baseline series
                                                                 (the SAME file run_analysis.py uses)
        Data/<CaseName>/climate_adjustment/<scenario>.csv     <- climate numbers for one scenario
        Data/Templates/                                       <- shipped example templates
        Module/run_climate_adjustment.py                      <- this file
        Module/floodfreq/                                     <- the package
        Output/<CaseName>/Climate_Adjustment/<scenario>/      <- CSV + summary.txt written here
        Plot/<CaseName>/Climate_Adjustment/<scenario>/        <- PNG plot written here

Setup: copy a climate-adjustment template from Data/Templates/ into
Data/<CaseName>/climate_adjustment/<scenario>.csv (one file per scenario, e.g.
rcp45.csv, rcp85.csv) and edit the four delta1/delta2/tau1/tau2 values. The
baseline series is the same Data/<CaseName>/<CaseName>.csv the single-station
tool uses. See the README's "Climate-informed adjustment (CIFAM)" section.

Usage:
    python run_climate_adjustment.py LomPangar rcp85
    python run_climate_adjustment.py LomPangar rcp85 --distribution lognormal2
    python run_climate_adjustment.py LomPangar rcp85 --method monte-carlo --n-sim 50000
    python run_climate_adjustment.py LomPangar rcp85 --confidence-level 90 --no-plot
    python run_climate_adjustment.py LomPangar rcp85 --delta1 0.30 --delta2 0.18 --tau1 0.10 --tau2 0.18

If <CaseName> and/or <scenario> are omitted you'll be prompted to pick them.
CLI flags override the corresponding values read from the inputs file, so a
one-off variation can be run without editing the CSV.
"""
from __future__ import annotations
import argparse
import sys
import datetime
import platform
from pathlib import Path

import numpy as np
import pandas as pd

import floodfreq
from floodfreq.io_utils import read_series, resolve_climate_case, load_climate_inputs
from floodfreq.climate_adjustment import (
    climate_adjusted_quantiles,
    mc_climate_adjusted_quantiles,
    CLOSED_FORM_DISTRIBUTIONS,
    MC_SUPPORTED_DISTRIBUTIONS,
)
from floodfreq.plots import save_climate_adjustment_plot


def build_provenance_header() -> str:
    """Short header stamping how this run was produced (mirrors run_analysis.py),
    so a summary.txt found later is self-documenting."""
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    command = " ".join([Path(sys.argv[0]).name] + sys.argv[1:])
    return "\n".join([
        "=" * 70,
        "RUN PROVENANCE",
        "=" * 70,
        f"Generated:     {now}",
        f"Tool version:  floodfreq {floodfreq.__version__}",
        f"Command:       python {command}",
        f"Python:        {platform.python_version()}",
        "=" * 70,
        "",
    ])


def _csv_hint(path: Path) -> str:
    import datetime as _dt
    try:
        with open(path) as f:
            n_rows = max(sum(1 for _ in f) - 1, 0)
    except Exception:
        n_rows = "?"
    mtime = _dt.datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    return f"{n_rows} rows, modified {mtime}"


def prompt_case_name(project_root: Path) -> str:
    """
    List Data/ case folders that have a climate_adjustment/ subfolder and let
    the user pick one. Used when case_name isn't passed on the command line.
    Reserved folders (Templates/, Regional/) are skipped.
    """
    data_dir = project_root / "Data"
    reserved = {"Templates", "Regional"}
    cases = sorted(
        d for d in data_dir.iterdir()
        if d.is_dir() and d.name not in reserved
        and (d / "climate_adjustment").is_dir()
        and any((d / "climate_adjustment").glob("*.csv"))
    ) if data_dir.is_dir() else []
    if not cases:
        sys.exit(f"ERROR: no cases with a climate_adjustment/ folder found in "
                 f"{data_dir}. Create Data/<CaseName>/climate_adjustment/"
                 f"<scenario>.csv (copy a template from Data/Templates/) first.")

    if len(cases) == 1:
        only = cases[0]
        raw = input(f"Found one case with scenarios: {only.name} — use it? "
                    f"[Y/n]: ").strip().lower()
        if raw in ("", "y", "yes"):
            return only.name

    print(f"\nCases with climate scenarios in {data_dir}:")
    for i, d in enumerate(cases, start=1):
        n = len(list((d / "climate_adjustment").glob("*.csv")))
        print(f"  {i}. {d.name}  ({n} scenario{'s' if n != 1 else ''})")
    while True:
        raw = input(f"Choose a case [1-{len(cases)}] or type a case name: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(cases):
            return cases[int(raw) - 1].name
        if (data_dir / raw / "climate_adjustment").is_dir():
            return raw
        print(f"  Not a valid choice. Enter a number 1-{len(cases)}, or an exact case name.")


def prompt_scenario(project_root: Path, case_name: str) -> str:
    """
    List the scenario CSVs in Data/<CaseName>/climate_adjustment/ and let the
    user pick one. Used when the scenario isn't passed on the command line.
    """
    scen_dir = project_root / "Data" / case_name / "climate_adjustment"
    if not scen_dir.is_dir():
        sys.exit(f"ERROR: {scen_dir} does not exist. Copy a climate-adjustment "
                 f"template from Data/Templates/ into "
                 f"Data/{case_name}/climate_adjustment/<scenario>.csv first.")
    csvs = sorted(scen_dir.glob("*.csv"))
    if not csvs:
        sys.exit(f"ERROR: no scenario .csv files in {scen_dir}. Add one "
                 f"(copy a template from Data/Templates/) and edit the four values.")

    if len(csvs) == 1:
        only = csvs[0]
        raw = input(f"Found one scenario: {only.stem} "
                    f"({_csv_hint(only)}) — use it? [Y/n]: ").strip().lower()
        if raw in ("", "y", "yes"):
            return only.stem

    print(f"\nScenarios in {scen_dir}:")
    for i, f in enumerate(csvs, start=1):
        print(f"  {i}. {f.stem}  ({_csv_hint(f)})")
    while True:
        raw = input(f"Choose a scenario [1-{len(csvs)}] or type its name: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(csvs):
            return csvs[int(raw) - 1].stem
        candidate = raw[:-4] if raw.lower().endswith(".csv") else raw
        if (scen_dir / f"{candidate}.csv").exists():
            return candidate
        print(f"  Not a valid choice. Enter a number 1-{len(csvs)}, or an exact scenario name.")


def _override(cli_value, file_value):
    """A CLI flag (when passed, i.e. not None) overrides the inputs-file value."""
    return cli_value if cli_value is not None else file_value


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("case_name", nargs="?", default=None,
                   help="Case name; reads Data/<CaseName>/<CaseName>.csv (baseline). "
                        "If omitted, you'll be prompted to pick one.")
    p.add_argument("scenario", nargs="?", default=None,
                   help="Scenario name; reads "
                        "Data/<CaseName>/climate_adjustment/<scenario>.csv (climate inputs) "
                        "and writes to Output/<CaseName>/Climate_Adjustment/<scenario>/. "
                        "If omitted, you'll be prompted to pick one.")
    p.add_argument("--value-col", default=None,
                   help="Column name of the flood series in the baseline CSV (auto-detected if omitted)")
    p.add_argument("--year-col", default=None,
                   help="Column name of the year index in the baseline CSV (auto-detected if omitted)")
    p.add_argument("--variable-name", default=None,
                   help="Axis label for the plotted variable (default: 'Flood magnitude')")
    p.add_argument("--units", default=None, help="Units suffix for the axis label, e.g. m3/s")
    p.add_argument("--short-name", default=None, help="Short name used in the plot title (default: 'Flood')")
    # climate inputs -- CLI overrides for the values in the inputs file
    p.add_argument("--delta1", type=float, default=None, help="Override: fractional increase in the mean")
    p.add_argument("--delta2", type=float, default=None, help="Override: ensemble sd of delta1")
    p.add_argument("--tau1", type=float, default=None, help="Override: fractional increase in the std dev")
    p.add_argument("--tau2", type=float, default=None, help="Override: ensemble sd of tau1")
    p.add_argument("--distribution", default=None,
                   help=f"Override the distribution. Closed form: {', '.join(CLOSED_FORM_DISTRIBUTIONS)}. "
                        f"With --method monte-carlo also: "
                        f"{', '.join(d for d in MC_SUPPORTED_DISTRIBUTIONS if d not in CLOSED_FORM_DISTRIBUTIONS)}.")
    p.add_argument("--confidence-level", type=float, default=None,
                   help="Override the two-sided CI level in percent (e.g. 90)")
    p.add_argument("--pmf", type=float, default=None, help="Override the PMF reference value")
    p.add_argument("--return-periods", default=None,
                   help="Override return periods, comma-separated (e.g. 100,1000,10000)")
    p.add_argument("--method", choices=["closed-form", "monte-carlo"], default="closed-form",
                   help="closed-form (default; Gumbel/Log-Normal/Pearson III only) or monte-carlo (any distribution)")
    p.add_argument("--n-sim", type=int, default=30000, help="Monte Carlo simulations (--method monte-carlo only)")
    p.add_argument("--random-state", type=int, default=0, help="Monte Carlo RNG seed (reproducibility)")
    p.add_argument("--no-plot", action="store_true", help="Skip PNG plot generation")
    p.add_argument("--return-period-axis", action="store_true",
                   help="Plot against a log return-period axis instead of the Gumbel reduced variate")
    return p.parse_args()


def main():
    args = parse_args()
    case_name = args.case_name
    scenario = args.scenario

    # Resolve paths first (needs a project root); prompt for case/scenario if missing.
    from floodfreq.io_utils import project_root_from
    project_root = project_root_from(__file__)
    if case_name is None:
        case_name = prompt_case_name(project_root)
    if scenario is None:
        scenario = prompt_scenario(project_root, case_name)

    paths = resolve_climate_case(case_name, scenario, __file__)

    # -- baseline series --
    if not paths.baseline_csv.exists():
        sys.exit(f"ERROR: baseline series not found: {paths.baseline_csv}\n"
                 f"CIFAM reuses the single-station input; create it (or run "
                 f"run_analysis.py on this case) first.")
    Q, years = read_series(paths.baseline_csv, value_col=args.value_col, year_col=args.year_col)

    # -- climate inputs (from file, then CLI overrides) --
    ci = load_climate_inputs(paths.climate_csv)
    delta1 = _override(args.delta1, ci["delta1"])
    delta2 = _override(args.delta2, ci["delta2"])
    tau1 = _override(args.tau1, ci["tau1"])
    tau2 = _override(args.tau2, ci["tau2"])
    distribution = _override(args.distribution, ci["distribution"])
    confidence_level = _override(args.confidence_level, ci["confidence_level"])
    pmf = _override(args.pmf, ci["pmf"])
    if args.return_periods is not None:
        return_periods = tuple(float(x) for x in args.return_periods.split(",") if x.strip())
    else:
        return_periods = ci["return_periods"]

    variable_name = args.variable_name or "Flood magnitude"
    units = args.units
    short_name = args.short_name or "Flood"

    print(f"\nClimate-informed flood adjustment (CIFAM) for case '{case_name}', scenario '{scenario}'")
    print(f"  baseline series : {paths.baseline_csv}  (N = {Q.size})")
    print(f"  climate inputs  : {paths.climate_csv}")
    print(f"  distribution    : {distribution}   method: {args.method}")
    print(f"  delta1={delta1:g} delta2={delta2:g} tau1={tau1:g} tau2={tau2:g}  "
          f"CI={confidence_level:g}%")

    # -- run --
    if args.method == "closed-form":
        if distribution not in CLOSED_FORM_DISTRIBUTIONS:
            sys.exit(f"ERROR: '{distribution}' has no closed form. Use one of "
                     f"{CLOSED_FORM_DISTRIBUTIONS}, or add --method monte-carlo "
                     f"(supported: {', '.join(MC_SUPPORTED_DISTRIBUTIONS)}).")
        result = climate_adjusted_quantiles(
            Q, distribution, return_periods, delta1, delta2, tau1, tau2,
            confidence_level=confidence_level)
    else:
        if distribution not in MC_SUPPORTED_DISTRIBUTIONS:
            sys.exit(f"ERROR: Monte Carlo climate adjustment is not defined for "
                     f"'{distribution}'. Supported: {', '.join(MC_SUPPORTED_DISTRIBUTIONS)}.")
        result = mc_climate_adjusted_quantiles(
            Q, distribution, return_periods, delta1, delta2, tau1, tau2,
            confidence_level=confidence_level, n_sim=args.n_sim,
            random_state=args.random_state)

    # -- outputs --
    table = result.to_frame()
    csv_path = paths.output_dir / "climate_adjustment_table.csv"
    table.to_csv(csv_path, index=False)
    print(f"\nResults table written to {csv_path}")

    summary = build_provenance_header()
    summary += _format_summary(case_name, scenario, result, args.method, delta1, delta2,
                               tau1, tau2, pmf, Q, distribution)
    summary_path = paths.output_dir / "summary.txt"
    summary_path.write_text(summary)
    print(f"Summary written to {summary_path}")

    if not args.no_plot:
        plot_path = paths.plot_dir / "climate_adjustment.png"
        save_climate_adjustment_plot(
            result, plot_path, pmf=pmf, variable_name=variable_name,
            units=units, short_name=short_name,
            reduced_variate_axis=not args.return_period_axis)
        print(f"Plot written to {plot_path}")

    print("\nDone.")


def _format_summary(case_name, scenario, result, method, delta1, delta2, tau1, tau2,
                    pmf, Q, distribution) -> str:
    from floodfreq.distributions import DISTRIBUTIONS
    lines = [
        "CLIMATE-INFORMED FLOOD ADJUSTMENT (CIFAM)",
        "Grijsen & Lino, ICOLD 2026",
        "=" * 70,
        f"Case:            {case_name}",
        f"Scenario:        {scenario}",
        f"Distribution:    {DISTRIBUTIONS[distribution]['label']}",
        f"Method:          {method}",
        f"Baseline record: N = {Q.size}, mean = {Q.mean():.1f}, std = {Q.std():.1f}",
        "",
        "Climate-change inputs (flow-space, fractions):",
        f"  delta1 (mean shift)          = {delta1:.3f}  (+/- {delta2:.3f})",
        f"  tau1   (std-dev shift)       = {tau1:.3f}  (+/- {tau2:.3f})",
        f"  confidence level             = {result.confidence_level:g}%",
    ]
    if pmf is not None:
        lines.append(f"  PMF reference                = {pmf:g}")
    lines += [
        "",
        "The 'baseline' columns reflect sampling uncertainty only (classical",
        "FFA). The 'climate' columns add climate-change uncertainty: a shifted",
        "central estimate and a confidence interval combining both sources.",
        "",
        "-" * 70,
    ]
    df = result.to_frame()
    header = (f"{'T':>9}  {'base_pt':>10}  {'base_lo':>10}  {'base_hi':>10}  "
              f"{'clim_pt':>10}  {'clim_lo':>10}  {'clim_hi':>10}")
    lines.append(header)
    lines.append("-" * len(header))
    for _, r in df.iterrows():
        lines.append(f"{r['T']:>9.0f}  {r['baseline_point']:>10.1f}  "
                     f"{r['baseline_lower']:>10.1f}  {r['baseline_upper']:>10.1f}  "
                     f"{r['climate_point']:>10.1f}  {r['climate_lower']:>10.1f}  "
                     f"{r['climate_upper']:>10.1f}")
    if pmf is not None:
        long_T = df.iloc[-1]
        lines += [
            "-" * len(header),
            "",
            f"At the longest return period (T = {long_T['T']:.0f} yr):",
            f"  climate {result.confidence_level:g}% upper = {long_T['climate_upper']:.0f}"
            f"  vs  PMF = {pmf:g}"
            f"  ({100 * (long_T['climate_upper'] / pmf - 1):+.0f}% relative to PMF)",
        ]
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    main()
