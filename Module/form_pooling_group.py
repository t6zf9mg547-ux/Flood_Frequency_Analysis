#!/usr/bin/env python3
"""
Propose a pooling group for regional flood/rainfall frequency analysis:
the "region of influence" approach (Burn, 1990) -- given a target site and
a catalog of candidate stations described by numeric descriptors (e.g.
catchment area, mean annual precipitation, a soil/permeability index for
streamflow gauges; geographic coordinates, elevation, point climatology
for rain gauges), rank candidates by similarity and propose which ones to
pool, instead of assembling the group by hand.

Do you need this? Only if you're choosing among a LARGER POOL of candidate
stations and want a data-driven suggestion for which ones to pool. If you
already know which stations belong together, skip this script entirely --
just build Data/Regional/<RegionName>/ directly (e.g. by copying
Data/Regional/Template/) and run run_regional_analysis.py on it.

This script is descriptor-SOURCE-agnostic: it doesn't extract descriptors
itself (e.g. from BasinATLAS/HydroSHEDS) -- it expects you've already
built a candidate catalog CSV with a `station` column and whatever numeric
descriptor columns you choose. See floodfreq/pooling.py's docstring for
the expected schema.

This is a PROPOSAL step, not a substitute for `run_regional_analysis.py`'s
own discordancy/heterogeneity screening -- always run that on the
resulting group before trusting it.

Usage:
    # Dry run: just see the ranking and proposed group
    python form_pooling_group.py --catalog candidates.csv \\
        --descriptors area_km2 mean_annual_precip_mm bfihost \\
        --target-station GAUGE_042 --n-stations 8 --region-name MyRegion

    # Ungauged target site, stop once 250 station-years are pooled
    python form_pooling_group.py --catalog candidates.csv \\
        --descriptors area_km2 mean_annual_precip_mm \\
        --target-descriptors area_km2=180 mean_annual_precip_mm=1050 \\
        --min-years 250 --region-name MyRegion

    # Actually copy the chosen stations' CSVs into Data/Regional/MyRegion/
    python form_pooling_group.py --catalog candidates.csv \\
        --descriptors area_km2 mean_annual_precip_mm \\
        --target-station GAUGE_042 --n-stations 8 --region-name MyRegion \\
        --station-data-dir /path/to/all_candidate_series/ --apply
"""
from __future__ import annotations
import argparse
import shutil
import sys
from pathlib import Path

from floodfreq.io_utils import resolve_region
from floodfreq.pooling import (
    read_candidate_catalog, resolve_target, similarity_ranking, propose_pooling_group,
)


def _parse_kv_list(pairs, what: str) -> dict:
    """Parse ['name=1.5', 'other=2'] into {'name': 1.5, 'other': 2.0}."""
    out = {}
    for p in pairs or []:
        if "=" not in p:
            sys.exit(f"ERROR: {what} entries must look like name=value (got '{p}').")
        name, _, value = p.partition("=")
        try:
            out[name] = float(value)
        except ValueError:
            sys.exit(f"ERROR: {what} value for '{name}' isn't numeric (got '{value}').")
    return out


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog", required=True,
                    help="CSV of candidate stations: a 'station' column, numeric descriptor "
                         "columns, and optionally an 'n_years' column (record length, needed "
                         "for --min-years).")
    p.add_argument("--descriptors", nargs="+", required=True,
                    help="Which columns in --catalog to use as descriptors, e.g. "
                         "area_km2 mean_annual_precip_mm bfihost")
    p.add_argument("--target-station", default=None,
                    help="Use an existing row in --catalog as the target (leave-one-out style: "
                         "'if I were building a region around this gauge, who looks similar?'). "
                         "Mutually exclusive with --target-descriptors.")
    p.add_argument("--target-descriptors", nargs="+", default=None,
                    help="Descriptor values for a target site NOT in the catalog (the usual "
                         "ungauged-target case), as name=value pairs matching --descriptors, "
                         "e.g. --target-descriptors area_km2=180 mean_annual_precip_mm=1050")
    p.add_argument("--weights", nargs="+", default=None,
                    help="Optional per-descriptor weights as name=value pairs (default: equal "
                         "weight, 1.0, for every descriptor). A larger weight makes that "
                         "descriptor count for more in the similarity ranking.")
    p.add_argument("--n-stations", type=int, default=None,
                    help="Stopping rule: take the N most similar candidates. Mutually exclusive "
                         "with --min-years.")
    p.add_argument("--min-years", type=int, default=None,
                    help="Stopping rule: accumulate the most similar candidates until their "
                         "combined record length reaches at least this many station-years "
                         "(Hosking & Wallis, 1997 suggest ~5x the design return period as a "
                         "rule of thumb). Requires an 'n_years' column in --catalog. Mutually "
                         "exclusive with --n-stations.")
    p.add_argument("--region-name", required=True,
                    help="Name for the proposed region; results are written to "
                         "Output/Regional/<RegionName>/pooling_*.csv")
    p.add_argument("--station-data-dir", default=None,
                    help="Directory containing each candidate's own annual-maximum series CSV "
                         "(named <station>.csv), needed only if --apply is given.")
    p.add_argument("--apply", action="store_true",
                    help="Actually copy the proposed stations' CSVs from --station-data-dir "
                         "into Data/Regional/<RegionName>/. Without this flag, the script only "
                         "reports the ranking/proposal (a dry run) -- nothing is written to "
                         "Data/.")
    return p.parse_args()


def main():
    args = parse_args()

    if (args.target_station is None) == (args.target_descriptors is None):
        sys.exit("ERROR: give exactly one of --target-station or --target-descriptors.")
    if (args.n_stations is None) == (args.min_years is None):
        sys.exit("ERROR: give exactly one of --n-stations or --min-years.")
    if args.apply and not args.station_data_dir:
        sys.exit("ERROR: --apply requires --station-data-dir (where to copy station CSVs from).")

    try:
        catalog = read_candidate_catalog(args.catalog)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"ERROR reading --catalog: {e}")

    target_descriptors = (_parse_kv_list(args.target_descriptors, "--target-descriptors")
                          if args.target_descriptors else None)
    weights = _parse_kv_list(args.weights, "--weights") if args.weights else None

    try:
        target = resolve_target(catalog, args.descriptors, target_station=args.target_station,
                                 target_descriptors=target_descriptors)
        ranking = similarity_ranking(target, catalog, args.descriptors, weights=weights)
        proposed = propose_pooling_group(ranking, n_stations=args.n_stations,
                                          min_total_years=args.min_years)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    paths = resolve_region(args.region_name, __file__)
    display_cols = ["station"] + args.descriptors + ["distance"]
    if "n_years" in ranking.columns:
        display_cols.append("n_years")

    print(f"Target: {target.name}  (descriptors: {target.descriptors})")
    print(f"Candidate pool: {len(catalog)} station(s), {len(ranking)} ranked "
          f"(excluding the target itself)")
    if ranking.attrs.get("dropped_missing_descriptors"):
        print(f"  Dropped for missing descriptor value(s): "
              f"{ranking.attrs['dropped_missing_descriptors']}")
    print()
    print("Full ranking (most similar first):")
    print(ranking[display_cols].to_string(index=False))
    print()
    print(f"Proposed pooling group ({len(proposed)} station(s)):")
    print(proposed[display_cols].to_string(index=False))
    print()
    print("This is a PROPOSAL, not a substitute for discordancy/heterogeneity screening --")
    print("run run_regional_analysis.py on the resulting group before trusting it.")

    ranking[display_cols].to_csv(paths.output_dir / "pooling_ranking_full.csv", index=False)
    proposed[display_cols].to_csv(paths.output_dir / "pooling_group_proposed.csv", index=False)
    print(f"\nWrote {paths.output_dir / 'pooling_ranking_full.csv'}")
    print(f"Wrote {paths.output_dir / 'pooling_group_proposed.csv'}")

    if args.apply:
        src_dir = Path(args.station_data_dir)
        paths.data_dir.mkdir(parents=True, exist_ok=True)
        copied, missing = [], []
        for name in proposed["station"]:
            src = src_dir / f"{name}.csv"
            if src.exists():
                shutil.copy2(src, paths.data_dir / f"{name}.csv")
                copied.append(name)
            else:
                missing.append(name)
        print(f"\nCopied {len(copied)} station CSV(s) into {paths.data_dir}")
        if missing:
            print(f"WARNING: no CSV found (looked for <station>.csv in --station-data-dir) for: "
                  f"{missing} -- add these manually.")
    else:
        print(f"\nDry run only -- nothing written to Data/. Re-run with --apply "
              f"--station-data-dir <path> to copy the proposed stations' CSVs into "
              f"{paths.data_dir}.")


if __name__ == "__main__":
    main()
