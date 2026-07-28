"""
Automatic pooling-group formation: the "region of influence" approach
(Burn, 1990) to proposing which stations to pool for a regional (index-
flood) analysis, based on similarity in catchment/site descriptors --
rather than assembling the group by hand.

This module is deliberately descriptor-SOURCE-agnostic: it doesn't know
or care whether a descriptor came from BasinATLAS/HydroSHEDS (catchment
area, mean annual precipitation, a soil/permeability index, land cover,
slope -- the natural choice for streamflow gauges, since these describe
how a drainage basin turns rain into runoff), a gridded climatology
sampled at a point plus geographic coordinates and elevation (the natural
choice for rain gauges, which don't have an upstream catchment for
basin-level attributes to mean anything), or something else entirely.
It just standardizes whatever numeric descriptor columns you give it and
ranks candidates by weighted distance to a target site. Extracting the
descriptors themselves is out of scope here -- see the module docstring
below (`read_candidate_catalog`) for the expected input schema.

Workflow:
    1. `read_candidate_catalog()`     -- load a pool of candidate stations
                                          and their descriptors
    2. `resolve_target()`             -- get the target site's own
                                          descriptor values (either an
                                          existing gauged candidate, for a
                                          leave-one-out-style check, or an
                                          ungauged site's descriptors
                                          supplied directly)
    3. `similarity_ranking()`         -- rank every OTHER candidate by
                                          weighted standardized distance
                                          to the target
    4. `propose_pooling_group()`      -- take the ranking and apply a
                                          stopping rule (top N stations,
                                          or accumulate until a minimum
                                          total station-years is reached)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd


def read_candidate_catalog(path, station_col: str = "station") -> pd.DataFrame:
    """
    Read a candidate catalog CSV: one row per candidate station, a
    `station` column (must match the station names/CSV stems used
    elsewhere, e.g. in `Data/Regional/<RegionName>/`), and one column per
    numeric descriptor -- whatever you choose to include, e.g. for
    streamflow gauges: `area_km2`, `mean_annual_precip_mm`, a soil/
    permeability index, mean basin slope, urban extent fraction; for rain
    gauges: `lat`, `lon`, `elevation_m`, `mean_annual_precip_mm` (sampled
    directly from a climate grid at the point, rather than via a basin
    polygon that has no real relationship to a rain gauge).

    An optional `n_years` column, if present, is used by
    `propose_pooling_group()`'s min-total-years stopping rule. Any other
    non-numeric columns (notes, source, etc.) are preserved but ignored
    for distance calculations.
    """
    df = pd.read_csv(path)
    if station_col not in df.columns:
        raise ValueError(f"Candidate catalog {path} has no '{station_col}' column. "
                          f"Available columns: {list(df.columns)}")
    if df[station_col].duplicated().any():
        dupes = df.loc[df[station_col].duplicated(), station_col].tolist()
        raise ValueError(f"Candidate catalog {path} has duplicate station name(s): {dupes}")
    return df


@dataclass
class Target:
    """The site a pooling group is being formed around, described by the
    same descriptor columns as the candidate catalog."""
    name: str
    descriptors: dict


def resolve_target(candidates: pd.DataFrame, descriptor_cols: Sequence[str],
                    target_station: str | None = None, target_descriptors: dict | None = None,
                    station_col: str = "station") -> Target:
    """
    Get the target site's descriptor values, either:
      - `target_station`: an existing row in `candidates` (e.g. a gauged
        station you want to check "if this were pooled, who would its
        region look like" -- a leave-one-out-style use), or
      - `target_descriptors`: a dict of descriptor values for a site NOT
        in the catalog (the normal ungauged-target case).
    Exactly one of the two must be given.
    """
    if (target_station is None) == (target_descriptors is None):
        raise ValueError("Give exactly one of target_station or target_descriptors.")
    if target_station is not None:
        row = candidates[candidates[station_col] == target_station]
        if row.empty:
            raise ValueError(f"target_station '{target_station}' not found in the candidate "
                              f"catalog. Available: {candidates[station_col].tolist()}")
        missing = [c for c in descriptor_cols if pd.isna(row.iloc[0][c])]
        if missing:
            raise ValueError(f"target_station '{target_station}' is missing a value for "
                              f"descriptor(s) {missing}.")
        return Target(name=target_station,
                       descriptors={c: float(row.iloc[0][c]) for c in descriptor_cols})
    else:
        missing = [c for c in descriptor_cols if c not in target_descriptors]
        if missing:
            raise ValueError(f"target_descriptors is missing value(s) for {missing}.")
        return Target(name="(ungauged target)",
                       descriptors={c: float(target_descriptors[c]) for c in descriptor_cols})


def similarity_ranking(target: Target, candidates: pd.DataFrame, descriptor_cols: Sequence[str],
                        weights: dict | None = None, station_col: str = "station") -> pd.DataFrame:
    """
    Rank every candidate by weighted Euclidean distance to `target` in
    STANDARDIZED descriptor space -- the "region of influence" approach
    (Burn, 1990): smaller distance = more similar = a better pooling
    candidate. Standardization (z-score: subtract mean, divide by std)
    uses the candidate pool's own mean/std for each descriptor, so
    descriptors on very different natural scales (e.g. catchment area in
    km^2 vs. a 0-1 soil index) contribute comparably; the target itself
    is projected onto that same scale rather than folded into computing
    it, so the ranking doesn't shift depending on where the target
    happens to sit relative to the pool.

    `weights`: optional {descriptor_name: weight}, applied after
    standardization (default: equal weight, 1.0, for every descriptor).
    A larger weight makes that descriptor count for more in the distance.

    If `target.name` matches a station in `candidates` (the leave-one-out
    case), that row is excluded from its own ranking.

    Returns candidates sorted by ascending distance, with a `distance`
    column added.
    """
    df = candidates[candidates[station_col] != target.name].copy()
    if df.empty:
        raise ValueError("No candidates left to rank (catalog only contained the target itself?).")

    missing_cols = [c for c in descriptor_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Candidate catalog is missing descriptor column(s): {missing_cols}")

    incomplete = df[descriptor_cols].isna().any(axis=1)
    if incomplete.any():
        dropped = df.loc[incomplete, station_col].tolist()
        df = df.loc[~incomplete].copy()
        if df.empty:
            raise ValueError("Every remaining candidate is missing at least one descriptor value.")
    else:
        dropped = []

    means = df[descriptor_cols].mean()
    stds = df[descriptor_cols].std(ddof=0)
    zero_std = stds[stds == 0].index.tolist()
    if zero_std:
        raise ValueError(f"Descriptor(s) {zero_std} have zero variance across the candidate "
                          f"pool (every candidate has the same value) and can't be standardized. "
                          f"Drop them from descriptor_cols.")

    w = pd.Series({c: (weights.get(c, 1.0) if weights else 1.0) for c in descriptor_cols})

    z_candidates = (df[descriptor_cols] - means) / stds
    z_target = pd.Series({c: (target.descriptors[c] - means[c]) / stds[c] for c in descriptor_cols})

    sq_diff = (z_candidates - z_target) ** 2 * w
    df["distance"] = np.sqrt(sq_diff.sum(axis=1))

    result = df.sort_values("distance").reset_index(drop=True)
    result.attrs["dropped_missing_descriptors"] = dropped
    result.attrs["target_name"] = target.name
    result.attrs["descriptor_cols"] = list(descriptor_cols)
    return result


def propose_pooling_group(ranking: pd.DataFrame, n_stations: int | None = None,
                           min_total_years: int | None = None, years_col: str = "n_years",
                           station_col: str = "station") -> pd.DataFrame:
    """
    Apply a stopping rule to an already-computed `similarity_ranking()`
    table and return the proposed pooling group (a prefix of the ranking,
    most-similar first).

    Exactly one of:
      - `n_stations`: take the N most similar candidates.
      - `min_total_years`: take the most similar candidates one at a time
        until their combined record length reaches at least this many
        station-years (Hosking & Wallis, 1997 suggest, as a rule of
        thumb, at least about 5T station-years of record to estimate a
        T-year event with reasonable confidence via pooling). Requires
        `years_col` to be present in `ranking` (from the candidate
        catalog's optional `n_years` column).

    Neither guarantees the result is actually homogeneous -- this is a
    proposal to try, not a substitute for running discordancy/
    heterogeneity on the resulting group (via `regional.py`) once you've
    assembled their station CSVs.
    """
    if (n_stations is None) == (min_total_years is None):
        raise ValueError("Give exactly one of n_stations or min_total_years.")

    if n_stations is not None:
        if n_stations < 1:
            raise ValueError("n_stations must be >= 1.")
        return ranking.head(n_stations).reset_index(drop=True)

    if years_col not in ranking.columns:
        raise ValueError(f"min_total_years requires a '{years_col}' column in the candidate "
                          f"catalog (record length per candidate station), which wasn't found. "
                          f"Available columns: {list(ranking.columns)}")
    cum_years = ranking[years_col].cumsum()
    n_needed = int(np.searchsorted(cum_years.to_numpy(), min_total_years) + 1)
    n_needed = min(n_needed, len(ranking))
    chosen = ranking.head(n_needed).reset_index(drop=True)
    if cum_years.iloc[min(n_needed, len(cum_years)) - 1] < min_total_years:
        import warnings
        warnings.warn(
            f"Only {cum_years.iloc[-1]:.0f} total station-years available across all "
            f"{len(ranking)} candidates, short of the requested min_total_years="
            f"{min_total_years}. Returning every available candidate.", stacklevel=2)
    return chosen
