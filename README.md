# Flood Frequency Analysis

[![Tests](https://github.com/<OWNER>/<REPO>/actions/workflows/tests.yml/badge.svg)](https://github.com/<OWNER>/<REPO>/actions/workflows/tests.yml)
<!-- Replace <OWNER>/<REPO> above with this repo's actual GitHub path once pushed, e.g. jsmith/Flood_Frequency_Analysis -->

A Python tool for flood frequency analysis — and more generally, frequency analysis of any annual-maximum series (rainfall included) — using [uv](https://docs.astral.sh/uv/) for dependency management. It has three main objectives:

1. **[Single-station analysis](#how-to-run-an-analysis)** — fits 9 candidate distributions (Normal, LogNormal 2p/3p, Gumbel, GEV, Exponential, Gamma, Pearson III, Log-Pearson III) to one station's annual-maximum series via method of moments / maximum likelihood / probability-weighted moments, ranks them by AIC/BIC/KS/Anderson-Darling, and reports a best-fit recommendation (plus an Akaike-weighted multi-model average) with bootstrap confidence intervals.
2. **[Regional (pooled) analysis](#regional-pooled-flood-frequency-analysis)** — the index-flood method (Dalrymple, 1960; Hosking & Wallis, 1997): pools several hydrologically-similar stations so the *shape* of the flood distribution is estimated from many station-years combined, with full discordancy/heterogeneity/distribution-selection diagnostics, per-station data quality checks, and Monte-Carlo confidence intervals.
3. **[Automatic pooling-group formation](#forming-a-pooling-group-automatically)** — the "region of influence" approach (Burn, 1990): given a target site and a catalog of candidate stations described by numeric descriptors (catchment area, precipitation, etc. for streamflow; coordinates, elevation, climatology for rainfall), ranks candidates by similarity and proposes which ones to pool, rather than requiring the group to be assembled by hand.

All three work for streamflow, rainfall, or any other generic annual-maximum variable (see [Beyond streamflow](#beyond-streamflow-any-annual-maximum-series)).

## What's included
```
Flood_Frequency_Analysis/
├── .github/workflows/tests.yml   # GitHub Actions: runs the pytest suite on every push/PR
├── Data/            # input data, one CSV per case: Data/<CaseName>.csv (not tracked in git)
│   ├── Templates/                # example pooling-group candidate catalogs (NOT a single-station case -- kept out
│   │                              # of Data/*.csv so run_analysis.py's picker never lists them by mistake)
│   │   ├── candidate_descriptors_discharge_TEMPLATE.csv
│   │   ├── candidate_descriptors_rainfall_TEMPLATE.csv
│   │   ├── discharge_station_data/     # matching annual-maximum series, one CSV per candidate
│   │   └── rainfall_station_data/      # matching annual-maximum series, one CSV per candidate
│   └── Regional/<RegionName>/   # station CSVs for one regional (pooled) analysis group
├── Module/
│   ├── run_analysis.py           # single-station CLI entry point
│   ├── run_regional_analysis.py  # regional (pooled) CLI entry point
│   ├── form_pooling_group.py     # propose a pooling group from a candidate descriptor catalog
│   ├── floodfreq/                # the package
│   └── tests/                    # pytest regression suite (200+ tests)
├── Output/          # generated CSVs, one subfolder per case: Output/<CaseName>/ (not tracked in git)
│   └── Regional/<RegionName>/    # regional analysis CSV outputs
├── Plot/            # generated PNG plots, one subfolder per case: Plot/<CaseName>/ (not tracked in git)
│   └── Regional/<RegionName>/    # regional analysis PNG outputs
├── pyproject.toml   # project metadata + dependencies (uv-managed)
└── .gitignore       # excludes venv, cache files, Data/Output/Plot, OS junk, etc.
```

## How to set up

```bash
cd Flood_Frequency_Analysis
uv sync
```

## How to run the tests

```bash
uv run pytest Module/tests/ -v
```

This same command also runs automatically on every push and pull request via GitHub Actions (`.github/workflows/tests.yml`), against the committed `uv.lock` for a reproducible environment. Check the Actions tab on GitHub for results, or the badge at the top of this file.

## How to run an analysis

Drop your annual-maximum series into `Data/<CaseName>.csv`, then either:

```bash
uv run python Module/run_analysis.py <CaseName>
```
or run it with no arguments to pick interactively from whatever's in `Data/`
(useful when running from an editor's "Run" button rather than a terminal):
```bash
uv run python Module/run_analysis.py
```

For example, with the included `Data/P1009.csv`:
```bash
uv run python Module/run_analysis.py P1009
```

This will:
1. Run data quality checks (stationarity/trend, outliers, basic input validation)
2. Fit all 9 distributions using their recommended method (MLE for Normal/LogNormal, PWM for the extreme-value/skewed families, MOM for Gamma) — plus, if `--regional-skew` is given, weighted-skew variants of Pearson III/LP3 as additional candidates
3. Rank them by AIC/BIC, Kolmogorov-Smirnov, and Anderson-Darling (tail-weighted) goodness-of-fit
4. Prompt you for a confidence level (e.g. `95`) unless `--confidence-level` is passed
5. Write CSVs to `Output/<CaseName>/`:
   - `descriptive_stats.csv`, `data_quality.csv` (with a plain-language `meaning` column), `goodness_of_fit.csv` (with a `fit_assessment` column)
   - `quantile_table.csv` up to T=10,000 years, `bootstrap_ci_<distribution>.csv`
   - `akaike_weights.csv`, `model_averaged_quantiles.csv` — Akaike-weighted multi-model-averaged design flood, an alternative to betting everything on one "best" distribution
   - `summary.txt` — the full recommendation, with Good/Caution/Warning labels throughout and an extrapolation-risk warning for large T
6. Write PNG plots to `Plot/<CaseName>/`:
   - `data_quality_timeseries.png` — annual maxima with trend line and any flagged outliers highlighted
   - `input_data_histogram.png` — histogram of the raw series with the recommended distribution's fitted PDF overlaid
   - `flood_frequency_curve.png` — all fitted distributions vs. observed data, extrapolated out to T=10,000 years
   - `moment_ratio_diagram.png` — L-moment ratio diagram for distribution-family diagnosis
   - `bestfit_<distribution>_ci.png` — the recommended distribution's quantile curve with a bootstrap CI band
   - `dashboard.png` — all of the above combined into a single one-page overview, plus a text panel with the key numbers

Useful flags:
```bash
uv run python Module/run_analysis.py P1009 --confidence-level 90
uv run python Module/run_analysis.py P1009 --plotting-position gringorten
uv run python Module/run_analysis.py P1009 --n-boot 2000 --xlsx-report --pdf-report
uv run python Module/run_analysis.py P1009 --regional-skew 0.0 --regional-skew-mse 0.15
uv run python Module/run_analysis.py P1009 --no-plots
```

## Beyond streamflow: any annual-maximum series

The statistics here (GEV, Gumbel, LogNormal, Pearson III, etc.) are fully generic — nothing requires the input to be discharge. To analyze, e.g., annual maximum rainfall depth instead of flood flow, relabel the output with:
```bash
uv run python Module/run_analysis.py RainStation --variable-name "Rainfall depth" --units mm --short-name Rainfall
```
This changes every plot axis, chart title, `summary.txt` header, and Excel/PDF report title to match (e.g. "Rainfall frequency curve", "RAINFALL FREQUENCY ANALYSIS", axis label "Rainfall depth (mm)"). Leaving these flags unset reproduces the exact original "Flood magnitude" / "Flood Frequency Analysis" wording — this is purely additive, nothing changes unless you opt in.

**One caveat that does NOT get handled by these flags:** `--regional-skew` (Bulletin 17B weighted skew) uses MSE coefficients empirically calibrated from streamflow data specifically. Avoid it for non-streamflow variables regardless of `--short-name`.

`--pdf-report` bundles the text summary, the dashboard, and every individual full-size plot into one `Output/<CaseName>/<CaseName>_report.pdf` — the single file to actually hand to a colleague, instead of the separate CSVs/PNGs.

## Regional (pooled) flood frequency analysis

For a group of hydrologically-similar stations, you can get a better estimate of the flood distribution's *shape* by pooling them — the classical index-flood method (Dalrymple, 1960), using the L-moment-based regional statistics of Hosking & Wallis (1997) (see also Rao & Hamed, 2000, Ch. 9). Each station keeps its own scale (its "index flood," the station's own mean), while the shape of the growth curve is estimated from all the pooled station-years combined.

Drop one CSV per station into `Data/Regional/<RegionName>/` (same `year,Q`-style format as single-station cases), then:
```bash
uv run python Module/run_regional_analysis.py <RegionName>
```
or with no arguments to pick interactively from whatever region folders exist under `Data/Regional/`. For example, with your own `Data/Regional/Template/` populated with station CSVs:
```bash
uv run python Module/run_regional_analysis.py Template
```

This will:
1. Check each station's data quality BEFORE pooling — the same Mann-Kendall stationarity test, Grubbs' outlier test, and basic input validation the single-station tool runs, applied to each station's own record (a trending or outlier-contaminated station can quietly bias the pooled growth curve). Station CSVs with a detectable year column (auto-detected the same way as single-station cases, or via `--year-col`) get the trend test keyed to real calendar time, Sen's slope reported in per-calendar-year units, and a check for missing years within a station's nominal span; stations without one fall back to record order, same as the single-station tool's own fallback.
2. Compute each station's sample L-moments/ratios (L-CV, L-skewness, L-kurtosis)
3. Screen for discordant stations — a discordancy measure D_i flags any station statistically inconsistent with the rest of the group, against Hosking & Wallis's tabulated critical values
4. Test whether the group is homogeneous enough to pool — the heterogeneity measure H1 compares the observed dispersion of at-site L-CV against Monte-Carlo-simulated synthetic homogeneous regions (H1 < 1 acceptably homogeneous, 1–2 possibly heterogeneous, ≥ 2 definitely heterogeneous)
5. Pick a regional distribution family — a Z-statistic compares the regional-average L-kurtosis against 5 candidates (GLO, GEV, GNO, PE3, GPA); `--family` overrides this to force a specific one
6. Fit the regional growth curve to the pooled, record-length-weighted L-moments, then compute each station's design flood (growth curve × that station's own index flood) out to the requested return periods
7. Compute a Monte-Carlo confidence interval for every station's design flood (skip with `--no-ci`) — combines regional growth-curve shape uncertainty (simulated synthetic homogeneous regions, same machinery as step 4) with each station's own index-flood sampling uncertainty (a bootstrap of that station's own record); the regional counterpart of the single-station tool's `bootstrap_ci_<distribution>.csv`
8. Write CSVs to `Output/Regional/<RegionName>/`: `station_lmoments.csv`, `data_quality.csv`, `discordancy.csv`, `heterogeneity.csv`, `zstatistics.csv`, `quantile_table.csv` (one column per station), `growth_curve_quantiles_ci.csv` (long format: station, T, design flood, CI bounds — unless `--no-ci`), and `summary.txt` with the full plain-language walkthrough (same RUN PROVENANCE header convention as the single-station tool)
9. Write PNG plots to `Plot/Regional/<RegionName>/`:
   - `regional_moment_ratio_diagram.png` — every station's (L-skewness, L-kurtosis) point, the regional average, and the GLO/GEV/GNO/PE3/GPA reference curves, with discordant stations outlined in red — the visual homogeneity check
   - `regional_growth_curve.png` — the dimensionless regional growth curve (a "Dalrymple plot"), with every station's own data rescaled by its own index flood so all stations land on one common curve, plus a light-grey diagnostic overlay showing all stations' growth factors pooled and re-ranked together as one combined record (subtle by design — see `regional_pooled_vs_stations.png` below for a much clearer side-by-side version of the same comparison)
   - `regional_pooled_vs_stations.png` — two panels at equal visual weight: left is the per-station view (what the curve is actually fit to); right is every station's data pooled into one array and re-ranked as if it were a single combined record (`sum(n_i)` years), colored uniformly since station identity is gone once concatenated — the clearest way to see how much further out along the T axis naive concatenation appears to reach, and why the index-flood method doesn't do that
   - `regional_discordancy.png` — bar chart of each station's D_i against the group's critical value
   - `regional_station_series.png` — small multiples of each station's raw annual-maximum series (plotted against calendar year when available, observation order otherwise), for a quick visual screen before trusting the pooling
   - `station_<StationName>_design_flood.png` — one per station: that station's own observed data against its regional design-flood curve (growth curve × that station's index flood) — the file to actually hand to a client for that site
   - `regional_dashboard.png` — growth curve, L-moment ratio diagram, discordancy chart, and a text summary panel combined into one page

**Note on how stations are actually pooled:** the growth curve is fit to a record-length-weighted *average of each station's own L-moment ratios* (`t_R = Σ(n_i·t_i)/Σ(n_i)`, etc.), computed after each station's data is rescaled by its own index flood — not by concatenating every station's raw (or rescaled) values into one long series and computing L-moments once on that. Naive concatenation would let whichever station has the most years dominate the shape estimate; the weighted-average-of-ratios approach lets each station contribute in proportion to its own record length while still keeping its own shape estimate distinct, which is also what makes the discordancy/heterogeneity diagnostics meaningful station-by-station. The grey "pooled rank" overlay on `regional_growth_curve.png`, and the much clearer `regional_pooled_vs_stations.png`, both show what naive concatenation would look like, for comparison only.

Useful flags:
```bash
uv run python Module/run_regional_analysis.py Template --n-sim 1000 --seed 42
uv run python Module/run_regional_analysis.py Template --family gev
uv run python Module/run_regional_analysis.py Template --candidates gev gno pe3
uv run python Module/run_regional_analysis.py Template --return-periods 10 50 100 500
uv run python Module/run_regional_analysis.py Template --confidence-level 90
uv run python Module/run_regional_analysis.py Template --no-ci
```

Needs at least 4 stations (discordancy requires inverting a 3×3 covariance matrix). Regional analysis is separate from the single-station tool above and doesn't currently read per-station `.toml` settings files.

### Forming a pooling group automatically

**Do you actually need this?** Only if you're choosing among a *larger pool* of candidate stations and want a data-driven suggestion for which ones to pool. If you already know which stations belong together — the common case if you're working from a short, hand-picked list — skip this entirely and just build `Data/Regional/<RegionName>/` directly (copy `Data/Regional/Template/`, swap in real data). There's nothing here that's required before running `run_regional_analysis.py`.

If you do have a larger candidate pool, the four-step recipe:

1. **Duplicate the descriptor catalog** — one row per candidate station:
   ```bash
   cp Data/Templates/candidate_descriptors_discharge_TEMPLATE.csv Data/Templates/my_candidates.csv
   ```
   Replace the example rows with your real candidates' real descriptor values.

2. **Duplicate the matching station-data folder** — each candidate's own annual-maximum series:
   ```bash
   cp -r Data/Templates/discharge_station_data Data/Templates/my_station_data
   ```
   Replace the example CSVs with your real stations' `year,Q` files. **The filename must exactly match the `station` value in your catalog** — a catalog row with `station=GAUGE_014` needs a file named `GAUGE_014.csv`.

3. **Run it** against your own files:
   ```bash
   uv run python Module/form_pooling_group.py \
       --catalog Data/Templates/my_candidates.csv \
       --descriptors area_km2 mean_annual_precip_mm \
       --target-station GAUGE_014 --n-stations 6 --region-name MyNewRegion \
       --station-data-dir Data/Templates/my_station_data --apply
   ```

4. **Run the actual analysis** on the result:
   ```bash
   uv run python Module/run_regional_analysis.py MyNewRegion
   ```

The rest of this section covers the details behind that recipe: descriptor sources, the ranking math, and every flag.

Rather than hand-picking which stations to pool, `form_pooling_group.py` implements the "region of influence" approach (Burn, 1990): given a target site and a catalog of candidate stations described by numeric descriptors, it ranks candidates by similarity and proposes a group.

This is deliberately **descriptor-source-agnostic** — it doesn't extract descriptors itself, it just consumes whatever numeric columns you give it in a candidate catalog CSV (`station` column + descriptor columns + an optional `n_years` column). Which descriptors make sense depends on what you're pooling:

- **Streamflow gauges**: catchment area, mean annual precipitation, a soil/permeability index, channel slope, urban extent — i.e. attributes of the drainage basin. A natural source is [HydroSHEDS/BasinATLAS](https://www.hydrosheds.org/), which provides these at the sub-basin level (snap each gauge's coordinates to its containing sub-basin to look them up).
- **Rain gauges** (or any other generic annual-maximum variable): geographic coordinates, elevation, and a point climatology (e.g. mean annual precipitation sampled directly from a gridded product like WorldClim/CHELSA at the station's coordinates) — basin-level attributes like soil permeability or channel slope don't mean anything for a point measurement with no upstream catchment.

Producing the candidate catalog CSV (i.e. actually looking up descriptor values for your stations) is outside this script's scope — bring your own GIS workflow. Two example catalogs are bundled so you can try the CLI immediately, no editing required: `Data/Templates/candidate_descriptors_discharge_TEMPLATE.csv` (10 synthetic gauges: `area_km2`, `mean_annual_precip_mm`, `bfihost`, `mean_slope_m_per_km`, `urban_extent_frac`, `n_years`) and `Data/Templates/candidate_descriptors_rainfall_TEMPLATE.csv` (10 synthetic rain gauges: `lat`, `lon`, `elevation_m`, `mean_annual_precip_mm`, `n_years`). Column names in both are illustrative — rename them to match your own extraction, or just pass whatever names your own catalog actually uses via `--descriptors`.

Each catalog has a matching folder of synthetic annual-maximum series, one CSV per candidate, with each station's record length matching its catalog row's `n_years` and its magnitude scaled (loosely, illustratively) with its own `area_km2`/`mean_annual_precip_mm` — enough to run the full pipeline (`--apply` → `run_regional_analysis.py`) end to end using only bundled content: `Data/Templates/discharge_station_data/` (columns `year,Q`) and `Data/Templates/rainfall_station_data/` (columns `year,Rainfall_mm`).

```bash
# Try it right now with the bundled discharge template (dry run)
uv run python Module/form_pooling_group.py \
    --catalog Data/Templates/candidate_descriptors_discharge_TEMPLATE.csv \
    --descriptors area_km2 mean_annual_precip_mm bfihost mean_slope_m_per_km urban_extent_frac \
    --target-station GAUGE_001 --n-stations 4 --region-name MyRegion

# ...or the bundled rainfall template
uv run python Module/form_pooling_group.py \
    --catalog Data/Templates/candidate_descriptors_rainfall_TEMPLATE.csv \
    --descriptors lat lon elevation_m mean_annual_precip_mm \
    --target-station RAINGAUGE_001 --n-stations 4 --region-name MyRegion

# Full pipeline using ONLY bundled content: form the group, copy the station
# CSVs, then actually run the regional analysis on the result
uv run python Module/form_pooling_group.py \
    --catalog Data/Templates/candidate_descriptors_discharge_TEMPLATE.csv \
    --descriptors area_km2 mean_annual_precip_mm bfihost mean_slope_m_per_km urban_extent_frac \
    --target-station GAUGE_001 --n-stations 5 --region-name MyRegion \
    --station-data-dir Data/Templates/discharge_station_data --apply
uv run python Module/run_regional_analysis.py MyRegion

# Dry run against your own catalog
uv run python Module/form_pooling_group.py \
    --catalog candidates.csv \
    --descriptors area_km2 mean_annual_precip_mm bfihost \
    --target-station GAUGE_042 --n-stations 8 --region-name MyRegion

# Ungauged target site, stop once 250 station-years are pooled
uv run python Module/form_pooling_group.py \
    --catalog candidates.csv \
    --descriptors area_km2 mean_annual_precip_mm \
    --target-descriptors area_km2=180 mean_annual_precip_mm=1050 \
    --min-years 250 --region-name MyRegion

# Actually copy the proposed stations' CSVs into Data/Regional/MyRegion/
uv run python Module/form_pooling_group.py \
    --catalog candidates.csv \
    --descriptors area_km2 mean_annual_precip_mm \
    --target-station GAUGE_042 --n-stations 8 --region-name MyRegion \
    --station-data-dir /path/to/all_candidate_series/ --apply
```

This ranks every other candidate by weighted Euclidean distance to the target in standardized descriptor space (z-scored using the candidate pool's own mean/std, so descriptors on very different scales — e.g. area in km² vs. a 0-1 soil index — contribute comparably), then applies one of two stopping rules: `--n-stations` (take the N most similar) or `--min-years` (accumulate the most similar candidates until their combined record length reaches a target, per Hosking & Wallis's rule-of-thumb of roughly 5× the design return period in station-years). Writes `pooling_ranking_full.csv` and `pooling_group_proposed.csv` to `Output/Regional/<RegionName>/`; with `--apply`, also copies the proposed stations' own CSVs (from `--station-data-dir`) into `Data/Regional/<RegionName>/`, ready for `run_regional_analysis.py`.

**This is a proposal step, not a substitute for discordancy/heterogeneity screening** — always run `run_regional_analysis.py` on the resulting group afterward; a statistically similar-looking descriptor set doesn't guarantee a homogeneous region.

## Per-case settings file (optional)

Instead of retyping flags every time, save a case's settings in `Data/<CaseName>.toml`:
```toml
confidence_level = 90
regional_skew = 0.0
regional_skew_mse = 0.15
n_boot = 2000
pdf_report = true
```
Then `uv run python Module/run_analysis.py P1009` picks these up automatically — no flags needed. Any CLI flag you *do* pass still overrides the file for that one run (e.g. `--confidence-level 95` overrides the `90` above without touching the saved file). Recognized keys mirror the CLI flag names with underscores instead of dashes: `value_col`, `year_col`, `variable_name`, `units`, `short_name`, `plotting_position`, `descriptive_plotting_position`, `n_boot`, `confidence_level`, `regional_skew`, `regional_skew_mse`, `no_plots`, `xlsx_report`, `pdf_report`.

A case with no `.toml` file behaves exactly as before — this is entirely optional.

**Note on git tracking:** `.gitignore` currently excludes everything under `Data/` (including `.toml` files), matching your instruction to keep `Data/`'s contents out of version control. If you'd rather have these small settings files tracked in git (they're not bulky/sensitive like the CSV data, and arguably worth sharing/preserving — "this is how station X should be analyzed"), that needs a small `.gitignore` change; let me know if you want that.

Run `uv run python Module/run_analysis.py --help` for the full list.

## Default dependencies

- numpy
- pandas
- scipy
- matplotlib
- openpyxl
- lmoments3 (probability-weighted-moment / L-moment fitting)

Add or remove packages as needed:
```bash
uv add <package>
uv remove <package>
```

## Notes

- `[tool.uv] package = false` in `pyproject.toml` marks this as a scripts project rather than an installable package — required so `uv sync`/`uv add` don't try (and fail) to build a wheel.
- `Data/`, `Output/`, and `Plot/` are excluded from git by default, since they typically hold large, regenerated, or locally-specific files — with explicit exceptions for `Data/Template.csv`, `Data/Templates/candidate_descriptors_{discharge,rainfall}_TEMPLATE.csv`, `Data/Templates/{discharge,rainfall}_station_data/*.csv`, and `Data/Regional/Template/*.csv` so a fresh clone has something to try immediately. Adjust `.gitignore` per-project if you want anything else tracked (e.g. your own real station data, which most people will want to keep local/private rather than committed).
- Reading legacy `.xls` (not `.xlsx`) input directly would additionally need `xlrd` (`uv add xlrd`) — not included by default since the expected input is `Data/<CaseName>.csv`.
- The "plotting position" formula (Weibull, Hazen, Cunnane, Gringorten, Hosking, Blom) affects where empirical points are drawn on probability plots and the descriptive-statistics summary — it does **not** change PWM-fitted parameters, which always use the standard unbiased L-moment estimator. See `summary.txt` for the full explanation, generated with each run.
- Quantiles are reported out to T=10,000 years by default, but treat anything much beyond ~2-3x your record length as increasingly uncertain (`summary.txt` flags this per-row and adds an explicit warning when T exceeds 10x the record length) — that far out, the design flood depends more on which distribution's tail you trust than on the data itself.
- Keep `Data/` free of stray/example CSVs you don't intend to analyze. The interactive picker (when running with no case name) shows row count and last-modified time for each file to help tell them apart, but it's easy to select the wrong one by accident, and every result and plot title is derived from that file's name — including the bundled `Data/Template.csv` demo, which will show up in that picker alongside your own real cases.
