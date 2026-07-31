# Flood Frequency Analysis

[![Tests](https://github.com/<OWNER>/<REPO>/actions/workflows/tests.yml/badge.svg)](https://github.com/<OWNER>/<REPO>/actions/workflows/tests.yml)
<!-- Replace <OWNER>/<REPO> above with this repo's actual GitHub path once pushed, e.g. jsmith/Flood_Frequency_Analysis -->

A Python tool for flood frequency analysis — and more generally, frequency analysis of any annual-maximum series (rainfall included) — using [uv](https://docs.astral.sh/uv/) for dependency management. It has four main objectives:

1. **[Single-station analysis](#how-to-run-an-analysis)** — fits 9 candidate distributions (Normal, LogNormal 2p/3p, Gumbel, GEV, Exponential, Gamma, Pearson III, Log-Pearson III) to one station's annual-maximum series via method of moments / maximum likelihood / probability-weighted moments, ranks them by AIC/BIC/KS/Anderson-Darling, and reports a best-fit recommendation (plus an Akaike-weighted multi-model average) with bootstrap confidence intervals.
2. **[Regional (pooled) analysis](#regional-pooled-flood-frequency-analysis)** — the index-flood method (Dalrymple, 1960; Hosking & Wallis, 1997): pools several hydrologically-similar stations so the *shape* of the flood distribution is estimated from many station-years combined, with full discordancy/heterogeneity/distribution-selection diagnostics, per-station data quality checks, and Monte-Carlo confidence intervals.
3. **[Automatic pooling-group formation](#forming-a-pooling-group-automatically)** — the "region of influence" approach (Burn, 1990): given a target site and a catalog of candidate stations described by numeric descriptors (catchment area, precipitation, etc. for streamflow; coordinates, elevation, climatology for rainfall), ranks candidates by similarity and proposes which ones to pool, rather than requiring the group to be assembled by hand.
4. **[Climate-informed adjustment (CIFAM)](#climate-informed-adjustment-cifam)** — the rapid Climate-Informed Flood Assessment Methodology (Grijsen & Lino, ICOLD 2026): adjusts a single-station fit for projected climate change by shifting the distribution's first two moments and widening its confidence intervals to combine sampling uncertainty with climate-change uncertainty, for Gumbel / Log-Normal / Pearson III (closed form) or any distribution (Monte Carlo).

All four work for streamflow, rainfall, or any other generic annual-maximum variable (see [Beyond streamflow](#beyond-streamflow-any-annual-maximum-series)).

## What's included
```
Flood_Frequency_Analysis/
├── .github/workflows/tests.yml   # GitHub Actions: runs the pytest suite on every push/PR
├── Data/            # input data, case-first: Data/<CaseName>/ (only Data/Templates/ is tracked)
│   ├── Templates/                # ALL shipped example/template files live here (the only tracked
│   │   │                         # part of Data/); a real case under Data/<CaseName>/ is never tracked
│   │   ├── Template.csv           # single-station baseline series template
│   │   ├── Climate_Adjustment/    # climate-adjustment scenario template(s)
│   │   └── Regional/              # regional + pooling-group templates:
│   │       ├── candidate_descriptors_discharge_TEMPLATE.csv   # pooling-group candidate catalogs
│   │       ├── candidate_descriptors_rainfall_TEMPLATE.csv
│   │       ├── discharge_station_data/     # matching annual-maximum series, one CSV per candidate
│   │       ├── rainfall_station_data/      # matching annual-maximum series, one CSV per candidate
│   │       └── Template/                    # regional-analysis station CSVs (STN_*.csv)
│   ├── <CaseName>/              # one folder per case (not tracked)
│   │   ├── <CaseName>.csv        #   baseline annual-maximum series (named after the folder)
│   │   ├── <CaseName>.toml       #   optional per-case settings
│   │   └── climate_adjustment/   #   climate-informed (CIFAM) inputs, one CSV per scenario
│   │       ├── rcp45.csv
│   │       └── rcp85.csv
│   └── Regional/<RegionName>/   # station CSVs for one regional (pooled) analysis group (not tracked)
├── Module/
│   ├── run_analysis.py           # single-station CLI entry point
│   ├── run_regional_analysis.py  # regional (pooled) CLI entry point
│   ├── form_pooling_group.py     # propose a pooling group from a candidate descriptor catalog
│   ├── run_climate_adjustment.py # climate-informed (CIFAM) adjustment CLI entry point
│   ├── floodfreq/                # the package
│   └── tests/                    # pytest regression suite (200+ tests)
├── Output/          # generated CSVs, one subfolder per case: Output/<CaseName>/ (not tracked in git)
│   ├── <CaseName>/Climate_Adjustment/<scenario>/   # climate-adjustment CSV outputs (per scenario)
│   └── Regional/<RegionName>/    # regional analysis CSV outputs
├── Plot/            # generated PNG plots, one subfolder per case: Plot/<CaseName>/ (not tracked in git)
│   ├── <CaseName>/Climate_Adjustment/<scenario>/   # climate-adjustment PNG outputs (per scenario)
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

Drop your annual-maximum series into `Data/<CaseName>/<CaseName>.csv`, then either:

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
or with no arguments to pick interactively from whatever region folders exist under `Data/Regional/`. For example, with a `Data/Regional/<RegionName>/` folder populated with station CSVs (you can copy the shipped regional example from `Data/Templates/` to start):
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

**Do you actually need this?** Only if you're choosing among a *larger pool* of candidate stations and want a data-driven suggestion for which ones to pool. If you already know which stations belong together — the common case if you're working from a short, hand-picked list — skip this entirely and just build `Data/Regional/<RegionName>/` directly (copy the shipped regional example from `Data/Templates/`, swap in real data). There's nothing here that's required before running `run_regional_analysis.py`.

If you do have a larger candidate pool, the four-step recipe:

1. **Duplicate the descriptor catalog** — one row per candidate station:
   ```bash
   cp Data/Templates/Regional/candidate_descriptors_discharge_TEMPLATE.csv Data/Templates/Regional/my_candidates.csv
   ```
   Replace the example rows with your real candidates' real descriptor values.

2. **Duplicate the matching station-data folder** — each candidate's own annual-maximum series:
   ```bash
   cp -r Data/Templates/Regional/discharge_station_data Data/Templates/Regional/my_station_data
   ```
   Replace the example CSVs with your real stations' `year,Q` files. **The filename must exactly match the `station` value in your catalog** — a catalog row with `station=GAUGE_014` needs a file named `GAUGE_014.csv`.

3. **Run it** against your own files:
   ```bash
   uv run python Module/form_pooling_group.py \
       --catalog Data/Templates/Regional/my_candidates.csv \
       --descriptors area_km2 mean_annual_precip_mm \
       --target-station GAUGE_014 --n-stations 6 --region-name MyNewRegion \
       --station-data-dir Data/Templates/Regional/my_station_data --apply
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

Producing the candidate catalog CSV (i.e. actually looking up descriptor values for your stations) is outside this script's scope — bring your own GIS workflow. Two example catalogs are bundled so you can try the CLI immediately, no editing required: `Data/Templates/Regional/candidate_descriptors_discharge_TEMPLATE.csv` (10 synthetic gauges: `area_km2`, `mean_annual_precip_mm`, `bfihost`, `mean_slope_m_per_km`, `urban_extent_frac`, `n_years`) and `Data/Templates/Regional/candidate_descriptors_rainfall_TEMPLATE.csv` (10 synthetic rain gauges: `lat`, `lon`, `elevation_m`, `mean_annual_precip_mm`, `n_years`). Column names in both are illustrative — rename them to match your own extraction, or just pass whatever names your own catalog actually uses via `--descriptors`.

Each catalog has a matching folder of synthetic annual-maximum series, one CSV per candidate, with each station's record length matching its catalog row's `n_years` and its magnitude scaled (loosely, illustratively) with its own `area_km2`/`mean_annual_precip_mm` — enough to run the full pipeline (`--apply` → `run_regional_analysis.py`) end to end using only bundled content: `Data/Templates/Regional/discharge_station_data/` (columns `year,Q`) and `Data/Templates/Regional/rainfall_station_data/` (columns `year,Rainfall_mm`).

```bash
# Try it right now with the bundled discharge template (dry run)
uv run python Module/form_pooling_group.py \
    --catalog Data/Templates/Regional/candidate_descriptors_discharge_TEMPLATE.csv \
    --descriptors area_km2 mean_annual_precip_mm bfihost mean_slope_m_per_km urban_extent_frac \
    --target-station GAUGE_001 --n-stations 4 --region-name MyRegion

# ...or the bundled rainfall template
uv run python Module/form_pooling_group.py \
    --catalog Data/Templates/Regional/candidate_descriptors_rainfall_TEMPLATE.csv \
    --descriptors lat lon elevation_m mean_annual_precip_mm \
    --target-station RAINGAUGE_001 --n-stations 4 --region-name MyRegion

# Full pipeline using ONLY bundled content: form the group, copy the station
# CSVs, then actually run the regional analysis on the result
uv run python Module/form_pooling_group.py \
    --catalog Data/Templates/Regional/candidate_descriptors_discharge_TEMPLATE.csv \
    --descriptors area_km2 mean_annual_precip_mm bfihost mean_slope_m_per_km urban_extent_frac \
    --target-station GAUGE_001 --n-stations 5 --region-name MyRegion \
    --station-data-dir Data/Templates/Regional/discharge_station_data --apply
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

## Climate-informed adjustment (CIFAM)

The `floodfreq.climate_adjustment` module implements the rapid Climate-Informed Flood Assessment Methodology (CIFAM) of Grijsen & Lino (ICOLD 2026): it takes a baseline single-station fit and produces climate-adjusted flood quantiles whose confidence intervals combine the usual sampling uncertainty with a second, independent climate-change uncertainty. The core assumption is that climate change acts on a distribution only through its first two moments (the mean and the inter-annual standard deviation); skewness is either held fixed (Gumbel) or moves only as a consequence of the mean and standard deviation shifting (Log-Normal, Pearson III).

You supply four numbers — projected changes in **flow space**, as decimals — that the tool does *not* derive for you (the same boundary as pooling-group formation not doing your GIS work): `δ1` = projected fractional increase in the mean of annual-maximum flows, `δ2` = its standard deviation across the climate-model ensemble, and `τ1`/`τ2` = the same pair for the inter-annual standard deviation. (If your projections are in precipitation space, translate them to flow space with a precipitation elasticity first; e.g. in the paper's Lom Pangar case a 20% ± 12% mean-precipitation change times an elasticity of 1.5 gives δ1 = 30%, δ2 = 18%.)

### Setting up a climate-adjustment case

The climate numbers live inside the case folder, one file per scenario, alongside the reused baseline series:

```
Data/<CaseName>/<CaseName>.csv                       # baseline annual-maximum series (the SAME
                                                     #   file the single-station tool uses)
Data/<CaseName>/climate_adjustment/<scenario>.csv    # climate numbers (+ options), one per scenario
Output/<CaseName>/Climate_Adjustment/<scenario>/     # CSV outputs (per scenario)
Plot/<CaseName>/Climate_Adjustment/<scenario>/       # PNG plot (per scenario)
```

Everything is **case-first** — a case owns one `Data/<CaseName>/` folder holding its baseline series, its optional `.toml` settings, and a `climate_adjustment/` subfolder — mirroring the case-first `Output/<CaseName>/` and `Plot/<CaseName>/` trees `run_analysis.py` already uses. Each climate **scenario** (e.g. `rcp45`, `rcp85`, `wet`, `dry`) is its own CSV, and its outputs land in a scenario-named subfolder so scenarios never overwrite each other. (Regional analysis stays `Data/Regional/<RegionName>/` — a region isn't a single case, whereas a climate adjustment is an extension of one.)

To set up a scenario, **copy a climate-adjustment template from `Data/Templates/`** into `Data/<CaseName>/climate_adjustment/<scenario>.csv` (using the same `<CaseName>` as your existing `Data/<CaseName>/<CaseName>.csv` baseline) and edit the values. The template is a self-documenting long-format table:

```
parameter,value,units,description
delta1,0.30,fraction,Projected increase in the MEAN of annual-maximum flows (flow-space; 0.30 = +30%)
delta2,0.18,fraction,Ensemble standard deviation of delta1
tau1,0.10,fraction,Projected increase in the inter-annual STANDARD DEVIATION of annual-maximum flows
tau2,0.18,fraction,Ensemble standard deviation of tau1
confidence_level,95,percent,Two-sided confidence level for the combined interval
return_periods,"2,5,10,...,100000",years,Return periods to evaluate (comma-separated)
distribution,gumbel,name,Closed-form: gumbel | lognormal2 | pearson3
pmf,,m3/s,Optional PMF reference value (leave blank if none)
```

Only `delta1/delta2/tau1/tau2` are required; the rest fall back to sensible defaults. **Values are fractions, not percents** (enter `0.30`, not `30` — the loader rejects magnitudes above 1.5 to catch that mistake). The shipped template is pre-filled with the paper's Lom Pangar numbers as a worked example.

### Running the adjustment

Run it like the other tools, passing the case name and the scenario (you'll be prompted for either if omitted):

```bash
uv run python Module/run_climate_adjustment.py MyCase rcp85
```

This reads `Data/MyCase/MyCase.csv` (baseline series) and `Data/MyCase/climate_adjustment/rcp85.csv` (the four climate numbers), then writes:

- `Output/MyCase/Climate_Adjustment/rcp85/climate_adjustment_table.csv` — per return period: baseline point + CI (sampling only) and climate point + combined CI,
- `Output/MyCase/Climate_Adjustment/rcp85/summary.txt` — a provenance-stamped summary with the inputs and the table,
- `Plot/MyCase/Climate_Adjustment/rcp85/climate_adjustment.png` — a Figure-4-style plot (baseline vs climate central lines with both confidence bands, on the Gumbel reduced-variate axis, with an optional PMF line).

Common options (all override the values in the inputs file, so a one-off variation needs no CSV edit):

```bash
uv run python Module/run_climate_adjustment.py MyCase rcp85 --distribution lognormal2
uv run python Module/run_climate_adjustment.py MyCase rcp85 --method monte-carlo --n-sim 50000
uv run python Module/run_climate_adjustment.py MyCase rcp85 --confidence-level 90 --no-plot
uv run python Module/run_climate_adjustment.py MyCase rcp85 --delta1 0.30 --delta2 0.18 --tau1 0.10 --tau2 0.18
uv run python Module/run_climate_adjustment.py MyCase rcp85 --variable-name "Peak inflow" --units m3/s --short-name Inflow
```

`--method closed-form` (default) covers Gumbel, Log-Normal (2p) and Pearson III. `--method monte-carlo` additionally supports GEV, Log-Pearson III, Normal, Gamma (2p) and Exponential — for the three-parameter skewed families (GEV, Log-Pearson III) the shape/skew is held fixed at its baseline-fitted value and only location and scale are shifted, consistent with the CIFAM assumption that climate change acts on the first two moments.

Note that the Log-Normal closed form uses an exact delta-method variance derived for this implementation rather than the paper's published Log-Normal variance term, which overestimates the true (Monte-Carlo) variance by ~45% at realistic coefficients of variation; see the module docstring and CHANGELOG [0.6.0] for details and for the Lom Pangar validation.

### Calling it from Python

The same computation is available as a library call if you'd rather script it:

```python
from floodfreq.io_utils import resolve_climate_case, load_climate_inputs, read_series
from floodfreq.climate_adjustment import (
    climate_adjusted_quantiles,      # closed form: Gumbel, Log-Normal (2p), Pearson III
    mc_climate_adjusted_quantiles,   # Monte Carlo: distribution-agnostic
)

paths = resolve_climate_case("MyCase", "rcp85", __file__)   # case + scenario; from a script in Module/
Q, _ = read_series(paths.baseline_csv)             # NB: read_series returns (values, years)
ci = load_climate_inputs(paths.climate_csv)        # the four numbers + options

result = climate_adjusted_quantiles(
    Q, ci["distribution"], ci["return_periods"],
    ci["delta1"], ci["delta2"], ci["tau1"], ci["tau2"],
    confidence_level=ci["confidence_level"],
)
result.to_frame()   # T, baseline point + CI, climate point + combined CI
```

The Monte Carlo path covers GEV, Log-Pearson III, Normal, Gamma (2p) and Exponential in addition to the three closed-form families; extending it to the remaining candidate distributions (e.g. 3-parameter log-normal) is the residual planned work.

## Per-case settings file (optional)

Instead of retyping flags every time, save a case's settings in `Data/<CaseName>/<CaseName>.toml`. A minimal file is just the settings you care about — for example, "always use a 95% CI and always write the PDF report":
```toml
confidence_level = 95
pdf_report = true
```
Then `uv run python Module/run_analysis.py P1009` picks these up automatically — no flags needed. A case with no `.toml` file, or a file that omits a key, simply uses the defaults; nothing here is required. Any CLI flag you *do* pass still overrides the file for that one run (e.g. `--confidence-level 90` overrides the `95` above without touching the saved file).

**Save it as plain text.** On macOS, TextEdit defaults to Rich Text (RTF) and will silently wrap the file in formatting markup, which fails to parse. Either use `Format → Make Plain Text` before saving, or create it from the terminal:
```bash
printf 'confidence_level = 95\npdf_report = true\n' > Data/P1009/P1009.toml
```

All keys are optional and mirror the CLI flag names with underscores instead of dashes. The most useful for a single-station run:

| key | default | what it does |
|-----|---------|--------------|
| `confidence_level` | 95 | CI level (%) for the bootstrap intervals |
| `pdf_report` | false | also write the bundled `<CaseName>_report.pdf` |
| `n_boot` | 1000 | bootstrap resamples for the CIs (higher = smoother bands, slower) |
| `no_plots` | false | skip PNG plot generation |
| `xlsx_report` | false | also write an `.xlsx` report |
| `value_col` / `year_col` | auto | column names, if auto-detection picks wrong |
| `variable_name` / `units` / `short_name` | flow defaults | axis labels and titles (e.g. for rainfall) |
| `plotting_position` / `descriptive_plotting_position` | — | plotting-position formula for probability plots |
| `regional_skew` / `regional_skew_mse` | — | weighted-skew inputs (specialized Log-Pearson III workflows) |

You almost certainly only need the first two or three. `n_boot` is worth knowing: it's how many times the tool resamples your record (with replacement), refits the distribution, and reads off the quantile spread to build the confidence interval — the default of 1000 is fine for a single-station report, and you'd only raise it (e.g. to 2000) if you want slightly smoother CI bands and don't mind a slower run.

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
- `Data/`, `Output/`, and `Plot/` are excluded from git by default, since they typically hold large, regenerated, or locally-specific files — with a single exception subtree, `Data/Templates/`, which holds all the shipped example/template files so a fresh clone has something to try immediately. A real case (`Data/<CaseName>/…`) and regional data (`Data/Regional/<RegionName>/…`) are never tracked. Adjust `.gitignore` per-project if you want anything else tracked (e.g. your own real station data, which most people will want to keep local/private rather than committed).
- Reading legacy `.xls` (not `.xlsx`) input directly would additionally need `xlrd` (`uv add xlrd`) — not included by default since the expected input is `Data/<CaseName>/<CaseName>.csv`.
- The "plotting position" formula (Weibull, Hazen, Cunnane, Gringorten, Hosking, Blom) affects where empirical points are drawn on probability plots and the descriptive-statistics summary — it does **not** change PWM-fitted parameters, which always use the standard unbiased L-moment estimator. See `summary.txt` for the full explanation, generated with each run.
- Quantiles are reported out to T=10,000 years by default, but treat anything much beyond ~2-3x your record length as increasingly uncertain (`summary.txt` flags this per-row and adds an explicit warning when T exceeds 10x the record length) — that far out, the design flood depends more on which distribution's tail you trust than on the data itself.
- Keep `Data/` tidy: each case is its own folder `Data/<CaseName>/` with a `<CaseName>.csv` baseline inside. The interactive picker (when running with no case name) lists those case folders — with row count and last-modified time — so it's easy to see what's there, though it's still worth naming cases clearly, since every result and plot title is derived from the folder name.