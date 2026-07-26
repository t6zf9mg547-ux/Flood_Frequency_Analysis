# Flood Frequency Analysis

A Python tool for flood frequency analysis on annual-maximum inflow series, using [uv](https://docs.astral.sh/uv/) for dependency management. Fits 9 candidate distributions (Normal, LogNormal 2p/3p, Gumbel, GEV, Exponential, Gamma, Pearson III, Log-Pearson III) via method of moments / maximum likelihood / probability-weighted moments, ranks them by AIC/BIC/KS/Anderson-Darling, and reports a best-fit recommendation (plus an Akaike-weighted multi-model average) with bootstrap confidence intervals.

## What's included
```
Flood_Frequency_Analysis/
├── Data/            # input data, one CSV per case: Data/<CaseName>.csv (not tracked in git)
├── Module/
│   ├── run_analysis.py  # CLI entry point
│   ├── floodfreq/       # the package
│   └── tests/           # pytest regression suite (81+ tests)
├── Output/          # generated CSVs, one subfolder per case: Output/<CaseName>/ (not tracked in git)
├── Plot/            # generated PNG plots, one subfolder per case: Plot/<CaseName>/ (not tracked in git)
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
- `Data/`, `Output/`, and `Plot/` are excluded from git by default, since they typically hold large or regenerated files. Adjust `.gitignore` per-project if you want any of them tracked.
- Reading legacy `.xls` (not `.xlsx`) input directly would additionally need `xlrd` (`uv add xlrd`) — not included by default since the expected input is `Data/<CaseName>.csv`.
- The "plotting position" formula (Weibull, Hazen, Cunnane, Gringorten, Hosking, Blom) affects where empirical points are drawn on probability plots and the descriptive-statistics summary — it does **not** change PWM-fitted parameters, which always use the standard unbiased L-moment estimator. See `summary.txt` for the full explanation, generated with each run.
- Quantiles are reported out to T=10,000 years by default, but treat anything much beyond ~2-3x your record length as increasingly uncertain (`summary.txt` flags this per-row and adds an explicit warning when T exceeds 10x the record length) — that far out, the design flood depends more on which distribution's tail you trust than on the data itself.
- Keep `Data/` free of stray/example CSVs you don't intend to analyze (e.g. an old `Template.csv`). The interactive picker (when running with no case name) shows row count and last-modified time for each file to help tell them apart, but it's easy to select the wrong one by accident, and every result and plot title is derived from that file's name.
