# Flood Frequency Analysis

A Python tool for flood frequency analysis on annual-maximum inflow series, using [uv](https://docs.astral.sh/uv/) for dependency management. Fits 9 candidate distributions (Normal, LogNormal 2p/3p, Gumbel, GEV, Exponential, Gamma, Pearson III, Log-Pearson III) via method of moments / maximum likelihood / probability-weighted moments, ranks them by AIC/BIC/KS, and reports a best-fit recommendation with bootstrap confidence intervals.

## What's included
```
Flood_Frequency_Analysis/
├── Data/            # input data, one CSV per case: Data/<CaseName>.csv (not tracked in git)
├── Module/          # the floodfreq package + run_analysis.py
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
1. Fit all 9 distributions using their recommended method (MLE for Normal/LogNormal, PWM for the extreme-value/skewed families, MOM for Gamma)
2. Rank them by AIC/BIC and a Kolmogorov-Smirnov goodness-of-fit check
3. Prompt you for a confidence level (e.g. `95`) unless `--confidence-level` is passed
4. Write CSVs to `Output/<CaseName>/` (descriptive stats, goodness-of-fit table, quantile table up to T=10,000 years, bootstrap CI, and a `summary.txt` with the recommendation and an extrapolation-risk warning for large T)
5. Write PNG plots to `Plot/<CaseName>/`:
   - `input_data_histogram.png` — histogram of the raw series with the recommended distribution's fitted PDF overlaid
   - `flood_frequency_curve.png` — all fitted distributions vs. observed data, extrapolated out to T=10,000 years
   - `moment_ratio_diagram.png` — L-moment ratio diagram for distribution-family diagnosis
   - `bestfit_<distribution>_ci.png` — the recommended distribution's quantile curve with a bootstrap CI band

Useful flags:
```bash
uv run python Module/run_analysis.py P1009 --confidence-level 90
uv run python Module/run_analysis.py P1009 --plotting-position gringorten
uv run python Module/run_analysis.py P1009 --n-boot 2000 --xlsx-report
uv run python Module/run_analysis.py P1009 --no-plots
```

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