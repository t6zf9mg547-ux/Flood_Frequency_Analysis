# Changelog

All notable changes to the Flood Frequency Analysis tool are documented here.

## [0.2.0] - 2026-07-26

### Added — Run provenance & per-case configuration
- Every `summary.txt` (and PDF report) now opens with a "RUN PROVENANCE" header: timestamp, tool version, the exact command used, Python version, and which config file (if any) supplied settings — so a result found later, or handed to someone else, is self-documenting
- `floodfreq.__version__`, stamped into that header
- Optional per-case settings file, `Data/<CaseName>.toml` — save `confidence_level`, `regional_skew`, `n_boot`, `pdf_report`, etc. once instead of retyping CLI flags every run. Precedence: explicit CLI flag > config file > hardcoded default. A case with no `.toml` behaves exactly as before

## [0.1.0] - 2026-07-26

### Added — Initial port
- Python re-implementation of the original Excel/VBA flood frequency analysis workbook (`Analyse_Frequentielle_V03.xls`)
- 9 candidate distributions — Normal, LogNormal (2p/3p), Gumbel, GEV, Exponential, Gamma, Pearson III, Log-Pearson III — fit via method of moments (MOM), maximum likelihood (MLE), or probability-weighted moments (PWM)
- 6 empirical plotting-position formulas (Weibull, Hazen, Cunnane, Gringorten, Hosking, Blom), matching the original workbook's selector
- Standard project layout (`Data/`, `Module/`, `Output/`, `Plot/`) with `run_analysis.py` as the CLI entry point; `uv`-managed dependencies
- Interactive input-file picker (with row count / last-modified hints) and interactive confidence-level prompt when not passed via CLI flags
- Per-distribution recommended (method, plotting-position) pairs, matching the workbook's own per-tab defaults, with global overrides available
- Design flood quantile table extended to T = 10,000 years, with per-row extrapolation-risk flagging and an explicit warning when T far exceeds the record length

### Added — Step 1: Data quality
- Mann-Kendall trend/stationarity test, with Sen's slope for the trend line
- Grubbs' outlier test (ASTM E178), applied in log-space
- Basic input validation warnings (missing years, short record, etc.)
- `data_quality_timeseries.png` — annual maxima with trend line and flagged outliers highlighted
- `data_quality.csv`, with a plain-language `meaning` column

### Added — Step 2: Fitting & extrapolation
- Anderson-Darling goodness-of-fit statistic (tail-weighted, complementing KS)
- Akaike-weighted multi-model averaging: `akaike_weights.csv`, `model_averaged_quantiles.csv`, and a model-averaged design flood table in `summary.txt`
- Bulletin 17B regional skew weighting for Pearson III / Log-Pearson III (`--regional-skew`, `--regional-skew-mse`), as additional candidates alongside the standard PWM fits
- Plain-language Good / Caution / Warning assessments on goodness-of-fit, model confidence, and model agreement, calibrated against real output rather than arbitrary thresholds

### Added — Step 3: Software robustness
- `Module/tests/`: pytest regression suite (88 tests), including a reference dataset embedded directly in `conftest.py` (independent of the gitignored `Data/` folder)
- Hard input validation in `FloodFrequencyAnalysis` (NaN, Inf, non-positive values, too-short records, mismatched years) with clear, actionable error messages
- Hardened `read_series()`: clear errors for missing files, empty files, missing columns, all-NaN columns

### Added — Step 4: Usability
- `dashboard.png` — one-page overview combining all 5 individual plots plus a text summary panel
- `--pdf-report` — bundles the text summary, dashboard, and every individual full-size plot into one PDF

### Fixed
- Mann-Kendall sign convention was inverted (an early version reported "decreasing" for a clearly increasing series and vice versa) — caught by testing against synthetic data with a known trend direction before shipping
- Log-likelihood was missing a Jacobian correction for log-transformed distributions (LogNormal, Log-Pearson III), making AIC/BIC not actually comparable across distributions
- `bootstrap_ci()` was not reusing the same plotting position as the original fit, so the CI band and point estimate could be computed slightly inconsistently
- An edit accidentally deleted `save_probability_plot`'s function definition while inserting the dashboard code, leaving its body as unreachable dead code — caught immediately by the test suite

### Changed
- Confidence interval level is now user-specified (`--confidence-level`, or an interactive prompt) instead of a hardcoded 95%
- Confirmed (and documented) that PWM parameter estimation always uses the standard unbiased L-moment estimator — the "plotting position" choice affects only where empirical points are drawn on plots and in the descriptive-statistics summary, not the fitted parameters themselves (an initial attempt to tie PWM fitting to the plotting-position choice was tested against the reference workbook, found to be a worse match, and reverted)