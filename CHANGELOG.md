# Changelog

All notable changes to the Flood Frequency Analysis tool are documented here.

## [0.4.4] - 2026-07-27

### Fixed — Dashboard rendering: legend/label collisions, plus a requested addition
- `regional_growth_curve_plot()`: the legend (which defaulted to matplotlib's automatic "best" placement) was landing on top of the y-axis label and tick numbers in the upper-left corner. Moved to a fixed `loc="lower right"`, which is reliably empty (at high T the curve rises well above the data cloud, leaving that corner clear) -- verified programmatically (no bounding-box overlap between the legend and the y-axis label/ticks) rather than by eye. Also shortened/shrank the y-axis label slightly for tighter dashboard embedding.
- `regional_moment_ratio_diagram()`: station name labels were overlapping when stations cluster tightly in (t3, t4) space (the common case for a well-pooling, homogeneous region -- ironically the case you most want to read clearly). Added `_label_points_staggered()`, a small no-dependency helper that orders points by y-value and gives each label a monotonically increasing vertical offset, alternating left/right of the point, connected back with a thin leader line -- verified with a bounding-box overlap check on the text glyphs themselves (isolated from the leader-line arrows, which inflate the naive bbox check).
- `regional_dashboard()`'s text summary panel now includes the T=10,000-year growth factor alongside the existing 10/100/1000-year rows.
- 3 new regression tests targeting these three fixes directly (bbox-overlap checks, not just "does it render") — 172 tests total.

## [0.4.3] - 2026-07-27

### Added — Clearer per-station-vs-pooled comparison plot
- The v0.4.2 pooled-rank overlay on `regional_growth_curve.png` was correct but too subtle: pale grey `x` markers at alpha=0.5, mostly hidden behind the colored per-station points in the dense part of the plot, only visibly distinct in the high-T tail where pooling's inflated apparent sample size reaches further than any single station's own record could.
- `floodfreq/plots.py`: `regional_pooled_vs_stations_plot()` / `save_regional_pooled_vs_stations_plot()` — a new two-panel side-by-side comparison at equal visual weight: left panel is the per-station view (what the growth curve is actually fit to); right panel is every station's data pooled into one array and re-ranked as a single combined record, uniformly colored (station identity is gone once concatenated, which is the point). Both panels share the same y-axis scale and the same fitted regional growth curve, making the difference in apparent reach along the T axis immediately visible rather than requiring the reader to spot faint grey markers.
- `run_regional_analysis.py` now writes `regional_pooled_vs_stations.png` to `Plot/Regional/<RegionName>/` alongside the existing plot set. The original subtle overlay on `regional_growth_curve.png` (`show_pooled_rank`) is kept as-is for anyone who wants it on the single combined plot.
- README's regional-analysis section updated with the new plot and a pointer to it from the existing pooling-methodology note.
- 3 new tests — 169 tests total.

## [0.4.2] - 2026-07-27

### Added — Per-station data quality, Monte-Carlo confidence intervals, and a pooled-rank diagnostic overlay
- `floodfreq/regional.py`:
  - `station_data_quality()` — per-station Mann-Kendall stationarity + Grubbs' outlier + basic input-validation checks, run BEFORE pooling (reuses `floodfreq.data_quality`, the same checks the single-station tool runs). `RegionalAnalysisResult.data_quality_df` now carries this for every station; on the bundled Piedmont demo data it correctly flagged a genuine significant trend in one station (STN_Birchford)
  - `_simulate_growth_curve_quantiles()` / `station_quantile_ci()` / `regional_quantile_ci()` — Monte-Carlo confidence intervals for regional design floods (the regional counterpart of the single-station tool's `bootstrap_ci()`), combining two independent simulated sources of uncertainty: (1) regional growth-curve *shape* uncertainty, via the same kappa-distribution-simulated synthetic homogeneous regions already used for the H1/Z-statistic (Hosking & Wallis, 1997, sec. 5.3's accuracy-assessment approach), refitting the growth curve to each replicate; and (2) each station's own *index-flood* (sample mean) sampling uncertainty, via a nonparametric bootstrap of that station's own record. `regional_quantile_ci()` shares one growth-curve simulation across all stations for efficiency.
- `run_regional_analysis.py`:
  - New `data_quality.csv` and `growth_curve_quantiles_ci.csv` (long format: station, T, design flood, CI bounds) outputs to `Output/Regional/<RegionName>/`
  - `summary.txt` gains a new "STEP 1: PER-STATION DATA QUALITY" section (renumbering the existing steps to 2-5) and, unless `--no-ci`, a "STEP 6: CONFIDENCE INTERVALS" section
  - New CLI flags: `--confidence-level` (default 95) and `--no-ci` (skip the Monte-Carlo CI step, which roughly doubles run time)
- `floodfreq/plots.py`: `regional_growth_curve_plot()` / `save_regional_growth_curve_plot()` gain a `show_pooled_rank` overlay (on by default) -- a light-grey diagnostic scatter of every station's dimensionless growth factor pooled and re-ranked together as one combined record. This is explicitly NOT what the growth curve is fit to (see the new README note on how pooling actually works -- a record-length-weighted average of each station's own L-moment ratios, not concatenation) -- it's a visual sanity check for whether the per-station view and a naive fully-pooled view tell the same story.
- 16 new tests (data quality, Monte-Carlo CI machinery, pooled-rank overlay toggle, expanded end-to-end CLI coverage) — 166 tests total

## [0.4.1] - 2026-07-27

### Added — Regional plot suite (parity with the single-station tool's plots)
- `floodfreq/plots.py`:
  - `regional_growth_curve_plot()` / `save_regional_growth_curve_plot()` — the regional analogue of `probability_plot`: the fitted dimensionless growth curve plotted against every pooled station's own data (rescaled by that station's index flood, a "Dalrymple plot")
  - `station_design_flood_plot()` / `save_station_design_flood_plot()` — the regional analogue of `probability_plot` for ONE station: its own observed data against the design-flood curve implied by the regional growth curve scaled by its index flood
  - `regional_station_series_plot()` (via `save_regional_station_series_plot()`) — small multiples of each pooled station's raw annual-maximum series, a visual screen before trusting the pooling
  - `regional_discordancy_plot()` / `save_regional_discordancy_plot()` — bar chart of each station's D_i against the group's critical value
  - `regional_dashboard()` / `save_regional_dashboard()` — one-page overview (growth curve + L-moment ratio diagram + discordancy chart + text summary panel), the regional analogue of the single-station `dashboard()`
- `run_regional_analysis.py` now writes the full plot set above to `Plot/Regional/<RegionName>/` (previously only `regional_moment_ratio_diagram.png`), including one `station_<StationName>_design_flood.png` per pooled station
- `RegionalAnalysisResult` now carries the raw per-station `station_data` (needed by the new plots) alongside the L-moment summaries
- 8 new tests covering the new plot functions, plus expanded end-to-end coverage in `test_run_regional_analysis.py` — 157 tests total

## [0.4.0] - 2026-07-26

### Added — Regional (pooled) flood frequency analysis
- New entry point `Module/run_regional_analysis.py` and module `floodfreq/regional.py` implementing the classical index-flood method (Dalrymple, 1960) with the Hosking & Wallis (1997) L-moment-based regional statistics (notation follows Rao & Hamed, Ch. 9)
- `floodfreq/regional.py`:
  - `station_lmoments()` — unbiased sample L-moments/ratios per station
  - `discordancy()` — D_i measure flagging stations statistically inconsistent with the group, against Hosking & Wallis's tabulated critical values
  - `heterogeneity()` — H1 measure via Monte-Carlo simulation of synthetic homogeneous regions (kappa distribution, GEV fallback), classifying a region as acceptably homogeneous / possibly heterogeneous / definitely heterogeneous
  - `zstatistics()` / `recommend_family()` — Z-statistic regional distribution selection across 5 candidate families (GLO, GEV, GNO, PE3, GPA), comparing regional-average L-kurtosis against each family's theoretical value (computed exactly via numerical integration of the quantile function, not polynomial approximations)
  - `fit_growth_curve()` / `RegionalGrowthCurve` — dimensionless regional growth curve fit to the pooled, record-length-weighted L-moments; `station_quantile(T, index_flood)` gives each station's design flood
  - `run_regional_analysis()` — ties all of the above into one six-step pipeline call
- New project layout: `Data/Regional/<RegionName>/*.csv` (one file per station), `Output/Regional/<RegionName>/`, `Plot/Regional/<RegionName>/`; `floodfreq.io_utils.resolve_region()` / `load_region_stations()` mirror the existing single-station `resolve_case()` / `read_series()` helpers
- `floodfreq/plots.py`: `regional_moment_ratio_diagram()` / `save_regional_moment_ratio_diagram()` — multi-station L-moment ratio diagram (extends the single-station version) plotting every station's (t3, t4), the record-length-weighted regional average, discordant stations outlined in red, against GLO/GPA (closed-form) and GEV/GNO/PE3 (numeric) reference curves
- `run_regional_analysis.py` CLI: writes `station_lmoments.csv`, `discordancy.csv`, `heterogeneity.csv`, `zstatistics.csv`, `quantile_table.csv`, and a plain-language `summary.txt` (same RUN PROVENANCE header convention as the single-station tool) to `Output/Regional/<RegionName>/`, plus `regional_moment_ratio_diagram.png` to `Plot/Regional/<RegionName>/`. Supports `--candidates`, `--family` (force a specific distribution), `--n-sim`, `--seed`, `--return-periods`, `--no-plots`, and an interactive region picker when run with no arguments
- Added `Data/Regional/Piedmont/`: a 6-station synthetic demo region (shared GEV shape, varying scale/record length) for trying the tool out of the box
- 11 new tests (regional-path coverage added to `test_io_utils.py` and `test_plots.py`, plus new `test_run_regional_analysis.py`), on top of the `test_regional.py` module suite — 150 tests total

## [0.3.0] - 2026-07-26

### Added — Generic (non-streamflow) variable support
- `--variable-name`, `--units`, `--short-name` — the tool can now be used for any annual-maximum series, not just discharge (e.g. annual maximum rainfall depth). Changes every plot axis, chart title, `summary.txt` header, and Excel/PDF report title to match. Fully backward compatible: leaving these unset reproduces the exact original "Flood magnitude" / "Flood Frequency Analysis" wording
- Recognized as `.toml` config keys too (`variable_name`, `units`, `short_name`)
- Documented caveat: `--regional-skew` (Bulletin 17B) uses MSE coefficients empirically calibrated from streamflow data specifically — not adjusted by these flags, and should be avoided for non-streamflow variables

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
