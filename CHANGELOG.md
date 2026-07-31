# Changelog

All notable changes to the Flood Frequency Analysis tool are documented here.

## [0.9.0] - 2026-07-30

### Added — Monte Carlo CIFAM for GEV and Log-Pearson III (plus Normal, Gamma, Exponential)
- The distribution-agnostic Monte Carlo path (`mc_climate_adjusted_quantiles`)
  now supports **GEV** and **Log-Pearson III**, and the previously-wired
  **Normal**, **Gamma (2p)** and **Exponential** are now covered explicitly. A
  new `MC_SUPPORTED_DISTRIBUTIONS` constant lists them. The three closed-form
  families (Gumbel, Log-Normal 2p, Pearson III) are unchanged.
- For the **three-parameter skewed families (GEV, Log-Pearson III)** the CIFAM
  assumption "climate acts on the first two moments" is realized by holding the
  **shape/skew parameter fixed** at its baseline-fitted value and moving only
  location/scale to hit each shifted (mean, standard deviation):
  - GEV: the shape `c` is fitted once (PWM) from the baseline; each realization
    solves `scale = sigma / s(c)`, `loc = mu - scale*m(c)` where `m(c), s(c)`
    are the standardized mean/sd of `genextreme(c)` — an exact moment match
    (verified to 1e-9).
  - Log-Pearson III: the skew of `log10(x)` is fitted once and held fixed; the
    shifted physical (mu, sigma) map to shifted log-space (mu_y, sigma_y), then
    to `pearson3(skew_fixed, mu_y, sigma_y)`.
- Fixed a latent bug: the earlier Monte Carlo moment-mapping produced only two
  parameters for `logpearson3` (a log-normal-style `(mu_y, sigma_y)`), which is
  invalid for the three-parameter `pearson3` scipy object — so `logpearson3`
  Monte Carlo could not have worked before this release. It now maps to the full
  three-parameter form.
- `run_climate_adjustment.py`: `--distribution` help now lists the extra
  Monte-Carlo-only families, and `--method monte-carlo` validates the
  distribution against `MC_SUPPORTED_DISTRIBUTIONS` with a clear error.
- Monte Carlo default raised from 20,000 to **30,000** simulations (both the
  `mc_climate_adjusted_quantiles` `n_sim` argument and the CLI `--n-sim`
  default). The single `n_sim` drives both the sampling-only baseline ensemble
  and the combined climate ensemble, so both confidence intervals get the extra
  precision. Rationale: with quantiles extrapolated out to T = 10,000+ years, a
  measurement of the seed-to-seed scatter of the T = 10,000 climate CI bound
  gave ~0.58% at 2,000 sims, ~0.18% at 20,000, and ~0.11% at 50,000; 30,000 is
  a conservative middle (~0.15%) at only ~1.5x the runtime of 20,000. Override
  per run with `--n-sim` when a final design figure warrants it.
- Regenerated the shipped climate-adjustment template so its `distribution`
  field documents the full set: closed-form gumbel/lognormal2/pearson3, plus
  gev/logpearson3/normal/gamma2/exponential via `--method monte-carlo`.

### Tests
- New tests in `test_climate_adjustment.py`: GEV and Log-Pearson III run through
  Monte Carlo with finite, positive, monotone quantiles; the combined CI is
  wider than the sampling-only baseline CI; zero climate change collapses the
  climate estimate onto the baseline; the GEV loc/scale rescale hits the target
  moments exactly with the shape held fixed; and the closed-form entry point
  still refuses GEV.

## [0.8.0] - 2026-07-30

### Changed (BREAKING) — case-first `Data/` layout, multiple climate scenarios per case
This release makes the `Data/` folder case-first, matching how `Output/` and
`Plot/` are already organised, and lets each case hold several climate-adjustment
scenarios. It is a **breaking change to input file locations** — existing data
must be moved (see "Migration" below). Regional analysis is unchanged.

- **Single-station baseline** moved from `Data/<CaseName>.csv` to
  `Data/<CaseName>/<CaseName>.csv` (baseline named after its folder, by
  convention). `resolve_case()` and `run_analysis.py`'s interactive picker now
  work in terms of case folders under `Data/` (skipping the reserved `Templates/`
  and `Regional/` folders).
- **Per-case settings** moved from `Data/<CaseName>.toml` to
  `Data/<CaseName>/<CaseName>.toml`.
- **Climate-adjustment inputs** moved from `Data/Climate_Adjustment/<CaseName>.csv`
  (a single file per case, v0.6.0–0.7.0) to
  `Data/<CaseName>/climate_adjustment/<scenario>.csv` — **one file per scenario**
  (e.g. `rcp45.csv`, `rcp85.csv`). `run_climate_adjustment.py` gains a second
  positional argument, the scenario: `run_climate_adjustment.py <CaseName> <scenario>`
  (either is prompted for if omitted). `resolve_climate_case()` now takes a
  `scenario` argument and `ClimatePaths` gains a `scenario` field.
- **Climate outputs** are now per-scenario:
  `Output/<CaseName>/Climate_Adjustment/<scenario>/` and
  `Plot/<CaseName>/Climate_Adjustment/<scenario>/`, so scenarios never overwrite
  each other.
- **Shipped examples/templates consolidated under `Data/Templates/`.** All
  example and template files now live in one tracked subtree, which lets
  `.gitignore` collapse its Data rules to essentially two lines: ignore all of
  `Data/`, except `Data/Templates/` (and `Data/.gitkeep`). The earlier scattered
  per-file exceptions (`Data/Template.csv`, `Data/Climate_Adjustment/Template.csv`,
  `Data/Regional/Template/*.csv`, the descriptor catalogs and station-data
  folders) are all replaced by that single subtree exception. Verified with an
  actual `git add -A` against a reconstructed tree: only `Data/Templates/**`
  tracks; every real case, regional data, and all outputs are ignored.

### Migration (do this once, manually)
For each existing case `<CaseName>`:
1. `mkdir Data/<CaseName>/`
2. move `Data/<CaseName>.csv` → `Data/<CaseName>/<CaseName>.csv`
3. if present, move `Data/<CaseName>.toml` → `Data/<CaseName>/<CaseName>.toml`
4. if you had a climate file `Data/Climate_Adjustment/<CaseName>.csv`, move it to
   `Data/<CaseName>/climate_adjustment/<scenario>.csv` (pick a scenario name)
Then move every shipped template/example under `Data/Templates/` and drop the old
`Data/Climate_Adjustment/` folder. Nothing in `Output/` or `Plot/` needs moving
(those are regenerated).

### Tests
- Updated `test_io_utils.py` (case-first `resolve_case`/`load_case_config` paths,
  new `resolve_climate_case(case, scenario, ...)` assertions, plus a
  scenarios-don't-collide test), `test_run_climate_adjustment.py` (case+scenario
  layout and CLI args), and `test_climate_adjustment.py` (end-to-end).

## [0.7.0] - 2026-07-30

### Added — CIFAM command-line entry point (`run_climate_adjustment.py`) + Figure-4-style plot
- New `Module/run_climate_adjustment.py`, so the climate adjustment is run the same way as the other tools rather than only as a library call: `uv run python Module/run_climate_adjustment.py <CaseName>`. It reads the baseline series from `Data/<CaseName>.csv` and the four climate numbers from `Data/Climate_Adjustment/<CaseName>.csv`, and writes three outputs:
  - `Output/<CaseName>/Climate_Adjustment/climate_adjustment_table.csv` — per return period, baseline point + CI (sampling only) and climate point + combined CI.
  - `Output/<CaseName>/Climate_Adjustment/summary.txt` — a provenance-stamped summary (matching `run_analysis.py`'s header style) with the inputs and the full table, plus a PMF comparison line when a PMF is given.
  - `Plot/<CaseName>/Climate_Adjustment/climate_adjustment.png` — a Figure-4-style plot.
- **Output layout changed to case-first** from the `Output/Climate_Adjustment/<CaseName>/` scheme introduced in v0.6.0 (which never shipped a CLI to populate it, so nothing downstream depended on it). Climate outputs now nest under the same `Output/<CaseName>/` and `Plot/<CaseName>/` that `run_analysis.py` already creates, so a case's baseline results and its climate adjustment live together. The inputs stay analysis-type-first at `Data/Climate_Adjustment/<CaseName>.csv`, so the shipped `Template.csv` and its single `.gitignore` exception are unchanged. (Regional analysis remains analysis-type-first on purpose — a region isn't a single case, whereas a climate adjustment is an extension of one.)
- Mirrors `run_analysis.py`'s conventions: an interactive case picker when no name is given, `--value-col`/`--year-col`/`--variable-name`/`--units`/`--short-name`, and a provenance header. Adds CIFAM-specific flags that override the inputs file for one-off scenarios: `--delta1/--delta2/--tau1/--tau2`, `--distribution`, `--confidence-level`, `--pmf`, `--return-periods`, `--method {closed-form,monte-carlo}`, `--n-sim`, `--random-state`, `--no-plot`, and `--return-period-axis`.
- New `plots.climate_adjustment_plot()` / `save_climate_adjustment_plot()`: baseline (black central line + grey sampling-only band) vs climate-adjusted (red central line + wider red combined band, dashed edges) on the Gumbel reduced-variate axis with return-period markers, and an optional PMF reference line — the Fig. 4 layout from Grijsen & Lino (2026).

### Fixed — value/year column unpacking in the CIFAM data path
- `io_utils.read_series()` returns `(values, years)` in that order; the initial CIFAM code (the README example, the end-to-end test, and the new run script) unpacked it as `years, values`, which silently fed the year labels in as the "flood series" and produced year-valued quantiles (~2000 m³/s) with an absurdly small standard deviation. Corrected everywhere. The bug had hidden in the end-to-end test because it compared two calls that both used the same mis-unpacked variable; the test now asserts the baseline mean lands in the flood-series range (not the ~1970 year range) and that the flood series, not the years, drives the quantiles, so a future swap fails loudly.

### Tests
- New `Module/tests/test_run_climate_adjustment.py` (5 tests): the `--flag`-overrides-file precedence helper (including the 0.0-is-not-falsy case), and two in-process end-to-end runs that exercise the whole CLI pipeline (CSV + summary.txt + plot written, sane values, combined CI wider than baseline CI, `--delta1` override and `--no-plot` honoured).

## [0.6.0] - 2026-07-29

### Added — CIFAM: climate-informed flood frequency adjustment (single-station), items 1 & 2
- New module `Module/floodfreq/climate_adjustment.py` implementing the rapid Climate-Informed Flood Assessment Methodology (CIFAM) of Grijsen & Lino (ICOLD 2026). It takes a baseline annual-maximum series plus four user-supplied, flow-space climate-change numbers — δ1/δ2 (projected % increase in the mean and its ensemble spread) and τ1/τ2 (same for the inter-annual standard deviation) — and returns, per return period, the baseline quantile + confidence interval (sampling uncertainty only, unchanged from classical FFA) alongside a climate-shifted central estimate + a *combined* confidence interval that adds climate uncertainty to sampling uncertainty. The four climate numbers are inputs, not something the tool derives (same boundary principle as pooling-group formation not doing your GIS extraction).
- **Two independent paths, cross-validated against each other:**
  - `climate_adjusted_quantiles(...)` — closed form for Gumbel, Log-Normal (2p), and Pearson III (the three distributions the paper gives formulas for).
  - `mc_climate_adjusted_quantiles(..., n_sim=20000)` — a distribution-agnostic Monte Carlo generalization (draws (μ,σ) realizations, maps to parameters via method of moments, takes quantile percentiles). This is the path that will later extend to the remaining candidate distributions, for which no closed form exists (items 3–4, not in this release).
- **The math was re-derived from first principles, not transcribed** (the source PDF had OCR damage in the Pearson III row of Table 1 and an ambiguous constant in eq. 8). Two transcription questions are now settled and documented in the module: (1) eq. 8's final `(tp − c1)²` uses Euler's constant c = 0.57722 — `c1` was an OCR artifact, not a separate constant; (2) the Pearson III `d(ln γ)`, `d(ln β)` and `∂xp/∂γ` formulas are correct once the dropped subscripts are restored (confirmed by finite-difference to ~1e-8, and the full `var(xp)` rebuilt from the verified partials).
- **A real correction to the paper's Log-Normal formula.** Cross-checking each closed-form climate variance against an exact Monte Carlo showed the paper's Log-Normal variance term overestimates the true variance by ~45% at realistic coefficients of variation — it relies on the small-σy simplification `σx ≈ μx·σy`. (Gumbel's and Pearson III's formulas, by contrast, match the Monte Carlo to <0.6%.) The Log-Normal closed form therefore uses an exact delta-method variance derived here, which agrees with the Monte Carlo; the paper-verbatim term is preserved in a comment for anyone wanting to reproduce the paper's exact figures.
- A subtlety worth recording: the Log-Normal and Pearson III quantiles use the standard-normal reduced variate Φ⁻¹(p), *not* the Gumbel reduced variate −ln(−ln p). The paper writes "tp" for all three, but the symbol means each family's own variate, and they differ substantially (at T = 10⁴: 9.21 vs 3.72). `reduced_variate(T, kind=...)` makes this explicit.
- **Validation against the paper's Lom Pangar (Cameroon) worked example** (Fig. 4). Since the paper gives no numeric baseline Gumbel parameters, the baseline was digitized from Fig. 4 (axes calibrated off the blue 10⁴/10⁵-yr reference lines and the y-gridlines; the red "with-CC" lines extracted and inverted through the CIFAM shift): u0 ≈ 1190, β0 ≈ 165 → μ0 ≈ 1285, σ0 ≈ 211, N ≈ 80. Driving the independent closed form with that baseline reproduces the paper's headline results — climate 97.5% upper ≈ 46–52% above the baseline point estimate ("broadly ~50%"), ≈ 32–42% above the baseline 97.5% CI (paper: 30–35%), climate 97.5% at the 100,000-yr level ≈ 4,499 m³/s against the stated PMF of 4,140, and a design/safety range of ≈ 3,440–4,500 against the stated 3,500–4,500.
- **Project-layout integration.** Added a `Climate_Adjustment/` namespace parallel to `Regional/`: the four climate numbers live in `Data/Climate_Adjustment/<CaseName>.csv` (a small, self-documenting long-format "parameter,value,units,description" file), while the baseline series is reused in place from `Data/<CaseName>.csv` — no data duplication. `Output/Climate_Adjustment/<CaseName>/` and `Plot/Climate_Adjustment/<CaseName>/` follow the existing per-case output convention.
  - New `Data/Climate_Adjustment/Template.csv` to duplicate per case, pre-filled with the paper's Lom Pangar numbers as a worked example.
  - New `io_utils.resolve_climate_case()` (mirrors `resolve_region()`: derives all four paths, creates the output/plot dirs) and `io_utils.load_climate_inputs()` (parses the template; requires delta1/delta2/tau1/tau2 as fractions and rejects magnitudes > 1.5 to catch the "entered 30 instead of 0.30" mistake; fills sensible defaults for confidence_level/distribution/return_periods/pmf).
  - `.gitignore`: added negation exceptions so `Data/Climate_Adjustment/Template.csv` is tracked while everything else under the namespace — a user's own `Data/Climate_Adjustment/<CaseName>.csv`, and all of `Output/`/`Plot/` — stays untracked, verified with an actual `git add -A` dry run against a reconstructed project tree (not just pattern reasoning). Same track-the-template/ignore-the-rest approach as the existing Regional and Templates exceptions.
- 25 new tests: 16 in `Module/tests/test_climate_adjustment.py` (hand-computable toy cases pinning each distribution's closed-form algebra, the Lom Pangar reproduction — asserted with loose tolerances, documented as figure-digitized/approximate — Monte-Carlo-vs-closed-form cross-validation for all three distributions, and an end-to-end template→result pipeline test), plus 9 in `Module/tests/test_io_utils.py` covering `resolve_climate_case`/`load_climate_inputs`.
- **Scope**: single-station only; closed forms for Gumbel/Log-Normal/Pearson III plus the Monte Carlo path, with the folder contract and inputs template in place. Extending the Monte Carlo path to the other candidate distributions, a `run_analysis.py`-style CLI that writes the CSV/PNG outputs, and a Fig-4-style plot are deliberately left for a later release (items 3–4).

## [0.5.6] - 2026-07-28

### Fixed — README gave the impression this tool only does single-station analysis
- Replaced the opening paragraph (which described only the 9-distribution single-station fitting) with a three-objective overview, giving equal footing to: (1) single-station analysis, (2) regional (pooled) analysis, and (3) automatic pooling-group formation -- each with a working jump-link to its full section.
- No code changes this round -- README.md only.

## [0.5.5] - 2026-07-28

### Fixed — Documentation gap: when do you need form_pooling_group.py, and what exactly do you edit?
- The README's pooling section and the CLI's own docstring were technically complete (every flag, the ranking math, working examples) but never actually answered the two questions that came up in practice: (1) do you need this tool at all, if you already know which stations to pool by hand, and (2) what precisely do you duplicate and edit to use it with your own data?
- Added a "Do you actually need this?" framing paragraph up front in both the README section and the CLI's `--help` text: skip this entirely if you already know your stations, just build `Data/Regional/<RegionName>/` directly (e.g. by copying `Data/Regional/Template/`).
- Added an explicit 4-step recipe to the top of the README section: (1) duplicate the descriptor catalog template, (2) duplicate the matching station-data folder, (3) run `form_pooling_group.py --apply` against your own files, (4) run `run_regional_analysis.py` on the result -- including the easy-to-miss requirement that a station CSV's filename must exactly match its `station` value in the catalog.
- No code changes this round -- README.md and `Module/form_pooling_group.py`'s docstring only.

## [0.5.4] - 2026-07-28

### Added — Bundled annual-maximum series matching the descriptor-catalog templates
- You noticed `Data/Templates/` only had the candidate *descriptor* catalogs (v0.5.1) with no matching input data, so `--apply` had nothing to actually copy without supplying your own `--station-data-dir`. Added `Data/Templates/discharge_station_data/` (10 CSVs, `GAUGE_001.csv`..`GAUGE_010.csv`, columns `year,Q`) and `Data/Templates/rainfall_station_data/` (10 CSVs, `RAINGAUGE_001.csv`..`RAINGAUGE_010.csv`, columns `year,Rainfall_mm`).
- Generated coherently with each catalog: every station's record length matches its catalog row's `n_years` exactly, and its synthetic magnitude is scaled (loosely, illustratively -- not a calibrated regionalization formula) with its own `area_km2`/`mean_annual_precip_mm`, so e.g. the largest discharge catchment (GAUGE_005, 455.8 km²) produces visibly larger flows than the smallest (GAUGE_008, 72.6 km²) -- confirmed directly (mean 76.4 vs. 41.7), not just asserted. Both sets use a single common GEV shape parameter within their group, so the bundled group reads as homogeneous by construction rather than needing a discordant-station workaround.
- Verified the full pipeline end to end using ONLY bundled content, not just unit tests: ran `form_pooling_group.py --apply` against each bundled catalog + station-data folder, then fed the result straight into `run_regional_analysis.py` -- got a real fitted growth curve out the other end for both discharge and rainfall.
- `.gitignore` updated with exceptions for both new subfolders; re-verified with an actual `git add -A` dry run against the real project tree.
- README's pooling section now documents the matching station data and adds a fully-bundled `--apply` → `run_regional_analysis.py` example using no external files at all.
- 2 new parametrized end-to-end tests (`test_bundled_station_data_matches_catalog_and_is_fully_runnable`) -- checks record-length coherence against the catalog for every station, runs a real `--apply`, and confirms the resulting region is runnable -- 215 tests total.

## [0.5.3] - 2026-07-27

### Removed
- The bundled Piedmont synthetic demo region (`Data/Regional/Piedmont/`, introduced in [0.4.0]) has been removed at the user's request, in favor of their own `Data/Regional/Template/` folder. `.gitignore`'s tracked-exception list dropped `Data/Regional/Piedmont/*.csv` accordingly, keeping only `Data/Regional/Template/*.csv` (plus the unrelated `Data/Template.csv` and `Data/Templates/candidate_descriptors_*_TEMPLATE.csv`, which are unaffected by this change).
- README's regional-analysis walkthrough no longer references Piedmont or claims a bundled 6-station demo region; CLI examples now use `Template` as the placeholder region name, consistent with the single-station tool's own `Template.csv`/`Output/Template/`/`Plot/Template/` convention. (What actually lives in `Data/Regional/Template/` is up to you -- unlike Piedmont, this isn't a dataset shipped or verified by me.)
- `test_io_utils.py`'s `resolve_region()` path tests used "Piedmont" only as an arbitrary string against an isolated `tmp_path` fixture, not the real bundled folder -- renamed to "Template" for consistency; no test behavior depended on the real Piedmont data, confirmed by the full suite still passing (213/213) after removal.
- Older changelog entries mentioning Piedmont (e.g. verification notes in [0.4.0] through [0.5.2]) are left as-is, since they're an accurate historical record of what was tested at the time.

## [0.5.2] - 2026-07-27

### Fixed — Descriptor templates were visible to the single-station picker
- You asked me to double-check for anything else incomplete after v0.5.1, and this is what turned up: `run_analysis.py`'s interactive case picker lists `Data/*.csv` (non-recursive) as candidate single-station cases. The two candidate-catalog templates added in v0.5.1 lived directly in `Data/`, so they showed up in that picker right alongside real annual-maximum cases. Verified this was a real problem, not just a theoretical one: fed `candidate_descriptors_discharge_TEMPLATE.csv` through `read_series()` directly, and it silently picked `area_km2` as the "flood series" with no year detected and no error -- someone picking the wrong file from the interactive list would have gotten a nonsense analysis with no warning that anything was wrong.
- Moved both templates to `Data/Templates/` (a subfolder, which the picker's non-recursive glob never sees): `Data/Templates/candidate_descriptors_discharge_TEMPLATE.csv` and `Data/Templates/candidate_descriptors_rainfall_TEMPLATE.csv`.
- Updated `.gitignore`'s exceptions, README (folder tree, usage examples, the Notes section), and `test_form_pooling_group.py` for the new path. Added `test_bundled_templates_are_not_visible_to_single_station_picker` as a regression test so this can't quietly regress.
- Re-verified the `.gitignore` fix from v0.5.1 still works correctly at the new path with an actual `git add -A` dry run against the real project tree.
- Also fixed two now-stale README passages caught in the same audit pass: the "Notes" section still described `Data/`/`Output/`/`Plot/` as unconditionally excluded from git (true before v0.5.1, no longer true given the tracked-exception list), and a piece of advice about keeping `Data/` free of "an old Template.csv" that read as self-contradictory now that `Template.csv` is an intentionally-shipped, git-tracked demo file rather than clutter.
- 1 new test -- 213 tests total.

## [0.5.1] - 2026-07-27

### Added — Bundled candidate-catalog templates for `form_pooling_group.py`
- `Data/candidate_descriptors_discharge_TEMPLATE.csv`: 10 synthetic streamflow gauges with `area_km2`, `mean_annual_precip_mm`, `bfihost`, `mean_slope_m_per_km`, `urban_extent_frac`, `n_years`.
- `Data/candidate_descriptors_rainfall_TEMPLATE.csv`: 10 synthetic rain gauges with `lat`, `lon`, `elevation_m`, `mean_annual_precip_mm`, `n_years`.
- These follow up on v0.5.0, which described the candidate-catalog schema in prose (README + docstrings) but never shipped an actual file to open, edit, or point the CLI at -- these close that gap. Column names are illustrative; real catalogs can use whatever names your own GIS extraction produces, since `form_pooling_group.py --descriptors` takes arbitrary column names.
- Verified both work out of the box, not just described: ran `form_pooling_group.py` against each directly (no editing), confirmed sane rankings (nearest-neighbor candidates rank first with small distances in both), and added a parametrized regression test (`test_bundled_template_catalog_works_out_of_the_box`) so this stays true going forward rather than silently rotting.
- README's pooling section now shows ready-to-run example commands against both bundled templates.

### Fixed
- `.gitignore` previously excluded everything under `Data/` with no exceptions beyond `.gitkeep` (`/Data/*`) -- meaning the Piedmont demo region (shipped in v0.4.0) and `Data/Template.csv` were **never actually trackable in git** despite being meant to ship with the repo, and the same problem would have applied to these two new templates. Added proper negation exceptions for `Data/Template.csv`, both new template CSVs, and `Data/Regional/Piedmont/*.csv`, while confirming (via an actual `git add -A` dry run against the real project tree, not just reasoning about the patterns) that everything else under `Data/`, `Output/`, and `Plot/` -- including any *other* regions someone creates locally -- still stays untracked as intended. This was flagged as a caveat back in the `git add -A` conversation but never actually fixed until now.

## [0.5.0] - 2026-07-27

### Added — Automatic pooling-group formation ("region of influence")
- New `floodfreq/pooling.py`: descriptor-source-agnostic implementation of the "region of influence" approach (Burn, 1990) to proposing a pooling group, rather than assembling one by hand.
  - `read_candidate_catalog()` -- load a pool of candidate stations, each with a `station` name and arbitrary numeric descriptor columns (plus an optional `n_years` for the min-total-years stopping rule).
  - `resolve_target()` -- get the target site's descriptor values, either from an existing candidate (leave-one-out style: "if I pooled around this gauge, who looks similar?") or from a supplied dict (the normal ungauged-target case).
  - `similarity_ranking()` -- rank every other candidate by weighted Euclidean distance to the target in z-scored descriptor space (standardized using the candidate pool's own mean/std, so descriptors on very different natural scales -- e.g. catchment area in km² vs. a 0-1 soil index -- contribute comparably).
  - `propose_pooling_group()` -- apply a stopping rule to the ranking: `n_stations` (take the N most similar) or `min_total_years` (accumulate the most similar candidates until a target total station-years is reached, per Hosking & Wallis's ~5× design-return-period rule of thumb), with a warning if the whole candidate pool falls short.
- New `Module/form_pooling_group.py` CLI: reads a candidate catalog, resolves a target (gauged or ungauged), ranks and proposes a group, writes `pooling_ranking_full.csv` / `pooling_group_proposed.csv` to `Output/Regional/<RegionName>/`, and (with `--apply` + `--station-data-dir`) copies the proposed stations' own annual-maximum CSVs into `Data/Regional/<RegionName>/`, ready for `run_regional_analysis.py`. Defaults to a dry run (proposal only, nothing written to `Data/`) unless `--apply` is explicitly given.
- Deliberately does NOT extract descriptors itself or bundle any GIS dependency -- it consumes whatever numeric columns you hand it. Since the tool supports both streamflow and rainfall (and other generic annual-maximum variables), the natural descriptor set differs by type: catchment area / mean annual precipitation / soil-permeability index / slope / urban extent for streamflow gauges (a natural source: HydroSHEDS/BasinATLAS, snapping each gauge to its containing sub-basin), vs. geographic coordinates / elevation / point climatology for rain gauges (which have no upstream catchment for basin-level attributes to mean anything -- BasinATLAS's own precipitation attribute is itself just WorldClim data aggregated over a basin polygon, so a rain gauge is better served sampling the same underlying climate grid directly at its point coordinates).
- Verified end-to-end, not just via unit tests: ran the actual CLI against a synthetic candidate catalog + station-CSV pool, confirmed the near-twin candidate ranks first with a small distance, confirmed `--apply` copies exactly the proposed stations, and confirmed the resulting `Data/Regional/<RegionName>/` folder is directly usable by `run_regional_analysis.py` (ran it, got a real fitted growth curve).
- 32 new tests: 25 in `test_pooling.py` (ranking math, both stopping rules, all error paths, standardization checked against hand-computed z-score distances) + 7 in `test_form_pooling_group.py` (subprocess-based CLI integration tests, including the dry-run/`--apply` distinction and cleanup so no test artifacts leak into the real project tree) -- 210 tests total.
- This is a proposal step, explicitly not a substitute for discordancy/heterogeneity screening -- the CLI's own output says as much, and the README documents it.

## [0.4.6] - 2026-07-27

### Added — Calendar years tracked per regional station
- `io_utils.load_region_stations()` now returns `(station_data, station_years)` instead of a single dict. `station_years[name]` is the station's real year array when its CSV has a detectable year column (auto-detected the same way as single-station cases, or via `--year-col`), or `None` otherwise. **This is a breaking change to `load_region_stations()`'s return signature** -- any code calling it directly needs to unpack the tuple.
- `regional.station_data_quality()`, `regional.run_regional_analysis()`, and `RegionalAnalysisResult` all gained a `station_years` parameter/field. When years are available for a station, the Mann-Kendall/Sen's-slope/validation checks use real calendar time instead of falling back to record order: Sen's slope is now reported in true per-calendar-year units (previously per-index-step, which is only meaningful if every year is actually present with no gaps), and the "year range spans N years but there are only M values" missing-year check (already used by the single-station tool) is now available for regional stations too. A new `years_available` boolean column in `data_quality.csv` flags which stations have real years vs. the record-order fallback.
- `run_regional_analysis.py`: console output now prints each station's year range when available; `summary.txt`'s Step 1 (data quality) section notes how many stations, if any, are missing a year column and therefore fell back to record order.
- `plots.save_regional_station_series_plot()`: the small-multiples station series plot now plots against real calendar year when available, falling back to observation order (1..n) otherwise -- confirmed on the bundled Piedmont demo data (all 6 stations have year columns, spanning e.g. 1933-2022 for STN_Dunmore).
- Verified end-to-end against the actual Piedmont demo data (not just unit tests): per-station year ranges print correctly, Sen's slope values changed from index-step to true per-year units, and the station-series plot's x-axis now shows real years.
- 8 new/updated tests across `test_io_utils.py`, `test_regional.py`, `test_plots.py`, and `test_run_regional_analysis.py` -- 178 tests total.

### Fixed
- This closes the gap flagged in the "possible improvements" memo: regional stations previously had no way to detect missing years within a nominal span, and Mann-Kendall's Sen's-slope magnitude was only ever reported per array-index rather than per calendar year. It's also a prerequisite for any future inter-site correlation work (deferred, see [0.4.4]/conversation), since estimating correlation between pooled stations requires knowing which of their years actually overlap.

## [0.4.5] - 2026-07-27

### Added — GitHub Actions CI
- `.github/workflows/tests.yml`: runs the full pytest suite (`uv run pytest Module/tests/ -v`) on every push and pull request to `main`, plus a manual `workflow_dispatch` trigger. Installs via `uv sync --locked` against the committed `uv.lock` for a reproducible environment (Python 3.12, matching `pyproject.toml`'s `requires-python`). Verified locally against the exact same commands/lockfile before committing: `uv sync --locked` resolves cleanly and all 172 tests pass.
- README: added a status badge placeholder (needs the actual `<OWNER>/<REPO>` GitHub path filled in once pushed) and a note in "How to run the tests" pointing to the workflow.
- Fixed a stale version-number drift: `pyproject.toml`'s `version` field had been stuck at `0.3.0` since before the regional-analysis work started, out of sync with `floodfreq.__version__`. Both are now `0.4.5`.

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
