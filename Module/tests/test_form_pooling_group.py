"""
Integration tests for form_pooling_group.py: invokes the actual CLI via
subprocess (the most faithful way to test an argparse entry point without
refactoring it), against a synthetic candidate catalog + station-CSV pool
written to a temp directory. Uses a distinctive, unlikely-to-collide
region name and cleans up the Output/Regional/<name> and
Data/Regional/<name> folders it creates in the real project tree
afterward, since resolve_region() (like the rest of this project's CLIs)
resolves paths relative to the actual project root, not an injectable one.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

MODULE_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = MODULE_DIR.parent
SCRIPT = MODULE_DIR / "form_pooling_group.py"
TEST_REGION = "_pytest_pooling_test_region"


@pytest.fixture
def candidate_pool(tmp_path):
    """A synthetic candidate catalog + one raw annual-max CSV per
    candidate, written to tmp_path (outside the real project tree)."""
    rng = np.random.default_rng(0)
    n = 10
    catalog = pd.DataFrame({
        "station": [f"CAND{i}" for i in range(n)],
        "area_km2": rng.uniform(50, 500, n),
        "precip_mm": rng.uniform(800, 1400, n),
        "n_years": rng.integers(20, 40, n),
    })
    catalog_path = tmp_path / "candidates.csv"
    catalog.to_csv(catalog_path, index=False)

    station_dir = tmp_path / "all_stations"
    station_dir.mkdir()
    for _, row in catalog.iterrows():
        years = np.arange(1980, 1980 + int(row["n_years"]))
        q = np.abs(rng.normal(100, 20, len(years))) + 1
        pd.DataFrame({"year": years, "Q": q}).to_csv(station_dir / f"{row['station']}.csv",
                                                      index=False)
    return catalog_path, station_dir, catalog


@pytest.fixture(autouse=True)
def _cleanup_region_dirs():
    """Remove any Data/Output/Plot/Regional/<TEST_REGION> folders this
    test file's subprocess calls create in the real project tree, before
    and after each test."""
    import shutil

    def _clean():
        for base in ("Data", "Output", "Plot"):
            d = PROJECT_ROOT / base / "Regional" / TEST_REGION
            if d.exists():
                shutil.rmtree(d)

    _clean()
    yield
    _clean()


def _run_cli(args):
    return subprocess.run([sys.executable, str(SCRIPT)] + args,
                           capture_output=True, text=True, cwd=MODULE_DIR)


def test_dry_run_reports_ranking_and_proposal(candidate_pool):
    catalog_path, station_dir, catalog = candidate_pool
    target = catalog.iloc[0]["station"]

    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2", "precip_mm",
        "--target-station", target,
        "--n-stations", "3",
        "--region-name", TEST_REGION,
    ])
    assert result.returncode == 0, result.stderr
    assert "Proposed pooling group (3 station(s))" in result.stdout
    assert "Dry run only" in result.stdout

    ranking_csv = PROJECT_ROOT / "Output" / "Regional" / TEST_REGION / "pooling_ranking_full.csv"
    proposed_csv = PROJECT_ROOT / "Output" / "Regional" / TEST_REGION / "pooling_group_proposed.csv"
    assert ranking_csv.exists()
    assert proposed_csv.exists()
    proposed = pd.read_csv(proposed_csv)
    assert len(proposed) == 3
    assert target not in set(proposed["station"])  # target excluded from its own ranking

    # dry run must NOT create the Data/ folder
    assert not (PROJECT_ROOT / "Data" / "Regional" / TEST_REGION).exists()


def test_apply_copies_proposed_stations_and_is_usable_by_regional_analysis(candidate_pool):
    catalog_path, station_dir, catalog = candidate_pool
    target = catalog.iloc[0]["station"]

    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2", "precip_mm",
        "--target-station", target,
        "--n-stations", "4",
        "--region-name", TEST_REGION,
        "--station-data-dir", str(station_dir),
        "--apply",
    ])
    assert result.returncode == 0, result.stderr
    assert "Copied 4 station CSV(s)" in result.stdout

    data_dir = PROJECT_ROOT / "Data" / "Regional" / TEST_REGION
    assert data_dir.is_dir()
    copied = sorted(p.stem for p in data_dir.glob("*.csv"))
    proposed = pd.read_csv(
        PROJECT_ROOT / "Output" / "Regional" / TEST_REGION / "pooling_group_proposed.csv")
    assert copied == sorted(proposed["station"])

    # the copied group should be directly usable by run_regional_analysis.py
    analysis = subprocess.run(
        [sys.executable, str(MODULE_DIR / "run_regional_analysis.py"), TEST_REGION,
         "--n-sim", "50", "--no-ci", "--no-plots"],
        capture_output=True, text=True, cwd=MODULE_DIR)
    assert analysis.returncode == 0, analysis.stderr
    assert "Done." in analysis.stdout


def test_min_years_stopping_rule(candidate_pool):
    catalog_path, station_dir, catalog = candidate_pool
    target = catalog.iloc[0]["station"]

    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2", "precip_mm",
        "--target-station", target,
        "--min-years", "60",
        "--region-name", TEST_REGION,
    ])
    assert result.returncode == 0, result.stderr
    proposed = pd.read_csv(
        PROJECT_ROOT / "Output" / "Regional" / TEST_REGION / "pooling_group_proposed.csv")
    assert proposed["n_years"].sum() >= 60


def test_ungauged_target_descriptors(candidate_pool):
    catalog_path, station_dir, catalog = candidate_pool

    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2", "precip_mm",
        "--target-descriptors", "area_km2=200", "precip_mm=1100",
        "--n-stations", "3",
        "--region-name", TEST_REGION,
    ])
    assert result.returncode == 0, result.stderr
    assert "(ungauged target)" in result.stdout


def test_mutually_exclusive_target_flags_rejected(candidate_pool):
    catalog_path, station_dir, catalog = candidate_pool
    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2",
        "--target-station", catalog.iloc[0]["station"],
        "--target-descriptors", "area_km2=1",
        "--n-stations", "2",
        "--region-name", TEST_REGION,
    ])
    assert result.returncode != 0
    assert "exactly one" in result.stdout or "exactly one" in result.stderr


def test_apply_without_station_data_dir_rejected(candidate_pool):
    catalog_path, station_dir, catalog = candidate_pool
    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2",
        "--target-station", catalog.iloc[0]["station"],
        "--n-stations", "2",
        "--region-name", TEST_REGION,
        "--apply",
    ])
    assert result.returncode != 0
    assert "--station-data-dir" in (result.stdout + result.stderr)


def test_apply_warns_about_missing_station_csvs(candidate_pool, tmp_path):
    catalog_path, station_dir, catalog = candidate_pool
    target = catalog.iloc[0]["station"]
    # an empty station-data-dir -- every proposed station's CSV will be "missing"
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", "area_km2", "precip_mm",
        "--target-station", target,
        "--n-stations", "3",
        "--region-name", TEST_REGION,
        "--station-data-dir", str(empty_dir),
        "--apply",
    ])
    assert result.returncode == 0, result.stderr
    assert "Copied 0 station CSV(s)" in result.stdout
    assert "WARNING" in result.stdout


# ---------------------------------------------------------------------------
# Bundled templates: Data/Templates/candidate_descriptors_{discharge,rainfall}_TEMPLATE.csv
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("template_name,descriptors,target", [
    ("candidate_descriptors_discharge_TEMPLATE.csv",
     ["area_km2", "mean_annual_precip_mm", "bfihost", "mean_slope_m_per_km", "urban_extent_frac"],
     "GAUGE_001"),
    ("candidate_descriptors_rainfall_TEMPLATE.csv",
     ["lat", "lon", "elevation_m", "mean_annual_precip_mm"],
     "RAINGAUGE_001"),
])
def test_bundled_template_catalog_works_out_of_the_box(template_name, descriptors, target):
    """The two example catalogs shipped in Data/Templates/ should be
    directly usable with form_pooling_group.py, no editing required --
    this is what someone reads the README and tries first."""
    template_path = PROJECT_ROOT / "Data" / "Templates" / template_name
    assert template_path.exists(), f"Expected bundled template at {template_path}"

    result = _run_cli([
        "--catalog", str(template_path),
        "--descriptors", *descriptors,
        "--target-station", target,
        "--n-stations", "4",
        "--region-name", TEST_REGION,
    ])
    assert result.returncode == 0, result.stderr
    assert "Proposed pooling group (4 station(s))" in result.stdout

    proposed = pd.read_csv(
        PROJECT_ROOT / "Output" / "Regional" / TEST_REGION / "pooling_group_proposed.csv")
    assert len(proposed) == 4
    assert target not in set(proposed["station"])


def test_bundled_templates_are_not_visible_to_single_station_picker():
    """Regression check: run_analysis.py's interactive case picker lists
    Data/*.csv (non-recursive) as candidate single-station cases. The
    pooling-group descriptor catalogs are NOT annual-maximum series --
    read_series() would silently pick a descriptor column (e.g. area_km2)
    as if it were the flood series, with no error and no year detected,
    producing a nonsense "analysis" with no warning. They must live in a
    subfolder (Data/Templates/) so run_analysis.py's glob never sees them."""
    data_dir = PROJECT_ROOT / "Data"
    top_level_csvs = {p.name for p in data_dir.glob("*.csv")}
    assert "candidate_descriptors_discharge_TEMPLATE.csv" not in top_level_csvs
    assert "candidate_descriptors_rainfall_TEMPLATE.csv" not in top_level_csvs
    # and they should still exist, just one level down
    assert (data_dir / "Templates" / "candidate_descriptors_discharge_TEMPLATE.csv").exists()
    assert (data_dir / "Templates" / "candidate_descriptors_rainfall_TEMPLATE.csv").exists()


# ---------------------------------------------------------------------------
# Bundled station data: Data/Templates/{discharge,rainfall}_station_data/
# (the annual-maximum series matching each candidate catalog's stations)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("catalog_name,station_dir_name,descriptors,target,value_col", [
    ("candidate_descriptors_discharge_TEMPLATE.csv", "discharge_station_data",
     ["area_km2", "mean_annual_precip_mm", "bfihost", "mean_slope_m_per_km", "urban_extent_frac"],
     "GAUGE_001", "Q"),
    ("candidate_descriptors_rainfall_TEMPLATE.csv", "rainfall_station_data",
     ["lat", "lon", "elevation_m", "mean_annual_precip_mm"],
     "RAINGAUGE_001", "Rainfall_mm"),
])
def test_bundled_station_data_matches_catalog_and_is_fully_runnable(
        catalog_name, station_dir_name, descriptors, target, value_col):
    """The whole point of these bundled folders: --apply should work using
    ONLY bundled content (no user-supplied files needed at all), and the
    resulting Data/Regional/<TEST_REGION>/ should be directly runnable by
    run_regional_analysis.py end to end."""
    templates_dir = PROJECT_ROOT / "Data" / "Templates"
    catalog_path = templates_dir / catalog_name
    station_dir = templates_dir / station_dir_name
    assert station_dir.is_dir(), f"Expected bundled station data at {station_dir}"

    catalog = pd.read_csv(catalog_path)
    # every candidate in the catalog should have a matching station CSV,
    # with a record length matching that catalog row's n_years
    for _, row in catalog.iterrows():
        station_csv = station_dir / f"{row['station']}.csv"
        assert station_csv.exists(), f"Missing bundled station data for {row['station']}"
        series = pd.read_csv(station_csv)
        assert len(series) == int(row["n_years"]), (
            f"{row['station']}: catalog says n_years={row['n_years']} but the bundled "
            f"CSV has {len(series)} rows")

    # --apply using ONLY bundled files (catalog + station_dir), no fixtures
    result = _run_cli([
        "--catalog", str(catalog_path),
        "--descriptors", *descriptors,
        "--target-station", target,
        "--n-stations", "5",
        "--region-name", TEST_REGION,
        "--station-data-dir", str(station_dir),
        "--apply",
    ])
    assert result.returncode == 0, result.stderr
    assert "Copied 5 station CSV(s)" in result.stdout
    assert "WARNING" not in result.stdout  # every proposed station's CSV should be found

    # and the resulting group should be a real, runnable regional analysis
    analysis = subprocess.run(
        [sys.executable, str(MODULE_DIR / "run_regional_analysis.py"), TEST_REGION,
         "--n-sim", "50", "--no-ci", "--no-plots", "--value-col", value_col],
        capture_output=True, text=True, cwd=MODULE_DIR)
    assert analysis.returncode == 0, analysis.stderr
    assert "Done." in analysis.stdout
