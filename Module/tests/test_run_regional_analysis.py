"""
Tests for run_regional_analysis.py: the helper functions directly, plus one
end-to-end run against a synthetic 4-station region written to a temporary
project skeleton (Data/Regional/<Region>/, Output/Regional/<Region>/,
Plot/Regional/<Region>/), mirroring how test_run_analysis.py covers
run_analysis.py's helpers without invoking the interactive CLI.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from run_regional_analysis import build_summary_text, build_provenance_header, _csv_hint
from floodfreq.io_utils import resolve_region, load_region_stations
from floodfreq.regional import run_regional_analysis, regional_quantile_ci


@pytest.fixture
def region_project(tmp_path):
    """A minimal <root>/Data/Regional/TestRegion/ project skeleton with 4
    synthetic, hydrologically-similar stations (shared GEV shape)."""
    import lmoments3.distr as ld
    (tmp_path / "Data").mkdir()
    (tmp_path / "Module").mkdir()
    region_dir = tmp_path / "Data" / "Regional" / "TestRegion"
    region_dir.mkdir(parents=True)

    rng = np.random.default_rng(20260726)
    shape_true = -0.1
    for i, (n, scale) in enumerate([(60, 50), (70, 80), (55, 65), (90, 40)]):
        loc = 3 * scale
        x = ld.gev(c=shape_true, loc=loc, scale=scale).rvs(size=n, random_state=rng)
        x = np.abs(x) + 1.0
        years = np.arange(1950, 1950 + n)
        pd.DataFrame({"year": years, "Q": x}).to_csv(region_dir / f"STN{i + 1}.csv", index=False)

    module_file = tmp_path / "Module" / "run_regional_analysis.py"
    return tmp_path, module_file


def test_build_provenance_header_contains_version_and_command():
    header = build_provenance_header()
    assert "RUN PROVENANCE" in header
    assert "floodfreq" in header
    assert "Python:" in header


def test_csv_hint_reports_row_count(tmp_path):
    p = tmp_path / "STN1.csv"
    pd.DataFrame({"year": [2000, 2001, 2002], "Q": [1.0, 2.0, 3.0]}).to_csv(p, index=False)
    hint = _csv_hint(p)
    assert "3 rows" in hint


def test_end_to_end_run_writes_expected_outputs(region_project):
    root, module_file = region_project
    paths = resolve_region("TestRegion", module_file)
    station_data, station_years = load_region_stations(paths)
    assert len(station_data) == 4
    # the fixture writes a year column for every station
    assert all(y is not None for y in station_years.values())

    result = run_regional_analysis("TestRegion", station_data, n_sim=100, seed=1,
                                    station_years=station_years)

    # -- CSV outputs (mirrors what main() writes) --
    result.stations_df.to_csv(paths.output_dir / "station_lmoments.csv", index=False)
    result.data_quality_df.to_csv(paths.output_dir / "data_quality.csv", index=False)
    result.discordancy_df.to_csv(paths.output_dir / "discordancy.csv", index=False)
    result.zstat_df.to_csv(paths.output_dir / "zstatistics.csv", index=False)
    result.quantile_table.to_csv(paths.output_dir / "quantile_table.csv", index=False)
    ci_df = regional_quantile_ci(result, return_periods=(10, 100, 1000), n_sim=100, seed=2)
    ci_df.to_csv(paths.output_dir / "growth_curve_quantiles_ci.csv", index=False)
    for fname in ("station_lmoments.csv", "data_quality.csv", "discordancy.csv",
                  "zstatistics.csv", "quantile_table.csv", "growth_curve_quantiles_ci.csv"):
        f = paths.output_dir / fname
        assert f.exists() and f.stat().st_size > 0

    # -- data_quality.csv content --
    assert len(result.data_quality_df) == len(station_data)
    assert set(result.data_quality_df["station"]) == set(station_data)
    for col in ("mann_kendall_trend", "mann_kendall_significant", "years_available",
                "sens_slope_per_year", "grubbs_high_outlier_flagged", "grubbs_low_outlier_flagged"):
        assert col in result.data_quality_df.columns
    # years were available for every station (per the fixture), so the
    # per-calendar-year Sen's slope should actually be populated
    assert result.data_quality_df["years_available"].all()

    # -- CI table content --
    assert set(ci_df["station"]) == set(station_data)
    assert (ci_df["CI_lower"] <= ci_df["Q_design"]).all()
    assert (ci_df["Q_design"] <= ci_df["CI_upper"]).all()

    # -- Summary text --
    summary = build_summary_text(result, build_provenance_header(), ci_df=ci_df, confidence_level=95.0)
    assert "TestRegion" in summary
    assert "PER-STATION DATA QUALITY" in summary
    assert "DISCORDANCY" in summary
    assert "HETEROGENEITY" in summary
    assert "REGIONAL DISTRIBUTION SELECTION" in summary
    assert "CONFIDENCE INTERVALS" in summary
    assert result.growth_curve.label in summary

    # -- Summary text without a CI table (e.g. --no-ci) should still build fine --
    summary_no_ci = build_summary_text(result, build_provenance_header())
    assert "CONFIDENCE INTERVALS" not in summary_no_ci

    # -- Plots (mirrors what main() writes) --
    from floodfreq.plots import (
        save_regional_moment_ratio_diagram,
        save_regional_growth_curve_plot,
        save_regional_pooled_vs_stations_plot,
        save_station_design_flood_plot,
        save_regional_station_series_plot,
        save_regional_discordancy_plot,
        save_regional_dashboard,
    )
    save_regional_moment_ratio_diagram(result.stations, paths.plot_dir / "regional_moment_ratio_diagram.png",
                                        region_name="TestRegion", discordancy_df=result.discordancy_df)
    save_regional_growth_curve_plot(result, paths.plot_dir / "regional_growth_curve.png")
    save_regional_pooled_vs_stations_plot(result, paths.plot_dir / "regional_pooled_vs_stations.png")
    save_regional_discordancy_plot(result, paths.plot_dir / "regional_discordancy.png")
    save_regional_station_series_plot(result, paths.plot_dir / "regional_station_series.png")
    for name in station_data:
        save_station_design_flood_plot(result, name, paths.plot_dir / f"station_{name}_design_flood.png")
    save_regional_dashboard(result, paths.plot_dir / "regional_dashboard.png")

    expected_plots = [
        "regional_moment_ratio_diagram.png", "regional_growth_curve.png",
        "regional_pooled_vs_stations.png",
        "regional_discordancy.png", "regional_station_series.png", "regional_dashboard.png",
    ] + [f"station_{name}_design_flood.png" for name in station_data]
    for fname in expected_plots:
        f = paths.plot_dir / fname
        assert f.exists() and f.stat().st_size > 0


def test_forced_family_is_reflected_in_summary(region_project):
    root, module_file = region_project
    paths = resolve_region("TestRegion", module_file)
    station_data, station_years = load_region_stations(paths)

    result = run_regional_analysis("TestRegion", station_data, n_sim=50, seed=1, family="pe3",
                                    station_years=station_years)
    assert result.chosen_family == "pe3"
    summary = build_summary_text(result, build_provenance_header())
    assert "Pearson Type III (PE3)" in summary
