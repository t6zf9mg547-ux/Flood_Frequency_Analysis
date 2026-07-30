"""
Tests for run_climate_adjustment.py: the CLI-override helper, and an
in-process end-to-end run that writes the CSV / summary.txt / plot and checks
the numbers are sane (this is the test that guards against the value/year
column-unpacking swap -- read_series returns (values, years), and getting
that order wrong silently produces year-valued "floods").
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from run_climate_adjustment import _override, main


def test_override_cli_wins():
    assert _override(0.25, 0.30) == 0.25


def test_override_file_value_when_cli_none():
    assert _override(None, 0.30) == 0.30


def test_override_cli_zero_is_respected_not_treated_as_falsy():
    # 0.0 is a legitimate override (e.g. tau1 = 0), must not be read as "unset"
    assert _override(0.0, 0.30) == 0.0


def _make_project(tmp_path, q):
    """Minimal <root>/Data + Module skeleton with a baseline series and a
    climate-inputs file, so project_root_from() resolves tmp_path."""
    (tmp_path / "Module").mkdir()
    (tmp_path / "Data").mkdir()
    (tmp_path / "Data" / "Climate_Adjustment").mkdir()
    # baseline series -- note year column FIRST, to exercise the auto-detect /
    # unpacking path that previously mis-picked the year column as the series
    pd.DataFrame({"year": np.arange(1950, 1950 + q.size), "Q": q}).to_csv(
        tmp_path / "Data" / "DemoDam.csv", index=False)
    pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1", "tau2",
                      "confidence_level", "return_periods", "distribution", "pmf"],
        "value": ["0.30", "0.18", "0.10", "0.18", "95", "100,1000,10000", "gumbel", "1200"],
    }).to_csv(tmp_path / "Data" / "Climate_Adjustment" / "DemoDam.csv", index=False)


def test_end_to_end_writes_outputs_with_sane_values(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    q = rng.gumbel(500, 150, size=45)   # ~hundreds-scale flood series
    _make_project(tmp_path, q)

    # run_climate_adjustment.main() locates the project via __file__ of the
    # run script; point that at the fake project by monkeypatching the module's
    # __file__ so project_root_from walks up to tmp_path.
    import run_climate_adjustment as rca
    fake_script = tmp_path / "Module" / "run_climate_adjustment.py"
    fake_script.write_text("# placeholder")
    monkeypatch.setattr(rca, "__file__", str(fake_script))

    monkeypatch.setattr(sys, "argv", ["run_climate_adjustment.py", "DemoDam"])
    rca.main()

    out_dir = tmp_path / "Output" / "DemoDam" / "Climate_Adjustment"
    plot_dir = tmp_path / "Plot" / "DemoDam" / "Climate_Adjustment"
    csv = out_dir / "climate_adjustment_table.csv"
    summary = out_dir / "summary.txt"
    png = plot_dir / "climate_adjustment.png"
    assert csv.exists() and summary.exists() and png.exists()

    df = pd.read_csv(csv)
    assert list(df["T"]) == [100, 1000, 10000]
    # the flood series is ~hundreds; if the year column had been used by
    # mistake, these would be ~2000. Guard the whole pipeline:
    assert df["baseline_point"].between(100, 3000).all()
    assert df["climate_point"].between(100, 3000).all()
    # climate central > baseline central (mean shifted up), and combined CI
    # is wider than the baseline CI everywhere
    assert (df["climate_point"] > df["baseline_point"]).all()
    base_w = df["baseline_upper"] - df["baseline_lower"]
    clim_w = df["climate_upper"] - df["climate_lower"]
    assert (clim_w > base_w).all()

    text = summary.read_text()
    assert "CLIMATE-INFORMED FLOOD ADJUSTMENT (CIFAM)" in text
    # baseline mean reported in the summary must be the flood series (~hundreds),
    # not the ~1970 year mean
    assert "mean = " in text


def test_end_to_end_cli_override_and_no_plot(tmp_path, monkeypatch):
    rng = np.random.default_rng(1)
    q = rng.gumbel(500, 150, size=40)
    _make_project(tmp_path, q)
    import run_climate_adjustment as rca
    fake_script = tmp_path / "Module" / "run_climate_adjustment.py"
    fake_script.write_text("# placeholder")
    monkeypatch.setattr(rca, "__file__", str(fake_script))

    # override delta1 on the CLI and skip the plot
    monkeypatch.setattr(sys, "argv", [
        "run_climate_adjustment.py", "DemoDam",
        "--delta1", "0.50", "--no-plot",
    ])
    rca.main()

    out_dir = tmp_path / "Output" / "DemoDam" / "Climate_Adjustment"
    plot_dir = tmp_path / "Plot" / "DemoDam" / "Climate_Adjustment"
    assert (out_dir / "climate_adjustment_table.csv").exists()
    # --no-plot: no PNG written
    assert not (plot_dir / "climate_adjustment.png").exists()
    # summary should record the overridden delta1 = 0.5
    assert "0.500" in (out_dir / "summary.txt").read_text()
