import pandas as pd
import pytest
from floodfreq.io_utils import (
    read_series, load_case_config, resolve_region, load_region_stations,
    resolve_climate_case, load_climate_inputs,
)


def test_reads_valid_csv(tmp_path):
    p = tmp_path / "station.csv"
    pd.DataFrame({"year": [2000, 2001, 2002, 2003, 2004],
                  "Q": [100.0, 110.0, 105.0, 120.0, 95.0]}).to_csv(p, index=False)
    values, years = read_series(p)
    assert len(values) == 5
    assert years is not None and len(years) == 5


def test_missing_file_raises_filenotfounderror(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_series(tmp_path / "does_not_exist.csv")


def test_empty_file_raises_clear_error(tmp_path):
    p = tmp_path / "empty.csv"
    p.write_text("year,Q\n")  # header only, no rows
    with pytest.raises(ValueError, match="no data|no rows"):
        read_series(p)


def test_no_numeric_column_raises(tmp_path):
    p = tmp_path / "text_only.csv"
    pd.DataFrame({"name": ["a", "b", "c"], "note": ["x", "y", "z"]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="[Nn]o numeric column"):
        read_series(p)


def test_explicit_value_col_not_found_raises(tmp_path):
    p = tmp_path / "station.csv"
    pd.DataFrame({"year": [2000, 2001, 2002], "Q": [100.0, 110.0, 105.0]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="value_col"):
        read_series(p, value_col="does_not_exist")


def test_all_missing_values_raises(tmp_path):
    p = tmp_path / "all_nan.csv"
    pd.DataFrame({"year": [2000, 2001, 2002], "Q": [None, None, None]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="[Nn]o data remains|dropping missing"):
        read_series(p)


def test_year_column_autodetected(tmp_path):
    p = tmp_path / "station.csv"
    pd.DataFrame({"year": [1990, 1991, 1992, 1993],
                  "flow_cms": [50.0, 55.0, 48.0, 60.0]}).to_csv(p, index=False)
    values, years = read_series(p)
    assert years is not None
    assert years.min() == 1990


# -- load_case_config -- #

def test_load_case_config_missing_file_returns_empty(tmp_path):
    (tmp_path / "Data").mkdir()
    assert load_case_config("NoSuchCase", tmp_path) == {}


def test_load_case_config_reads_valid_toml(tmp_path):
    (tmp_path / "Data" / "MyCase").mkdir(parents=True)
    (tmp_path / "Data" / "MyCase" / "MyCase.toml").write_text(
        'confidence_level = 90\nregional_skew = 0.0\npdf_report = true\n')
    config = load_case_config("MyCase", tmp_path)
    assert config == {"confidence_level": 90, "regional_skew": 0.0, "pdf_report": True}


def test_load_case_config_malformed_toml_raises_clear_error(tmp_path):
    (tmp_path / "Data" / "Broken").mkdir(parents=True)
    (tmp_path / "Data" / "Broken" / "Broken.toml").write_text("this is not valid toml {{{")
    with pytest.raises(ValueError, match="Broken.toml"):
        load_case_config("Broken", tmp_path)


# -- resolve_region / load_region_stations -- #

def _make_fake_project(tmp_path):
    """Minimal <root>/Data, <root>/Module skeleton so project_root_from()
    can find tmp_path as the project root."""
    (tmp_path / "Data").mkdir()
    (tmp_path / "Module").mkdir()
    return tmp_path / "Module" / "run_regional_analysis.py"  # need not exist


def test_resolve_region_builds_expected_paths_and_creates_output_dirs(tmp_path):
    module_file = _make_fake_project(tmp_path)
    paths = resolve_region("Template", module_file)
    assert paths.region_name == "Template"
    assert paths.project_root == tmp_path
    assert paths.data_dir == tmp_path / "Data" / "Regional" / "Template"
    assert paths.output_dir == tmp_path / "Output" / "Regional" / "Template"
    assert paths.plot_dir == tmp_path / "Plot" / "Regional" / "Template"
    # output/plot dirs are created eagerly (mirrors resolve_case's behavior)
    assert paths.output_dir.is_dir()
    assert paths.plot_dir.is_dir()
    # the data dir is NOT created -- the station CSVs must already exist there
    assert not paths.data_dir.exists()


def test_load_region_stations_reads_one_array_per_csv(tmp_path):
    module_file = _make_fake_project(tmp_path)
    paths = resolve_region("Template", module_file)
    paths.data_dir.mkdir(parents=True)
    pd.DataFrame({"year": [2000, 2001, 2002, 2003, 2004],
                  "Q": [100.0, 110.0, 105.0, 120.0, 95.0]}).to_csv(
        paths.data_dir / "STN_A.csv", index=False)
    pd.DataFrame({"year": [1990, 1991, 1992, 1993, 1994, 1995],
                  "Q": [50.0, 55.0, 48.0, 60.0, 52.0, 58.0]}).to_csv(
        paths.data_dir / "STN_B.csv", index=False)

    station_data, station_years = load_region_stations(paths)
    assert set(station_data) == {"STN_A", "STN_B"}
    assert len(station_data["STN_A"]) == 5
    assert len(station_data["STN_B"]) == 6
    assert set(station_years) == {"STN_A", "STN_B"}
    assert list(station_years["STN_A"]) == [2000, 2001, 2002, 2003, 2004]
    assert list(station_years["STN_B"]) == [1990, 1991, 1992, 1993, 1994, 1995]


def test_load_region_stations_year_is_none_without_year_column(tmp_path):
    module_file = _make_fake_project(tmp_path)
    paths = resolve_region("Template", module_file)
    paths.data_dir.mkdir(parents=True)
    # no year-like column at all -- just an unlabeled value series
    pd.DataFrame({"flood": [100.0, 110.0, 105.0, 120.0, 95.0]}).to_csv(
        paths.data_dir / "STN_C.csv", index=False)

    station_data, station_years = load_region_stations(paths)
    assert len(station_data["STN_C"]) == 5
    assert station_years["STN_C"] is None


def test_load_region_stations_missing_dir_raises(tmp_path):
    module_file = _make_fake_project(tmp_path)
    paths = resolve_region("NoSuchRegion", module_file)
    with pytest.raises(FileNotFoundError, match="Region data folder not found"):
        load_region_stations(paths)


def test_load_region_stations_empty_dir_raises(tmp_path):
    module_file = _make_fake_project(tmp_path)
    paths = resolve_region("EmptyRegion", module_file)
    paths.data_dir.mkdir(parents=True)
    with pytest.raises(ValueError, match="No station CSVs found"):
        load_region_stations(paths)

# -- resolve_climate_case / load_climate_inputs (CIFAM) -- #

def _climate_template_rows():
    """The canonical long-format rows, mirroring Data/Climate_Adjustment/Template.csv."""
    return pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1", "tau2", "confidence_level",
                      "return_periods", "distribution", "pmf"],
        "value": ["0.30", "0.18", "0.10", "0.18", "95",
                  "2,10,100,1000,10000", "gumbel", ""],
        "units": ["fraction", "fraction", "fraction", "fraction", "percent",
                  "years", "name", "m3/s"],
        "description": ["mean shift", "mean shift sd", "sd shift", "sd shift sd",
                        "CI level", "return periods", "dist", "optional PMF"],
    })


def test_resolve_climate_case_builds_expected_paths_and_creates_output_dirs(tmp_path):
    module_file = _make_fake_project(tmp_path)
    paths = resolve_climate_case("LomPangar", "rcp85", module_file)
    assert paths.case_name == "LomPangar"
    assert paths.scenario == "rcp85"
    assert paths.project_root == tmp_path
    # baseline series reused from the case-first single-station location
    assert paths.baseline_csv == tmp_path / "Data" / "LomPangar" / "LomPangar.csv"
    # the four climate numbers live in the case's climate_adjustment/ subfolder,
    # one file per scenario
    assert paths.climate_csv == (tmp_path / "Data" / "LomPangar"
                                 / "climate_adjustment" / "rcp85.csv")
    # outputs are case-first AND scenario-nested (scenarios don't collide)
    assert paths.output_dir == tmp_path / "Output" / "LomPangar" / "Climate_Adjustment" / "rcp85"
    assert paths.plot_dir == tmp_path / "Plot" / "LomPangar" / "Climate_Adjustment" / "rcp85"
    # output/plot dirs created eagerly (mirrors resolve_case / resolve_region)
    assert paths.output_dir.is_dir()
    assert paths.plot_dir.is_dir()
    # neither input file is created by the resolver -- they must be supplied
    assert not paths.baseline_csv.exists()
    assert not paths.climate_csv.exists()


def test_resolve_climate_case_scenarios_do_not_collide(tmp_path):
    module_file = _make_fake_project(tmp_path)
    p45 = resolve_climate_case("Dam", "rcp45", module_file)
    p85 = resolve_climate_case("Dam", "rcp85", module_file)
    assert p45.output_dir != p85.output_dir
    assert p45.climate_csv != p85.climate_csv
    assert p45.baseline_csv == p85.baseline_csv  # same baseline, different scenarios


def test_load_climate_inputs_parses_all_fields(tmp_path):
    p = tmp_path / "case.csv"
    _climate_template_rows().to_csv(p, index=False)
    ci = load_climate_inputs(p)
    assert ci["delta1"] == pytest.approx(0.30)
    assert ci["delta2"] == pytest.approx(0.18)
    assert ci["tau1"] == pytest.approx(0.10)
    assert ci["tau2"] == pytest.approx(0.18)
    assert ci["confidence_level"] == pytest.approx(95.0)
    assert ci["distribution"] == "gumbel"
    assert ci["return_periods"] == (2.0, 10.0, 100.0, 1000.0, 10000.0)
    assert ci["pmf"] is None


def test_load_climate_inputs_optional_defaults_when_only_required_given(tmp_path):
    p = tmp_path / "minimal.csv"
    pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1", "tau2"],
        "value": [0.30, 0.18, 0.10, 0.18],
    }).to_csv(p, index=False)
    ci = load_climate_inputs(p)
    assert ci["confidence_level"] == 95.0
    assert ci["distribution"] == "gumbel"
    assert ci["pmf"] is None
    assert ci["return_periods"][0] == 2 and ci["return_periods"][-1] == 100000


def test_load_climate_inputs_reads_pmf_when_present(tmp_path):
    p = tmp_path / "with_pmf.csv"
    pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1", "tau2", "pmf"],
        "value": [0.30, 0.18, 0.10, 0.18, 4140],
    }).to_csv(p, index=False)
    ci = load_climate_inputs(p)
    assert ci["pmf"] == pytest.approx(4140.0)


def test_load_climate_inputs_percent_mistake_raises(tmp_path):
    # entering 30 instead of 0.30 is the classic mistake -- must be caught
    p = tmp_path / "percent_mistake.csv"
    pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1", "tau2"],
        "value": [30, 18, 10, 18],
    }).to_csv(p, index=False)
    with pytest.raises(ValueError, match="looks like a percent"):
        load_climate_inputs(p)


def test_load_climate_inputs_missing_required_raises(tmp_path):
    p = tmp_path / "missing.csv"
    pd.DataFrame({
        "parameter": ["delta1", "delta2", "tau1"],  # tau2 missing
        "value": [0.30, 0.18, 0.10],
    }).to_csv(p, index=False)
    with pytest.raises(ValueError, match="tau2"):
        load_climate_inputs(p)


def test_load_climate_inputs_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_climate_inputs(tmp_path / "nope.csv")


def test_load_climate_inputs_wrong_columns_raises(tmp_path):
    p = tmp_path / "wrong.csv"
    pd.DataFrame({"name": ["delta1"], "val": [0.3]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="parameter.*value|value.*column"):
        load_climate_inputs(p)
