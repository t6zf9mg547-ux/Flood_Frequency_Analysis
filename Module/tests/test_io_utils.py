import pandas as pd
import pytest
from floodfreq.io_utils import read_series, load_case_config, resolve_region, load_region_stations


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
    (tmp_path / "Data").mkdir()
    (tmp_path / "Data" / "MyCase.toml").write_text(
        'confidence_level = 90\nregional_skew = 0.0\npdf_report = true\n')
    config = load_case_config("MyCase", tmp_path)
    assert config == {"confidence_level": 90, "regional_skew": 0.0, "pdf_report": True}


def test_load_case_config_malformed_toml_raises_clear_error(tmp_path):
    (tmp_path / "Data").mkdir()
    (tmp_path / "Data" / "Broken.toml").write_text("this is not valid toml {{{")
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
