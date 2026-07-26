import pandas as pd
import pytest
from floodfreq.io_utils import read_series, load_case_config


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
