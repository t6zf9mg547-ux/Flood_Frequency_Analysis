"""
Tests for the setting-resolution helpers in run_analysis.py: CLI flag >
config file > hardcoded default precedence.
"""
import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from run_analysis import resolve_setting, resolve_bool_setting


def test_cli_value_wins_over_config():
    config = {"n_boot": 2000}
    assert resolve_setting(500, config, "n_boot", 1000) == 500


def test_config_value_used_when_cli_not_passed():
    config = {"n_boot": 2000}
    assert resolve_setting(None, config, "n_boot", 1000) == 2000


def test_hardcoded_default_used_when_neither_given():
    assert resolve_setting(None, {}, "n_boot", 1000) == 1000


def test_config_value_of_zero_is_respected_not_treated_as_falsy():
    # regional_skew=0.0 is a legitimate value, not "unset" -- must not be
    # confused with None (this would be a classic falsy-value bug)
    config = {"regional_skew": 0.0}
    assert resolve_setting(None, config, "regional_skew", None) == 0.0


def test_bool_setting_true_from_cli_only():
    assert resolve_bool_setting(True, {}, "pdf_report") is True


def test_bool_setting_true_from_config_only():
    assert resolve_bool_setting(False, {"pdf_report": True}, "pdf_report") is True


def test_bool_setting_false_when_neither_set():
    assert resolve_bool_setting(False, {}, "pdf_report") is False


def test_bool_setting_true_when_both_set():
    assert resolve_bool_setting(True, {"pdf_report": True}, "pdf_report") is True
