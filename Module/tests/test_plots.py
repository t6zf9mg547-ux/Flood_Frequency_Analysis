"""
Smoke tests for floodfreq.plots: every save_* function should exist,
accept the expected arguments, and produce a non-empty file without
raising. These don't check visual correctness (that needs a human), but
they do catch import errors, missing functions, and broken signatures --
exactly the class of bug that slipped through once already during
development (an edit accidentally deleted save_probability_plot's
function definition while its body became dead code elsewhere; nothing
caught it until the CLI was run by hand). A quick smoke test per function
is cheap insurance against that happening silently again.
"""
import pandas as pd
import pytest
import matplotlib.pyplot as plt
from floodfreq.analysis import FloodFrequencyAnalysis
from floodfreq import plots as P


@pytest.fixture
def fitted_ffa(reference_data):
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years)
    ffa.fit_all()
    return ffa


def _assert_nonempty_file(path):
    assert path.exists()
    assert path.stat().st_size > 0


def test_save_probability_plot(fitted_ffa, tmp_path):
    out = tmp_path / "prob.png"
    P.save_probability_plot(fitted_ffa, out)
    _assert_nonempty_file(out)


def test_save_data_histogram(fitted_ffa, tmp_path):
    best_key, best_method = fitted_ffa.best_fit()
    out = tmp_path / "hist.png"
    P.save_data_histogram(fitted_ffa, out, dist_key=best_key, method=best_method)
    _assert_nonempty_file(out)


def test_save_data_quality_plot(fitted_ffa, tmp_path):
    out = tmp_path / "dq.png"
    P.save_data_quality_plot(fitted_ffa, out)
    _assert_nonempty_file(out)


def test_save_quantile_ci_plot(fitted_ffa, tmp_path):
    best_key, best_method = fitted_ffa.best_fit()
    out = tmp_path / "ci.png"
    P.save_quantile_ci_plot(fitted_ffa, best_key, best_method, out, n_boot=50)
    _assert_nonempty_file(out)


def test_save_moment_ratio_diagram(fitted_ffa, tmp_path):
    out = tmp_path / "lmr.png"
    P.save_moment_ratio_diagram(fitted_ffa, out)
    _assert_nonempty_file(out)


def test_save_dashboard(fitted_ffa, tmp_path):
    best_key, best_method = fitted_ffa.best_fit()
    out = tmp_path / "dashboard.png"
    P.save_dashboard(fitted_ffa, out, best_key, best_method, n_boot=50)
    _assert_nonempty_file(out)


def test_save_pdf_report(fitted_ffa, tmp_path):
    best_key, best_method = fitted_ffa.best_fit()
    text = fitted_ffa.generate_recommendation(ci_n_boot=50)
    out = tmp_path / "report.pdf"
    P.save_pdf_report(fitted_ffa, out, best_key, best_method, text, n_boot=50)
    _assert_nonempty_file(out)

    from pypdf import PdfReader
    reader = PdfReader(str(out))
    assert len(reader.pages) >= 5  # text page(s) + dashboard + individual plots


def test_probability_plot_uses_custom_variable_labels(reference_data):
    # confirms the variable_name/units/short_name generalization (added for
    # non-streamflow use, e.g. rainfall) actually reaches the rendered plot
    x, years = reference_data
    ffa = FloodFrequencyAnalysis(x, station_id="test", years=years,
                                  variable_name="Rainfall depth", units="mm", short_name="Rainfall")
    ffa.fit_all()
    fig, ax = plt.subplots()
    P.probability_plot(ffa, ax=ax)
    assert ax.get_ylabel() == "Rainfall depth (mm)"
    assert "Rainfall frequency curve" in ax.get_title()
    plt.close(fig)
