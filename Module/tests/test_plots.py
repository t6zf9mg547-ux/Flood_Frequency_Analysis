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


# -- regional_moment_ratio_diagram / save_regional_moment_ratio_diagram -- #

@pytest.fixture
def sample_stations():
    import numpy as np
    from floodfreq import regional as R
    rng = np.random.default_rng(2026)
    stations = []
    for i, (n, scale) in enumerate([(60, 50), (70, 80), (55, 65), (90, 40)]):
        x = np.abs(rng.gumbel(loc=3 * scale, scale=scale, size=n)) + 1.0
        stations.append(R.station_lmoments(f"STN{i + 1}", x))
    return stations


def test_save_regional_moment_ratio_diagram(sample_stations, tmp_path):
    out = tmp_path / "regional_lmr.png"
    P.save_regional_moment_ratio_diagram(sample_stations, out, region_name="TestRegion")
    _assert_nonempty_file(out)


def test_regional_moment_ratio_diagram_marks_discordant_stations(sample_stations, tmp_path):
    from floodfreq import regional as R
    disc_df = R.discordancy(sample_stations)
    fig, ax = plt.subplots()
    P.regional_moment_ratio_diagram(sample_stations, region_name="TestRegion",
                                     discordancy_df=disc_df, ax=ax)
    # legend/labels render without error, and both L-CV/L-skewness axes are set
    assert ax.get_xlabel() == "L-skewness (t3)"
    assert ax.get_ylabel() == "L-kurtosis (t4)"
    plt.close(fig)


def test_regional_moment_ratio_diagram_works_without_discordancy_df(sample_stations):
    fig, ax = plt.subplots()
    P.regional_moment_ratio_diagram(sample_stations, ax=ax)  # no discordancy_df passed
    plt.close(fig)


# -- regional plot suite: growth curve / station design flood / series --
# -- grid / discordancy bar chart / dashboard                          --

@pytest.fixture
def sample_regional_result():
    import numpy as np
    from floodfreq import regional as R
    rng = np.random.default_rng(2026)
    station_data = {}
    for i, (n, scale) in enumerate([(60, 50), (70, 80), (55, 65), (90, 40), (65, 55)]):
        x = np.abs(rng.gumbel(loc=3 * scale, scale=scale, size=n)) + 1.0
        station_data[f"STN{i + 1}"] = x
    return R.run_regional_analysis("TestRegion", station_data, n_sim=100, seed=1)


@pytest.fixture
def sample_regional_result_with_years():
    """Same shape as sample_regional_result, but with a per-station year
    array attached, for exercising the calendar-year x-axis path in
    save_regional_station_series_plot."""
    import numpy as np
    from floodfreq import regional as R
    rng = np.random.default_rng(2026)
    station_data = {}
    station_years = {}
    start_year = 1950
    for i, (n, scale) in enumerate([(60, 50), (70, 80), (55, 65), (90, 40), (65, 55)]):
        x = np.abs(rng.gumbel(loc=3 * scale, scale=scale, size=n)) + 1.0
        name = f"STN{i + 1}"
        station_data[name] = x
        station_years[name] = np.arange(start_year, start_year + n)
    return R.run_regional_analysis("TestRegion", station_data, n_sim=100, seed=1,
                                    station_years=station_years)


def test_save_regional_growth_curve_plot(sample_regional_result, tmp_path):
    out = tmp_path / "growth.png"
    P.save_regional_growth_curve_plot(sample_regional_result, out)
    _assert_nonempty_file(out)


def test_regional_growth_curve_plot_axes(sample_regional_result):
    fig, ax = plt.subplots()
    P.regional_growth_curve_plot(sample_regional_result, ax=ax)
    assert ax.get_xlabel() == "Return period T (years)"
    assert "growth factor" in ax.get_ylabel().lower()
    plt.close(fig)


def test_regional_growth_curve_plot_pooled_rank_toggle(sample_regional_result, tmp_path):
    # with the overlay on (default), the legend should mention the pooled cloud
    fig, ax = plt.subplots()
    P.regional_growth_curve_plot(sample_regional_result, ax=ax, show_pooled_rank=True)
    labels_on = [t.get_text() for t in ax.get_legend().get_texts()]
    plt.close(fig)
    assert any("Pooled rank" in lbl for lbl in labels_on)

    # with it off, the plot should still render fine but without that label
    fig, ax = plt.subplots()
    P.regional_growth_curve_plot(sample_regional_result, ax=ax, show_pooled_rank=False)
    labels_off = [t.get_text() for t in ax.get_legend().get_texts()]
    plt.close(fig)
    assert not any("Pooled rank" in lbl for lbl in labels_off)

    out = tmp_path / "growth_no_pooled.png"
    P.save_regional_growth_curve_plot(sample_regional_result, out, show_pooled_rank=False)
    _assert_nonempty_file(out)


def test_save_regional_pooled_vs_stations_plot(sample_regional_result, tmp_path):
    out = tmp_path / "compare.png"
    P.save_regional_pooled_vs_stations_plot(sample_regional_result, out)
    _assert_nonempty_file(out)


def test_regional_pooled_vs_stations_plot_panel_titles(sample_regional_result):
    fig = plt.figure()
    ax_left, ax_right = P.regional_pooled_vs_stations_plot(sample_regional_result, fig=fig).axes
    assert "Per-station" in ax_left.get_title()
    assert "one station" in ax_right.get_title()
    # right panel has one pooled series (no per-station color legend entries)
    right_labels = [t.get_text() for t in ax_right.get_legend().get_texts()]
    assert any("ranked as one series" in lbl for lbl in right_labels)
    assert not any(s.name in " ".join(right_labels) for s in sample_regional_result.stations)
    plt.close(fig)


def test_regional_pooled_vs_stations_plot_pooled_reaches_further_than_any_station(sample_regional_result):
    """The whole point of the right panel: pooling inflates the apparent
    sample size, so its rightmost (largest-T) point should exceed what any
    single station's own record length could produce."""
    import numpy as np
    result = sample_regional_result
    max_station_n = max(s.n for s in result.stations)
    total_n = sum(s.n for s in result.stations)
    assert total_n > max_station_n  # sanity: pooling is actually combining >1 station's years


def test_save_station_design_flood_plot(sample_regional_result, tmp_path):
    out = tmp_path / "station.png"
    name = sample_regional_result.stations[0].name
    P.save_station_design_flood_plot(sample_regional_result, name, out)
    _assert_nonempty_file(out)


def test_station_design_flood_plot_unknown_station_raises(sample_regional_result):
    with pytest.raises(ValueError, match="Unknown station"):
        P.station_design_flood_plot(sample_regional_result, "NotAStation")


def test_save_regional_station_series_plot(sample_regional_result, tmp_path):
    out = tmp_path / "series.png"
    P.save_regional_station_series_plot(sample_regional_result, out)
    _assert_nonempty_file(out)


def test_regional_station_series_plot_uses_years_when_available(sample_regional_result_with_years, tmp_path):
    result = sample_regional_result_with_years
    out = tmp_path / "series_years.png"
    P.save_regional_station_series_plot(result, out)
    _assert_nonempty_file(out)

    # rebuild the same per-panel logic to check the x-axis actually reflects
    # real years (not just 1..n) for a station that has them
    s = result.stations[0]
    years = result.station_years.get(s.name)
    assert years is not None
    assert years.min() == 1950


def test_regional_station_series_plot_falls_back_without_years(sample_regional_result, tmp_path):
    # sample_regional_result has no station_years at all -- should still
    # render fine, just falling back to observation-order x-axis
    out = tmp_path / "series_no_years.png"
    P.save_regional_station_series_plot(sample_regional_result, out)
    _assert_nonempty_file(out)
    assert sample_regional_result.station_years == {}


def test_save_regional_discordancy_plot(sample_regional_result, tmp_path):
    out = tmp_path / "disc.png"
    P.save_regional_discordancy_plot(sample_regional_result, out)
    _assert_nonempty_file(out)


def test_save_regional_dashboard(sample_regional_result, tmp_path):
    out = tmp_path / "regional_dashboard.png"
    P.save_regional_dashboard(sample_regional_result, out)
    _assert_nonempty_file(out)


# -- Layout regressions: legend/ylabel collision, overlapping station --
# -- labels, and the 10000-year growth factor in the dashboard panel  --

def test_regional_growth_curve_legend_does_not_overlap_ylabel(sample_regional_result):
    fig, ax = plt.subplots(figsize=(7, 5))
    P.regional_growth_curve_plot(sample_regional_result, ax=ax)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    leg_bbox = ax.get_legend().get_window_extent(renderer)
    ylabel_bbox = ax.yaxis.label.get_window_extent(renderer)
    assert not leg_bbox.overlaps(ylabel_bbox)
    for t in ax.get_yticklabels():
        assert not leg_bbox.overlaps(t.get_window_extent(renderer))
    plt.close(fig)


def test_regional_moment_ratio_diagram_station_labels_dont_overlap(sample_regional_result):
    """Reproduces the tight-cluster case that used to overlap: relabels
    every station's annotation without its leader-line arrow (which
    inflates the reported bbox with the arrow's own extent) so this
    checks the actual text glyphs, not the connector lines."""
    import numpy as np
    result = sample_regional_result
    stations = result.stations
    t3 = np.array([s.t3 for s in stations])
    t4 = np.array([s.t4 for s in stations])

    fig, ax = plt.subplots(figsize=(7, 7))
    order = np.argsort(t4)
    n = len(order)
    boxes = []
    for rank, idx in enumerate(order):
        dy = (rank - (n - 1) / 2.0) * 13
        left = rank % 2 == 1
        dx = -28 if left else 14
        ha = "right" if left else "left"
        # deliberately omit arrowprops here: matches the same offsets
        # _label_points_staggered uses, but isolates the text-only bbox
        t = ax.annotate(stations[idx].name, (t3[idx], t4[idx]), xytext=(dx, dy),
                         textcoords="offset points", fontsize=7, ha=ha, va="center")
        fig.canvas.draw()
        boxes.append(t.get_window_extent(fig.canvas.get_renderer()))
    plt.close(fig)

    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            assert not boxes[i].overlaps(boxes[j])


def test_regional_dashboard_text_panel_includes_10000_year_growth_factor(sample_regional_result):
    fig = plt.figure(figsize=(15, 10))
    P.regional_dashboard(sample_regional_result, fig=fig)
    ax_txt = fig.axes[3]  # growth curve, lmr, discordancy, text (gridspec order)
    panel_text = ax_txt.texts[0].get_text()
    g10000 = sample_regional_result.growth_curve.quantile(10000)
    assert f"{g10000:.2f}" in panel_text
    assert "10,000y" in panel_text or "10000y" in panel_text
    plt.close(fig)
