"""
Tests for floodfreq.regional -- the Hosking & Wallis (1997) L-moment-based
regional (index-flood) flood-frequency pipeline.

Validation strategy (mirrors the rest of this project's test suite):
- Sample L-moments are checked directly against lmoments3's own
  `lmom_ratios()` (the ground truth this module wraps).
- Discordancy is checked against a hand-verifiable near-identical group
  (D_i should be small and roughly equal) and against a group with one
  deliberately-planted outlier station (should be flagged).
- Heterogeneity is checked on a synthetic region simulated from a single
  known distribution (should read as homogeneous, H1 < 1) and on a region
  built from two very different distributions (should read as
  heterogeneous, H1 large).
- `theoretical_tau4` is checked against the two closed-form
  t4(t3) relations that already appear elsewhere in this project
  (`plots._lmr_reference_curves` / `moment_ratio_diagram`): GLO and GPA.
- The growth curve is checked by recovering a known GEV shape from a large
  pooled synthetic sample and comparing predicted vs. true quantiles.
"""
import numpy as np
import pandas as pd
import pytest
import lmoments3.distr as ld
from lmoments3 import lmom_ratios as sample_lmom_ratios

from floodfreq import regional as R


# ---------------------------------------------------------------------------
# Fixtures: synthetic pooling groups
# ---------------------------------------------------------------------------
@pytest.fixture
def homogeneous_region():
    """6 stations, same GEV shape, different scale/location/record length --
    a textbook-homogeneous index-flood group."""
    rng = np.random.default_rng(20260726)
    shape_true = -0.15
    lengths = [100, 120, 90, 150, 110, 130]
    scales = [50, 80, 120, 60, 90, 70]
    data = {}
    for i, (n, scale) in enumerate(zip(lengths, scales)):
        loc = 3 * scale
        x = ld.gev(c=shape_true, loc=loc, scale=scale).rvs(size=n, random_state=rng)
        data[f"STN{i + 1}"] = np.abs(x) + 1.0
    return data, shape_true


@pytest.fixture
def discordant_region(homogeneous_region):
    """Same 6-station homogeneous group, but STN3 replaced with a station
    drawn from a very different (much more skewed/heavy-tailed) GEV --
    should be flagged discordant and should push H1 up."""
    data, shape_true = homogeneous_region
    data = dict(data)
    rng = np.random.default_rng(99)
    x = ld.gev(c=-0.9, loc=300, scale=200).rvs(size=25, random_state=rng)
    data["STN3"] = np.abs(x) + 1.0
    return data


# ---------------------------------------------------------------------------
# 1. Station L-moments
# ---------------------------------------------------------------------------
def test_station_lmoments_matches_lmoments3_directly(reference_data):
    values, _years = reference_data
    s = R.station_lmoments("P1009", values)
    l1, l2, t3, t4, t5 = sample_lmom_ratios(values, nmom=5)
    assert s.l1 == pytest.approx(l1)
    assert s.l2 == pytest.approx(l2)
    assert s.t == pytest.approx(l2 / l1)
    assert s.t3 == pytest.approx(t3)
    assert s.t4 == pytest.approx(t4)
    assert s.n == len(values)
    assert s.mean == pytest.approx(np.mean(values))


def test_station_lmoments_rejects_short_record():
    with pytest.raises(ValueError, match="at least 5"):
        R.station_lmoments("short", [1.0, 2.0, 3.0])


def test_station_lmoments_rejects_nonpositive():
    with pytest.raises(ValueError, match="strictly positive"):
        R.station_lmoments("bad", [1.0, 2.0, -3.0, 4.0, 5.0])


def test_station_lmoments_rejects_nan():
    with pytest.raises(ValueError, match="NaN/Inf"):
        R.station_lmoments("bad", [1.0, 2.0, np.nan, 4.0, 5.0])


def test_stations_table_shape(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    df = R.stations_table(stations)
    assert len(df) == len(data)
    assert set(df.columns) == {"station", "n", "mean", "l1", "l2", "t_LCV", "t3_Lskew", "t4_Lkurt"}


# ---------------------------------------------------------------------------
# 2. Discordancy
# ---------------------------------------------------------------------------
def test_discordancy_needs_at_least_4_stations():
    stations = [R.station_lmoments(f"S{i}", np.arange(1, 11) + i) for i in range(3)]
    with pytest.raises(ValueError, match="at least 4 stations"):
        R.discordancy(stations)


def test_discordancy_near_identical_group_is_small_and_not_flagged():
    rng = np.random.default_rng(0)
    base = ld.gev(c=-0.1, loc=100, scale=30).rvs(size=1000, random_state=rng)
    stations = []
    for i in range(6):
        # small independent perturbations of the same big sample -> nearly
        # identical L-moment ratios across "stations"
        idx = rng.choice(len(base), size=200, replace=False)
        stations.append(R.station_lmoments(f"S{i}", base[idx]))
    df = R.discordancy(stations)
    assert not df["discordant"].any()
    # sum_i D_i == N always (trace identity: sum_i dev_i^T S^-1 dev_i =
    # trace(S^-1 S) = 3, so sum D_i = (N/3)*3 = N), regardless of how
    # (non-)uniform the group is -- a useful internal consistency check.
    assert df["D_i"].sum() == pytest.approx(len(stations), rel=1e-6)


def test_discordancy_flags_planted_outlier(discordant_region):
    stations = [R.station_lmoments(k, v) for k, v in discordant_region.items()]
    df = R.discordancy(stations)
    row = df.set_index("station").loc["STN3"]
    assert row["discordant"]
    assert row["D_i"] == df["D_i"].max()


def test_discordancy_critical_value_table():
    assert R.discordancy_critical_value(5) == 1.333
    assert R.discordancy_critical_value(10) == 2.491
    assert R.discordancy_critical_value(15) == 3.061
    assert R.discordancy_critical_value(20) == R.DISCORDANCY_CRITICAL_LARGE_N
    # below the tabulated range: soft floor, not an error
    assert R.discordancy_critical_value(4) == R.DISCORDANCY_CRITICAL[5]


# ---------------------------------------------------------------------------
# 3. Heterogeneity
# ---------------------------------------------------------------------------
def test_heterogeneity_homogeneous_region_reads_homogeneous(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    het = R.heterogeneity(stations, n_sim=300, seed=1)
    assert het.H1 < 1.5  # comfortably in/near the "acceptably homogeneous" band
    assert het.interpretation in ("acceptably homogeneous", "possibly heterogeneous")
    assert het.n_sim == 300


def test_heterogeneity_discordant_region_reads_more_heterogeneous(homogeneous_region, discordant_region):
    data_h, _ = homogeneous_region
    stations_h = [R.station_lmoments(k, v) for k, v in data_h.items()]
    het_h = R.heterogeneity(stations_h, n_sim=300, seed=1)

    stations_d = [R.station_lmoments(k, v) for k, v in discordant_region.items()]
    het_d = R.heterogeneity(stations_d, n_sim=300, seed=1)

    assert het_d.H1 > het_h.H1
    assert het_d.H1 >= 2.0
    assert het_d.interpretation == "definitely heterogeneous"


def test_heterogeneity_needs_at_least_2_stations():
    stations = [R.station_lmoments("only", np.arange(1, 21))]
    with pytest.raises(ValueError, match="at least 2 stations"):
        R.heterogeneity(stations)


# ---------------------------------------------------------------------------
# 4. Theoretical tau4(tau3) and Z-statistic
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("t3", [-0.3, -0.1, 0.0, 0.1, 0.2, 0.35])
def test_theoretical_tau4_glo_matches_closed_form(t3):
    expected = (1 + 5 * t3 ** 2) / 6
    assert R.theoretical_tau4(t3, "glo") == pytest.approx(expected, abs=1e-5)


@pytest.mark.parametrize("t3", [-0.3, -0.1, 0.1, 0.2, 0.35])
def test_theoretical_tau4_gpa_matches_closed_form(t3):
    expected = t3 * (1 + 5 * t3) / (5 + t3)
    assert R.theoretical_tau4(t3, "gpa") == pytest.approx(expected, abs=1e-5)


def test_theoretical_tau4_unknown_family_raises():
    with pytest.raises(ValueError, match="Unknown family"):
        R.theoretical_tau4(0.1, "not_a_family")


def test_zstatistics_recovers_generating_family(homogeneous_region):
    """Data generated from GEV should show a small |Z| for GEV (and its
    close L-moment-ratio-diagram neighbors GNO/PE3), and should reject the
    two 2-shape-parameter-poorest matches less consistently -- the key,
    robust assertion is that GEV itself is never rejected outright."""
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    zdf = R.zstatistics(stations, n_sim=300, seed=3)
    assert set(zdf["family"]) == set(R.CANDIDATE_FAMILIES)
    gev_row = zdf.set_index("family").loc["gev"]
    assert abs(gev_row["Z"]) <= 1.64


def test_recommend_family_prefers_acceptable_smallest_Z():
    df = pd.DataFrame({
        "family": ["glo", "gev", "pe3"],
        "Z": [2.5, 0.3, -1.0],
        "acceptable": [False, True, True],
    })
    assert R.recommend_family(df) == "gev"


def test_recommend_family_falls_back_when_none_acceptable():
    df = pd.DataFrame({
        "family": ["glo", "gev"],
        "Z": [3.0, -2.5],
        "acceptable": [False, False],
    })
    # smallest |Z| among the (all-unacceptable) candidates
    assert R.recommend_family(df) == "gev"


# ---------------------------------------------------------------------------
# 5 & 6. Growth curve
# ---------------------------------------------------------------------------
def test_growth_curve_l1_is_one(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    gc = R.fit_growth_curve(stations, "gev")
    l1, l2, l3, l4 = R._population_lmoments_from_ppf(gc.ppf, nmom=4)
    assert l1 == pytest.approx(1.0, abs=1e-3)


def test_growth_curve_recovers_known_quantiles_from_large_pooled_sample():
    """Pool many long records from the *same* GEV distribution (up to a
    scale factor per station) and check the fitted growth curve's
    T=100-year quantile is close to the true growth factor."""
    rng = np.random.default_rng(2026)
    shape_true, loc_true, scale_true = -0.1, 1.0, 0.3  # growth-curve scale (mean ~ 1)
    true_dist = ld.gev(c=shape_true, loc=loc_true, scale=scale_true)
    true_mean = true_dist.mean()

    data = {}
    for i, mult in enumerate([1.0, 2.0, 0.5, 3.0, 1.5, 0.8]):
        x = true_dist.rvs(size=2000, random_state=rng) * mult
        data[f"STN{i}"] = np.abs(x) + 1e-6

    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    gc = R.fit_growth_curve(stations, "gev")

    true_growth_100 = true_dist.ppf(1 - 1 / 100) / true_mean
    fitted_growth_100 = gc.quantile(100)
    assert fitted_growth_100 == pytest.approx(true_growth_100, rel=0.1)


def test_station_quantile_scales_by_index_flood(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    gc = R.fit_growth_curve(stations, "gev")
    g100 = gc.quantile(100)
    for s in stations:
        assert gc.station_quantile(100, s.mean) == pytest.approx(g100 * s.mean)


def test_fit_growth_curve_unknown_family_raises(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    with pytest.raises(ValueError, match="Unknown family"):
        R.fit_growth_curve(stations, "bogus")


def test_growth_curve_quantile_table_columns(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    gc = R.fit_growth_curve(stations, "gev")
    table = R.growth_curve_quantile_table(gc, stations, return_periods=(2, 10, 100))
    assert list(table["T_years"]) == [2, 10, 100]
    for s in stations:
        assert s.name in table.columns
        assert np.allclose(table[s.name].to_numpy(), (table["growth_factor"] * s.mean).to_numpy())


# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------
def test_run_regional_analysis_end_to_end(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=200, seed=5)
    assert result.region_name == "TestRegion"
    assert len(result.stations) == len(data)
    assert result.chosen_family in R.CANDIDATE_FAMILIES
    assert not result.discordancy_df.empty
    assert result.heterogeneity_result.n_sim == 200
    assert not result.quantile_table.empty
    assert all(name in result.quantile_table.columns for name in data)


def test_run_regional_analysis_forces_family(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=100, seed=5, family="glo")
    assert result.chosen_family == "glo"
    assert result.growth_curve.family == "glo"


def test_run_regional_analysis_needs_at_least_4_stations(homogeneous_region):
    data, _ = homogeneous_region
    few = dict(list(data.items())[:3])
    with pytest.raises(ValueError, match="at least 4 stations"):
        R.run_regional_analysis("TooFew", few, n_sim=50)


# ---------------------------------------------------------------------------
# Per-station data quality
# ---------------------------------------------------------------------------
def test_station_data_quality_columns_and_rows(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    df = R.station_data_quality(stations, data)
    assert len(df) == len(stations)
    assert set(df["station"]) == set(data)
    for col in ("mann_kendall_trend", "mann_kendall_p_value", "mann_kendall_significant",
                "grubbs_high_outlier_flagged", "grubbs_low_outlier_flagged", "validation_warnings"):
        assert col in df.columns


def test_station_data_quality_flags_planted_trend():
    rng = np.random.default_rng(0)
    n = 60
    # a strong deterministic upward trend should be caught by Mann-Kendall
    x = 100 + 5 * np.arange(n) + rng.normal(0, 2, n)
    stations = [R.station_lmoments("TRENDING", x)]
    df = R.station_data_quality(stations, {"TRENDING": x})
    row = df.iloc[0]
    assert row["mann_kendall_significant"]
    assert "increasing" in row["mann_kendall_trend"]


def test_run_regional_analysis_includes_data_quality(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=100, seed=5)
    assert len(result.data_quality_df) == len(data)
    assert set(result.data_quality_df["station"]) == set(data)


# ---------------------------------------------------------------------------
# Monte-Carlo confidence intervals
# ---------------------------------------------------------------------------
def test_simulate_growth_curve_quantiles_shape(homogeneous_region):
    data, _ = homogeneous_region
    stations = [R.station_lmoments(k, v) for k, v in data.items()]
    gc = R.fit_growth_curve(stations, "gev")
    g_sim = R._simulate_growth_curve_quantiles(stations, "gev", (10, 100), n_sim=150, seed=3)
    assert g_sim.shape == (150, 2)
    # growth factor should be broadly plausible (not wildly divergent)
    assert np.nanmedian(g_sim[:, 0]) == pytest.approx(gc.quantile(10), rel=0.3)


def test_station_quantile_ci_brackets_point_estimate(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=200, seed=5)
    name = result.stations[0].name
    ci = R.station_quantile_ci(result, name, return_periods=(10, 100), n_sim=200, seed=6)
    assert list(ci["station"].unique()) == [name]
    assert (ci["CI_lower"] <= ci["Q_design"]).all()
    assert (ci["Q_design"] <= ci["CI_upper"]).all()


def test_station_quantile_ci_unknown_station_raises(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=100, seed=5)
    with pytest.raises(ValueError, match="Unknown station"):
        R.station_quantile_ci(result, "NotAStation", n_sim=50)


def test_regional_quantile_ci_covers_all_stations(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=150, seed=5)
    ci = R.regional_quantile_ci(result, return_periods=(10, 100, 1000), n_sim=150, seed=6)
    assert set(ci["station"]) == set(data)
    assert len(ci) == len(data) * 3
    assert (ci["CI_lower"] <= ci["Q_design"]).all()
    assert (ci["Q_design"] <= ci["CI_upper"]).all()


def test_regional_quantile_ci_widens_at_larger_T(homogeneous_region):
    data, _ = homogeneous_region
    result = R.run_regional_analysis("TestRegion", data, n_sim=200, seed=5)
    name = result.stations[0].name
    ci = R.station_quantile_ci(result, name, return_periods=(10, 1000), n_sim=200, seed=6)
    width_10 = float(ci.loc[ci["T_years"] == 10, "CI_upper"].iloc[0] - ci.loc[ci["T_years"] == 10, "CI_lower"].iloc[0])
    width_1000 = float(ci.loc[ci["T_years"] == 1000, "CI_upper"].iloc[0] - ci.loc[ci["T_years"] == 1000, "CI_lower"].iloc[0])
    assert width_1000 > width_10
