"""
Tests for floodfreq.pooling -- automatic pooling-group formation via the
"region of influence" approach (Burn, 1990): rank candidate stations by
similarity to a target site's descriptors, then propose a pooling group.

Validation strategy:
- similarity_ranking() is checked against a hand-verifiable case: a
  candidate constructed to be a near-twin of the target should rank first
  with a small distance, and distances should increase monotonically with
  how "off" a candidate's descriptors are.
- Standardization is checked directly: with only two candidates whose
  descriptors differ by a known amount, the z-scored distance should
  match a hand-computed value.
- propose_pooling_group()'s two stopping rules (n_stations, min_total_years)
  are checked against hand-computed expectations, including the shortfall
  warning when the pool doesn't have enough total years.
- Error paths (missing columns, zero-variance descriptor, duplicate
  station names, ambiguous/missing target specification) are all covered.
"""
import numpy as np
import pandas as pd
import pytest

from floodfreq import pooling as PL


@pytest.fixture
def catalog():
    return pd.DataFrame({
        "station": ["A", "B", "C", "D", "E"],
        "area_km2": [100.0, 105.0, 400.0, 50.0, 110.0],
        "precip_mm": [1000.0, 1010.0, 700.0, 1400.0, 990.0],
        "n_years": [30, 25, 40, 20, 35],
    })


# ---------------------------------------------------------------------------
# read_candidate_catalog
# ---------------------------------------------------------------------------
def test_read_candidate_catalog(tmp_path):
    p = tmp_path / "candidates.csv"
    pd.DataFrame({"station": ["A", "B"], "area_km2": [1.0, 2.0]}).to_csv(p, index=False)
    df = PL.read_candidate_catalog(p)
    assert list(df["station"]) == ["A", "B"]


def test_read_candidate_catalog_missing_station_col_raises(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame({"name": ["A"], "area_km2": [1.0]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="station"):
        PL.read_candidate_catalog(p)


def test_read_candidate_catalog_duplicate_station_raises(tmp_path):
    p = tmp_path / "dupe.csv"
    pd.DataFrame({"station": ["A", "A"], "area_km2": [1.0, 2.0]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        PL.read_candidate_catalog(p)


# ---------------------------------------------------------------------------
# resolve_target
# ---------------------------------------------------------------------------
def test_resolve_target_from_existing_station(catalog):
    t = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    assert t.name == "A"
    assert t.descriptors == {"area_km2": 100.0, "precip_mm": 1000.0}


def test_resolve_target_from_descriptors(catalog):
    t = PL.resolve_target(catalog, ["area_km2", "precip_mm"],
                           target_descriptors={"area_km2": 123.0, "precip_mm": 456.0})
    assert t.name == "(ungauged target)"
    assert t.descriptors == {"area_km2": 123.0, "precip_mm": 456.0}


def test_resolve_target_requires_exactly_one_source(catalog):
    with pytest.raises(ValueError, match="exactly one"):
        PL.resolve_target(catalog, ["area_km2"])
    with pytest.raises(ValueError, match="exactly one"):
        PL.resolve_target(catalog, ["area_km2"], target_station="A",
                           target_descriptors={"area_km2": 1.0})


def test_resolve_target_unknown_station_raises(catalog):
    with pytest.raises(ValueError, match="not found"):
        PL.resolve_target(catalog, ["area_km2"], target_station="ZZZ")


def test_resolve_target_missing_descriptor_value_raises():
    df = pd.DataFrame({"station": ["A"], "area_km2": [1.0], "precip_mm": [np.nan]})
    with pytest.raises(ValueError, match="missing"):
        PL.resolve_target(df, ["area_km2", "precip_mm"], target_station="A")


def test_resolve_target_incomplete_descriptors_dict_raises(catalog):
    with pytest.raises(ValueError, match="missing"):
        PL.resolve_target(catalog, ["area_km2", "precip_mm"],
                           target_descriptors={"area_km2": 1.0})


# ---------------------------------------------------------------------------
# similarity_ranking
# ---------------------------------------------------------------------------
def test_similarity_ranking_near_twin_ranks_first(catalog):
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    # B (105, 1010) is a near-twin of A (100, 1000); C/D are far off
    assert ranking.iloc[0]["station"] == "B"
    assert (ranking["distance"].diff().dropna() >= 0).all()  # monotonically increasing


def test_similarity_ranking_excludes_target_itself(catalog):
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    assert "A" not in set(ranking["station"])
    assert len(ranking) == len(catalog) - 1


def test_similarity_ranking_matches_hand_computed_zscore_distance():
    """Two candidates, one descriptor: distance should be the exact
    z-scored Euclidean distance, computed by hand."""
    df = pd.DataFrame({"station": ["X", "Y"], "val": [0.0, 10.0]})
    target = PL.resolve_target(df, ["val"], target_descriptors={"val": 5.0})
    ranking = PL.similarity_ranking(target, df, ["val"])
    mean, std = df["val"].mean(), df["val"].std(ddof=0)  # 5.0, 5.0
    z_target = (5.0 - mean) / std  # 0.0
    for _, row in ranking.iterrows():
        z_cand = (row["val"] - mean) / std
        expected = abs(z_cand - z_target)
        assert row["distance"] == pytest.approx(expected)


def test_similarity_ranking_weights_change_the_ranking():
    """Two descriptors pulling in different directions: weighting one
    heavily should flip which candidate ranks closer."""
    df = pd.DataFrame({
        "station": ["NEAR_ON_X", "NEAR_ON_Y"],
        "x": [1.0, 100.0],
        "y": [100.0, 1.0],
    })
    # add a third row so std isn't degenerate with only 2 points... actually
    # 2 points give a well-defined nonzero std, so this is fine as-is, but
    # add a neutral third candidate to keep things realistic
    df = pd.concat([df, pd.DataFrame({"station": ["MID"], "x": [50.0], "y": [50.0]})],
                    ignore_index=True)
    target = PL.resolve_target(df, ["x", "y"], target_descriptors={"x": 1.0, "y": 1.0})

    ranking_x_heavy = PL.similarity_ranking(target, df, ["x", "y"], weights={"x": 100.0, "y": 1.0})
    ranking_y_heavy = PL.similarity_ranking(target, df, ["x", "y"], weights={"x": 1.0, "y": 100.0})

    assert ranking_x_heavy.iloc[0]["station"] == "NEAR_ON_X"
    assert ranking_y_heavy.iloc[0]["station"] == "NEAR_ON_Y"


def test_similarity_ranking_missing_descriptor_column_raises(catalog):
    target = PL.resolve_target(catalog, ["area_km2"], target_station="A")
    with pytest.raises(ValueError, match="missing descriptor column"):
        PL.similarity_ranking(target, catalog, ["area_km2", "not_a_column"])


def test_similarity_ranking_zero_variance_descriptor_raises():
    df = pd.DataFrame({"station": ["A", "B", "C"], "const": [5.0, 5.0, 5.0]})
    target = PL.resolve_target(df, ["const"], target_station="A")
    with pytest.raises(ValueError, match="zero variance"):
        PL.similarity_ranking(target, df, ["const"])


def test_similarity_ranking_drops_candidates_missing_descriptor_values(catalog):
    df = catalog.copy()
    df.loc[df["station"] == "C", "precip_mm"] = np.nan
    target = PL.resolve_target(df, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, df, ["area_km2", "precip_mm"])
    assert "C" not in set(ranking["station"])
    assert ranking.attrs["dropped_missing_descriptors"] == ["C"]


def test_similarity_ranking_empty_after_excluding_target_raises():
    df = pd.DataFrame({"station": ["A"], "area_km2": [1.0]})
    target = PL.resolve_target(df, ["area_km2"], target_station="A")
    with pytest.raises(ValueError, match="No candidates left"):
        PL.similarity_ranking(target, df, ["area_km2"])


# ---------------------------------------------------------------------------
# propose_pooling_group
# ---------------------------------------------------------------------------
def test_propose_pooling_group_n_stations(catalog):
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    group = PL.propose_pooling_group(ranking, n_stations=2)
    assert len(group) == 2
    assert list(group["station"]) == list(ranking["station"].iloc[:2])


def test_propose_pooling_group_n_stations_capped_at_available(catalog):
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    group = PL.propose_pooling_group(ranking, n_stations=1000)
    assert len(group) == len(ranking)


def test_propose_pooling_group_n_stations_requires_positive(catalog):
    target = PL.resolve_target(catalog, ["area_km2"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2"])
    with pytest.raises(ValueError, match=">= 1"):
        PL.propose_pooling_group(ranking, n_stations=0)


def test_propose_pooling_group_min_total_years(catalog):
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    # B=25, C=40, D=20, E=35 in some similarity order; find how many are
    # needed by hand from the actual ranking order
    cum = ranking["n_years"].cumsum()
    expected_n = int((cum < 50).sum()) + 1  # first index where cum >= 50, 1-indexed count
    group = PL.propose_pooling_group(ranking, min_total_years=50)
    assert len(group) == expected_n
    assert group["n_years"].sum() >= 50


def test_propose_pooling_group_min_total_years_requires_years_column():
    df = pd.DataFrame({"station": ["A", "B", "C"], "x": [1.0, 2.0, 3.0]})
    target = PL.resolve_target(df, ["x"], target_station="A")
    ranking = PL.similarity_ranking(target, df, ["x"])
    with pytest.raises(ValueError, match="n_years"):
        PL.propose_pooling_group(ranking, min_total_years=100)


def test_propose_pooling_group_min_total_years_shortfall_warns(catalog):
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    total_available = ranking["n_years"].sum()
    with pytest.warns(UserWarning, match="short of the requested"):
        group = PL.propose_pooling_group(ranking, min_total_years=int(total_available) + 1000)
    assert len(group) == len(ranking)  # returns everything available


def test_propose_pooling_group_requires_exactly_one_stopping_rule(catalog):
    target = PL.resolve_target(catalog, ["area_km2"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2"])
    with pytest.raises(ValueError, match="exactly one"):
        PL.propose_pooling_group(ranking)
    with pytest.raises(ValueError, match="exactly one"):
        PL.propose_pooling_group(ranking, n_stations=2, min_total_years=50)


# ---------------------------------------------------------------------------
# End-to-end: ranking + proposal feeds a real regional analysis
# ---------------------------------------------------------------------------
def test_pooling_output_is_usable_by_regional_analysis(catalog):
    """The whole point: the proposed station names should be usable
    directly as keys into a station_data dict for run_regional_analysis()."""
    import lmoments3.distr as ld
    from floodfreq import regional as R

    rng = np.random.default_rng(0)
    target = PL.resolve_target(catalog, ["area_km2", "precip_mm"], target_station="A")
    ranking = PL.similarity_ranking(target, catalog, ["area_km2", "precip_mm"])
    proposed = PL.propose_pooling_group(ranking, n_stations=4)

    station_data = {}
    for name in proposed["station"]:
        n = int(proposed.set_index("station").loc[name, "n_years"])
        x = np.abs(ld.gev(c=-0.1, loc=300, scale=80).rvs(size=n, random_state=rng)) + 1.0
        station_data[name] = x

    result = R.run_regional_analysis("FromPooling", station_data, n_sim=80, seed=1)
    assert set(s.name for s in result.stations) == set(proposed["station"])
