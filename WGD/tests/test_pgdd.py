"""Parse the tx_data PGDD mock and check it against the scenario truth.

Needs a tx_data checkout with the PGDD mock (PR #5). Set TX_DATA_REPO, or
keep tx_data next to darwins_analysis (../../tx_data from WGD/).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import polars as pl
import pytest

from wgd import pgdd

TX_DATA_REPO = Path(os.environ.get("TX_DATA_REPO", Path(__file__).resolve().parents[3] / "tx_data"))
PGDD_DIR = pgdd.pgdd_dir(TX_DATA_REPO / "nemo_mock")
SCENARIOS = TX_DATA_REPO / "synthetic" / "pgdd" / "pgdd_mock_scenarios.json"

pytestmark = pytest.mark.skipif(not PGDD_DIR.exists(), reason=f"PGDD mock not found at {PGDD_DIR}")


@pytest.fixture(scope="module")
def tables():
    return pgdd.read_pgdd(PGDD_DIR)


def _clusters(ec: pl.DataFrame, event_id: str) -> set[str]:
    return set(ec.filter(pl.col("GD_event_id") == event_id)["cluster_id"].drop_nulls())


def test_ids_canonical(tables):
    for t in ("tumours", "regions", "events"):
        assert tables[t]["tumour_id"].str.contains(r"^LTX\d{4}-Tumour\d+$").all()


def test_all_scenario_tumours_present(tables):
    truth = json.loads(SCENARIOS.read_text())
    assert set(tables["tumours"]["tumour_id"]) == {t["tumour_id"] for t in truth["tumours"]}


def test_comma_columns_are_strings(tables):
    assert tables["events"].schema["clusters"] == pl.Utf8
    assert tables["regions"].schema["gd_clusters"] == pl.Utf8


def test_multi_cluster_event_explodes(tables):
    ec = pgdd.event_clusters(tables["events"])
    assert _clusters(ec, "LTX0005-Tumour1_GD2") == {"3", "6"}
    assert ec.filter(pl.col("GD_event_id") == "LTX0005-Tumour1_GD2")["clusters_raw"].unique().to_list() == ["3,6"]


def test_every_event_kept_once_per_cluster(tables):
    ev = tables["events"]
    ec = pgdd.event_clusters(ev)
    assert set(ec["GD_event_id"]) == set(ev["GD_event_id"])
    # unsupported events (clonal / ploidy-only) keep a single null-cluster row
    unsupported = ev.filter(pl.col("clusters").is_null())["GD_event_id"]
    rows = ec.filter(pl.col("GD_event_id").is_in(unsupported.implode()))
    assert rows["cluster_id"].null_count() == rows.height == unsupported.len()


def test_event_regions_match_counts(tables):
    er = pgdd.event_regions(tables["events"])
    n = er.group_by("GD_event_id").agg(pl.len().alias("n"), pl.col("num_regions_gd_event").first())
    assert (n["n"] == n["num_regions_gd_event"]).all()


def test_non_tree_clusters_are_carried_through(tables):
    # pruning is a later step; parsing must keep calls on non-tree clusters
    ec = pgdd.event_clusters(tables["events"])
    assert _clusters(ec, "LTX0006-Tumour1_GD1") == {"7"}
    rc = pgdd.region_clusters(tables["regions"])
    assert "12" in set(rc.filter(pl.col("tumour_id") == "LTX0007-Tumour1")["cluster_id"])


def test_check_settings(tables):
    diff = pgdd.check_settings(tables)
    assert set(diff) <= {"check_mut_cpn_2_threshold"}  # mock before tx_data#6 used the 1.5 fallback
    bad = {**tables, "settings": tables["settings"].with_columns(pl.lit(100).alias("n_boot"))}
    assert pgdd.check_settings(bad)["n_boot"] == {"expected": 500, "found": [100]}
