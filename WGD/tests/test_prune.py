"""Step 02 pruning on the tx_data PGDD mock + built clone_proportions."""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from wgd import pgdd, prune

TX_DATA_REPO = Path(os.environ.get("TX_DATA_REPO", Path(__file__).resolve().parents[3] / "tx_data"))
PGDD_DIR = pgdd.pgdd_dir(TX_DATA_REPO / "nemo_mock")
CP = TX_DATA_REPO / "data" / "clone_proportions.parquet"

pytestmark = pytest.mark.skipif(not (PGDD_DIR.exists() and CP.exists()),
                                reason="needs tx_data PGDD mock and built clone_proportions")


@pytest.fixture(scope="module")
def result():
    events = pgdd.read_pgdd(PGDD_DIR)["events"]
    tree = prune.tree_clones(pl.read_parquet(CP))
    flagged = prune.flag_event_clusters(pgdd.event_clusters(events), tree)
    kept, dropped = prune.prune_events(events, flagged)
    return events, flagged, kept, dropped


def _row(df: pl.DataFrame, event_id: str) -> dict:
    rows = df.filter(pl.col("GD_event_id") == event_id).to_dicts()
    assert len(rows) == 1, event_id
    return rows[0]


def test_every_event_kept_or_dropped(result):
    events, _, kept, dropped = result
    assert kept.height + dropped.height == events.height
    assert set(kept["GD_event_id"]).isdisjoint(dropped["GD_event_id"])


def test_tumours_without_tree_dropped(result):
    _, _, kept, dropped = result
    assert set(dropped["tumour_id"]) == {"LTX0009-Tumour1", "LTX0010-Tumour1"}
    assert (dropped["drop_reason"] == "no_tree").all()
    assert not kept["tumour_id"].is_in(["LTX0009-Tumour1", "LTX0010-Tumour1"]).any()


def test_off_tree_cluster_removed_from_multi_cluster_event(result):
    r = _row(result[2], "LTX0005-Tumour1_GD2")
    assert (r["clusters"], r["clusters_pruned"], r["clusters_removed"]) == ("3,6", "3", "6")
    assert r["support"] == "mut_supported"


def test_all_clusters_off_tree_becomes_ploidy_only(result):
    r = _row(result[2], "LTX0006-Tumour1_GD1")
    assert r["clusters_pruned"] is None and r["clusters_removed"] == "7"
    assert r["support"] == "ploidy_only_pruned"
    # PGDD's own flag is left untouched; `support` carries the pruned state
    assert r["is_subclonal_mutation_supported"] is True


def test_on_tree_and_clonal_untouched(result):
    kept = result[2]
    r = _row(kept, "LTX0004-Tumour1_GD1")
    assert (r["clusters_pruned"], r["support"]) == ("3", "mut_supported")
    assert _row(kept, "LTX0002-Tumour1_GD1")["support"] == "clonal"


def test_support_consistent_with_counts(result):
    kept = result[2]
    sub = kept.filter(~pl.col("is_clonal"))
    assert ((sub["support"] == "mut_supported") == (sub["n_clusters_kept"] > 0)).all()
    assert (kept["n_clusters_kept"] <= kept["n_clusters_raw"]).all()
