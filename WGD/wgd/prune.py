"""Prune PGDD calls against the final tree.

The tree is the set of clones ALPACA has for each tumour: the clones in
its `cp_table.csv` (tx_data `clone_proportions`) or `tree_paths.json`. PGDD runs on
pre-tree mutation clusters, so:

* a tumour with no ALPACA clones has no tree -> all its events are dropped;
* a supporting cluster that is not on the tree is removed from its event;
* an event whose supporting clusters are all removed is KEPT as
  ploidy-only: the regions still carry the doubling, but the clone it sits
  on is no longer known (Piotr, 2026-11-06).

PGDD cluster "3" is tree clone "clone3".
"""

from __future__ import annotations

import polars as pl

SUPPORT_LEVELS = (
    "clonal",                   # clonal event (no subclonal support needed)
    "mut_supported",            # subclonal, >= 1 supporting cluster on the tree
    "ploidy_only",              # subclonal, PGDD found no supporting cluster
    "ploidy_only_pruned",       # subclonal, all PGDD supporting clusters were off-tree
)


def tree_clones(clone_proportions: pl.DataFrame) -> pl.DataFrame:
    """Distinct (tumour_id, clone) on the ALPACA tree."""
    return clone_proportions.select("tumour_id", "clone").unique().sort("tumour_id", "clone")


def flag_event_clusters(event_clusters: pl.DataFrame, tree: pl.DataFrame) -> pl.DataFrame:
    """Add `clone`, `tumour_has_tree` and `on_tree` to the event x cluster table (step 01)."""
    tumours = tree.select("tumour_id").unique().with_columns(pl.lit(True).alias("tumour_has_tree"))
    on_tree = tree.with_columns(pl.lit(True).alias("on_tree"))
    return (
        event_clusters
        .with_columns(pl.when(pl.col("cluster_id").is_not_null())
                      .then(pl.lit("clone") + pl.col("cluster_id")).alias("clone"))
        .join(tumours, on="tumour_id", how="left")
        .join(on_tree, on=["tumour_id", "clone"], how="left")
        .with_columns(
            pl.col("tumour_has_tree").fill_null(False),
            # null for events with no supporting cluster; False = cluster off the tree
            pl.when(pl.col("clone").is_not_null()).then(pl.col("on_tree").fill_null(False)).alias("on_tree"),
        )
    )


def prune_events(events: pl.DataFrame, flagged: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return (pruned events, dropped events).

    Pruned events: one row per kept event, with the original PGDD columns plus
    `clusters_pruned` (comma-joined on-tree clusters, null if none),
    `clusters_removed`, `n_clusters_raw`, `n_clusters_kept` and `support`.
    Dropped events: events of tumours without a tree, with `drop_reason`.
    """
    per_event = flagged.group_by("GD_event_id").agg(
        pl.col("tumour_has_tree").first(),
        pl.col("cluster_id").filter(pl.col("on_tree") == True).sort().alias("_kept"),  # noqa: E712
        pl.col("cluster_id").filter(pl.col("on_tree") == False).sort().alias("_removed"),  # noqa: E712
        pl.col("cluster_id").drop_nulls().len().alias("n_clusters_raw"),
    )
    ev = events.join(per_event, on="GD_event_id", how="left")

    dropped = (ev.filter(~pl.col("tumour_has_tree"))
                 .drop("_kept", "_removed", "tumour_has_tree")
                 .with_columns(pl.lit("no_tree").alias("drop_reason")))

    kept = ev.filter(pl.col("tumour_has_tree")).with_columns(
        pl.col("_kept").list.len().alias("n_clusters_kept"),
        pl.when(pl.col("_kept").list.len() > 0).then(pl.col("_kept").list.join(",")).alias("clusters_pruned"),
        pl.when(pl.col("_removed").list.len() > 0).then(pl.col("_removed").list.join(",")).alias("clusters_removed"),
    ).with_columns(
        pl.when(pl.col("is_clonal")).then(pl.lit("clonal"))
          .when(pl.col("n_clusters_kept") > 0).then(pl.lit("mut_supported"))
          .when(pl.col("n_clusters_raw") > 0).then(pl.lit("ploidy_only_pruned"))
          .otherwise(pl.lit("ploidy_only"))
          .cast(pl.Enum(SUPPORT_LEVELS)).alias("support")
    ).drop("_kept", "_removed", "tumour_has_tree")

    return kept.sort("tumour_id", "GD_event"), dropped.sort("tumour_id", "GD_event")
