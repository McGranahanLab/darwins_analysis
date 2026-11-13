"""Step 02: prune PGDD calls against the ALPACA tree.

Usage:
  pixi run python scripts/02_prune_tree.py --tx-data-dir /path/to/tx_data/data --out output/<name>

Reads step 01 output from `<out>/01_pgdd/` and tx_data `clone_proportions`
(ALPACA cp_table: the clones on the tree). Writes `<out>/02_pruned/`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import prune  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tx-data-dir", type=Path, required=True, help="tx_data built parquet dir (data/)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    step01 = args.out / "01_pgdd"
    out = args.out / "02_pruned"
    out.mkdir(parents=True, exist_ok=True)

    events = pl.read_parquet(step01 / "pgdd_events.parquet")
    event_clusters = pl.read_parquet(step01 / "pgdd_event_clusters.parquet")
    tree = prune.tree_clones(pl.read_parquet(args.tx_data_dir / "clone_proportions.parquet"))

    flagged = prune.flag_event_clusters(event_clusters, tree)
    kept, dropped = prune.prune_events(events, flagged)

    flagged.write_parquet(out / "event_clusters_flagged.parquet")
    kept.write_parquet(out / "events_pruned.parquet")
    dropped.write_parquet(out / "events_dropped.parquet")

    summary = {
        "tumours_with_tree": tree["tumour_id"].n_unique(),
        "events_in": events.height,
        "events_kept": kept.height,
        "events_dropped_no_tree": dropped.height,
        "tumours_dropped_no_tree": sorted(dropped["tumour_id"].unique().to_list()),
        "clusters_removed_off_tree": int(flagged.filter(pl.col("tumour_has_tree") & (pl.col("on_tree") == False)).height),  # noqa: E712
        "support": {k: v for k, v in kept.group_by("support").len().sort("support").iter_rows()},
    }
    (out / "run.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))
    with pl.Config(tbl_rows=50, tbl_cols=20, fmt_str_lengths=40, tbl_width_chars=200):
        print(kept.select("GD_event_id", "is_clonal", "clusters", "clusters_pruned", "clusters_removed", "support"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
