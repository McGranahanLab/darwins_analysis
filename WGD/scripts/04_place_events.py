"""Step 04: place PGDD GD events on ALPACA clones (port of the previous notebook's rules A-I).

Usage:
  pixi run place_events --tx-data-dir <tx_data>/data --out output/<name> \
      [--summary <pgdd_cohort_gd_summary.tsv>] [--qc-table <manual_qc.tsv> --release-dir <_RELEASE/release_tx842>] \
      [--alpaca-score ratio_A|max_AB]

Reads step 01 output (`<out>/01_pgdd/`) and tx_data `alpaca_wgd_ratios`, `clone_proportions`
and `alpaca` (for the reversion flag). `--summary` defaults to
`<out>/03_summary/pgdd_cohort_gd_summary.tsv` (used for the H flag). Without
`--qc-table`/`--release-dir` no tumour QC is applied.
`--alpaca-score` defaults to `ratio_A`, as in the notebook; `max_AB` is the 2026-11-09 decision.
Writes `<out>/04_placed/`: events_placed.tsv, tumours_placed.tsv, run.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import placement  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tx-data-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--summary", type=Path)
    ap.add_argument("--qc-table", type=Path)
    ap.add_argument("--release-dir", type=Path)
    ap.add_argument("--alpaca-score", choices=("ratio_A", "max_AB"), default="ratio_A")
    args = ap.parse_args()
    if (args.qc_table is None) != (args.release_dir is None):
        ap.error("give both --qc-table and --release-dir, or neither")

    s01 = args.out / "01_pgdd"
    tx = args.tx_data_dir
    summary_path = args.summary or args.out / "03_summary" / "pgdd_cohort_gd_summary.tsv"
    summary = None
    if summary_path.exists():
        summary = pl.read_csv(summary_path, separator="\t", null_values=["NA", ""], infer_schema_length=None)
    else:
        print(f"no summary at {summary_path}: H flag skipped")

    qc = None
    if args.qc_table:
        tx_qc, qc = placement.qc_pass(args.qc_table, args.release_dir)
        print(f"QC: {len(tx_qc)} tumours FINALIZED, {len(qc)} also pass truncal CCF > 90%")

    res = placement.place_events(
        events=pl.read_parquet(s01 / "pgdd_events.parquet"),
        tumours=pl.read_parquet(s01 / "pgdd_tumours.parquet"),
        mut_counts=pl.read_parquet(s01 / "pgdd_mut_counts.parquet"),
        ratios=pl.read_parquet(tx / "alpaca_wgd_ratios.parquet"),
        clone_proportions=pl.read_parquet(tx / "clone_proportions.parquet"),
        summary=summary,
        alpaca_cn=pl.read_parquet(tx / "alpaca.parquet") if (tx / "alpaca.parquet").exists() else None,
        qc_tumours=qc, score=args.alpaca_score)

    out = args.out / "04_placed"
    out.mkdir(parents=True, exist_ok=True)
    res.events.write_csv(out / "events_placed.tsv", separator="\t", null_value="NA")
    res.tumours.write_csv(out / "tumours_placed.tsv", separator="\t", null_value="NA")
    t = res.tumours
    run = {
        "alpaca_score": args.alpaca_score,
        "qc_applied": qc is not None,
        "tumours": t.height,
        "triage": dict(t.group_by("triage").len().sort("triage").iter_rows()),
        "status": dict(t.group_by("status").len().sort("status").iter_rows()),
        "tumours_without_alpaca_tree": t.filter(~pl.col("has_alpaca_tree"))["tumour_id"].to_list(),
        "events": res.events.height,
        "events_unplaced": int(res.events["placed_clone"].null_count()) if res.events.height else 0,
        "placed_as_kind": (dict(res.events.with_columns(pl.col("placed_as").fill_null("NA")
                                .str.replace(r"_.*", "").str.replace(r"^\d+$", "cluster")
                                .str.replace(r"^\d+(,\d+)+$", "multi_unresolved").alias("k"))
                                .group_by("k").len().sort("k").iter_rows()) if res.events.height else {}),
        "H_mismatch_tumours": t.filter("H_mismatch")["tumour_id"].to_list(),
        "log": res.log,
    }
    (out / "run.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps({k: v for k, v in run.items() if k != "log"}, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
