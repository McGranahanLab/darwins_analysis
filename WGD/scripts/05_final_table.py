"""Step 05: per-clone WGD table in the previous release format.

Usage:
  pixi run final_table --tx-data-dir <tx_data>/data --out output/<name> [--reject-below 0.3]

Reads `<out>/04_placed/`, `<out>/01_pgdd/pgdd_tumours.parquet` and tx_data
`alpaca_wgd_ratios`. Writes `<out>/05_table/`:
  wgd_events.tsv   one row per GD (unique wgd_event_id), with quality, factors and `rejected`
  wgd_calls.tsv    one row per ALPACA clone, all GDs
  wgd_calls_rejected_below_<Q>.tsv   with --reject-below Q: GDs with event_quality < Q dropped,
                   per-clone and per-tumour scores recalculated
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import final_table  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tx-data-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--reject-below", type=float, help="drop GDs with event_quality below this and recompute")
    args = ap.parse_args()

    s04 = args.out / "04_placed"
    read = lambda p: pl.read_csv(p, separator="\t", null_values=["NA"], infer_schema_length=None)  # noqa: E731
    events = read(s04 / "events_placed.tsv")
    tumours = read(s04 / "tumours_placed.tsv")
    pgdd_t = pl.read_parquet(args.out / "01_pgdd" / "pgdd_tumours.parquet")
    ratios = pl.read_parquet(args.tx_data_dir / "alpaca_wgd_ratios.parquet")
    ev_table = final_table.build_events(events, tumours, pgdd_t, ratios, args.reject_below)
    table = final_table.build_table(events, tumours, pgdd_t, ratios,
                                    event_table=ev_table.with_columns(pl.lit(False).alias("rejected")))
    out = args.out / "05_table"
    out.mkdir(parents=True, exist_ok=True)
    table.write_csv(out / "wgd_calls.tsv", separator="\t", null_value="NA")
    ev_table.write_csv(out / "wgd_events.tsv", separator="\t", null_value="NA")
    print(f"events: {ev_table.height}")
    if args.reject_below is not None:
        filt = final_table.build_table(events, tumours, pgdd_t, ratios, reject_below=args.reject_below,
                                       event_table=ev_table)
        name = f"wgd_calls_rejected_below_{args.reject_below:g}.tsv"
        filt.write_csv(out / name, separator="\t", null_value="NA")
        rej = ev_table.filter("rejected")
        print(f"rejected {rej.height} GDs (event_quality < {args.reject_below:g}) in "
              f"{rej['tumour_id'].n_unique()} tumours; wrote {out / name}")
    per_t = table.unique("tumour_id")
    print(f"tumours: {per_t.height}, clones: {table.height}")
    print(per_t.group_by("class", "status").len().sort("class", "status"))
    print(f"wrote {out / 'wgd_calls.tsv'}")


if __name__ == "__main__":
    main()
