"""Step 05: per-clone WGD table in the previous release format.

Usage:
  pixi run final_table --tx-data-dir <tx_data>/data --out output/<name> [--reject-below 0.3 [--reject-mutation-supported]]

Reads `<out>/04_placed/`, `<out>/01_pgdd/pgdd_tumours.parquet` and tx_data
`alpaca_wgd_ratios`. Writes `<out>/05_table/`:
  modified_tumours.tsv   one row per departure from the plain rules (ploidy prior, forced GDs),
                   with the strict and final solutions, for manual review
  wgd_events.tsv   one row per GD (unique wgd_event_id), with quality, factors and `rejected`
  wgd_calls.tsv    one row per ALPACA clone, all GDs; quality_pass = no GD on the clone below --reject-below; QC-failed tumours included with
                   status qc_failed and NA calls, QC columns on every row (from 04_placed/tumour_qc.tsv)
  wgd_calls_rejected_below_<Q>.tsv   with --reject-below Q: GDs with event_quality < Q dropped,
                   per-clone and per-tumour scores recalculated; mutation-supported GDs are kept
                   unless --reject-mutation-supported
  rejected_events_below_<Q>.tsv      report: one row per rejected GD, tumour score before/after
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
    ap.add_argument("--reject-mutation-supported", action="store_true",
                    help="also allow rejecting mutation-supported GDs (default: always kept)")
    args = ap.parse_args()

    s04 = args.out / "04_placed"
    read = lambda p: pl.read_csv(p, separator="\t", null_values=["NA"], infer_schema_length=None)  # noqa: E731
    events = read(s04 / "events_placed.tsv")
    tumours = read(s04 / "tumours_placed.tsv")
    pgdd_t = pl.read_parquet(args.out / "01_pgdd" / "pgdd_tumours.parquet")
    ratios = pl.read_parquet(args.tx_data_dir / "alpaca_wgd_ratios.parquet")
    qc_path = s04 / "tumour_qc.tsv"
    qc = read(qc_path) if qc_path.exists() else None
    if qc is None:
        print(f"no {qc_path} (step 04 from before 2026-11-16): QC columns left empty")
    ev_table = final_table.build_events(events, tumours, pgdd_t, ratios, args.reject_below,
                                        reject_mutation_supported=args.reject_mutation_supported)
    table = final_table.build_table(events, tumours, pgdd_t, ratios,
                                    event_table=ev_table.with_columns(pl.lit(False).alias("rejected")))
    mods = final_table.modifications(tumours, events, table)
    table = final_table.with_quality_pass(table, ev_table, args.reject_below)
    table = final_table.with_qc(table, qc, ratios)
    out = args.out / "05_table"
    out.mkdir(parents=True, exist_ok=True)
    table.write_csv(out / "wgd_calls.tsv", separator="\t", null_value="NA")
    ev_table.write_csv(out / "wgd_events.tsv", separator="\t", null_value="NA")
    print(f"events: {ev_table.height}")
    mods.write_csv(out / "modified_tumours.tsv", separator="\t", null_value="NA")
    print(f"modified tumours (prior / forced GDs, for manual review): {mods['tumour_id'].n_unique()} "
          f"-> {out / 'modified_tumours.tsv'}")
    if mods.height:
        print(mods.group_by("modification").len().sort("modification"))
    if args.reject_below is not None:
        filt = final_table.build_table(events, tumours, pgdd_t, ratios, reject_below=args.reject_below,
                                       event_table=ev_table)
        filt = final_table.with_quality_pass(filt, ev_table.filter(~pl.col("rejected")), args.reject_below)
        filt = final_table.with_qc(filt, qc, ratios)
        name = f"wgd_calls_rejected_below_{args.reject_below:g}.tsv"
        filt.write_csv(out / name, separator="\t", null_value="NA")
        report = final_table.rejection_report(ev_table, table, filt, args.reject_below)
        rname = f"rejected_events_below_{args.reject_below:g}.tsv"
        report.write_csv(out / rname, separator="\t", null_value="NA")
        kept_mut = ev_table.filter(pl.col("source").is_in(final_table.MUTATION_SOURCES)
                                   & (pl.col("event_quality") < args.reject_below) & ~pl.col("rejected")).height
        print(f"rejected {report.height} GDs (event_quality < {args.reject_below:g}) in "
              f"{report['tumour_id'].n_unique()} tumours; kept {kept_mut} low-quality mutation-supported GDs; "
              f"wrote {out / name} and {out / rname}")
    per_t = table.unique("tumour_id")
    n_qc = per_t.filter(pl.col("status") == final_table.QC_FAILED).height
    n_nt = per_t.filter(pl.col("status") == final_table.NO_TREE).height
    print(f"tumours without calls (NA): {n_qc} qc_failed, {n_nt} no_alpaca_tree")
    print(f"tumours: {per_t.height}, clones: {table.height}")
    print(per_t.group_by("class", "status").len().sort("class", "status"))
    print(f"wrote {out / 'wgd_calls.tsv'}")


if __name__ == "__main__":
    main()
