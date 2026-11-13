"""Correct ALPACA clone copy number with the agreed WGD calls (prototype).

Usage:
  pixi run correct --tx-data-dir <tx_data>/data --out output/<name> \
      [--wgd-calls <WGD output>/05_table/wgd_calls.tsv] [--min-quality 0.3]

Reads tx_data `alpaca` (clone CN per segment), `clone_proportions` (which clones have cells)
and `alpaca_wgd_ratios` (tree). WGD calls: `--wgd-calls` (the WGD pipeline's step 05
wgd_calls.tsv), else tx_data `wgd_calls` (placeholder until the WGD table is released).
Writes `<out>/`: alpaca_corrected.parquet (tx_data alpaca schema + pred_CN_*_original +
corrected), wgd_ratios_corrected.tsv (ALPACA's 12-column format), ratio_comparison.tsv,
corrections.tsv (one row per corrected clone), skipped.tsv, flags.tsv (ALPACA vs calls
disagreements left after correction: flagged, not fixed), run.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from apc import correct  # noqa: E402


def read_any(p: Path) -> pl.DataFrame:
    if p.suffix == ".parquet":
        return pl.read_parquet(p)
    return pl.read_csv(p, separator="\t" if p.suffix in (".tsv", ".txt") else ",", null_values=["NA"],
                       infer_schema_length=None)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tx-data-dir", type=Path, required=True)
    ap.add_argument("--wgd-calls", type=Path)
    ap.add_argument("--min-quality", type=float, help="correct only tumours whose WGD call quality is at least this")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    tx = args.tx_data_dir
    calls_path = args.wgd_calls or tx / "wgd_calls.parquet"
    res = correct.run(pl.read_parquet(tx / "alpaca.parquet"), pl.read_parquet(tx / "clone_proportions.parquet"),
                      read_any(calls_path), pl.read_parquet(tx / "alpaca_wgd_ratios.parquet"), args.min_quality)
    args.out.mkdir(parents=True, exist_ok=True)
    res["alpaca_corrected"].write_parquet(args.out / "alpaca_corrected.parquet")
    for k in ("wgd_ratios_corrected", "ratio_comparison", "corrections", "skipped", "flags"):
        res[k].write_csv(args.out / f"{k}.tsv", separator="\t", null_value="NA")
    cn = res["alpaca_corrected"]
    run = {"wgd_calls": str(calls_path), "min_quality": args.min_quality,
           "tumours": cn["tumour_id"].n_unique(), "tumours_corrected": res["corrections"]["tumour_id"].n_unique()
           if res["corrections"].height else 0,
           "clones_corrected": res["corrections"].height, "rows_changed": int(cn["corrected"].sum()),
           "skipped": res["skipped"].height,
           "flags": dict(res["flags"].group_by("flag").len().iter_rows()) if res["flags"].height else {}}
    (args.out / "run.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps(run, indent=2))
    if res["corrections"].height:
        print(res["corrections"].select("tumour_id", "clone", "children_with_gd", "bp_fraction_changed",
                                        "ratio_A_before", "ratio_A_after", "children_ratio_before",
                                        "children_ratio_after"))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
