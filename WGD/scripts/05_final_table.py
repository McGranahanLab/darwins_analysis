"""Step 05: per-clone WGD table in the previous release format.

Usage:
  pixi run final_table --tx-data-dir <tx_data>/data --out output/<name>

Reads `<out>/04_placed/`, `<out>/01_pgdd/pgdd_tumours.parquet` and tx_data
`alpaca_wgd_ratios`. Writes `<out>/05_table/wgd_calls.tsv`.
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
    args = ap.parse_args()

    s04 = args.out / "04_placed"
    read = lambda p: pl.read_csv(p, separator="\t", null_values=["NA"], infer_schema_length=None)  # noqa: E731
    table = final_table.build_table(
        read(s04 / "events_placed.tsv"), read(s04 / "tumours_placed.tsv"),
        pl.read_parquet(args.out / "01_pgdd" / "pgdd_tumours.parquet"),
        pl.read_parquet(args.tx_data_dir / "alpaca_wgd_ratios.parquet"))
    out = args.out / "05_table"
    out.mkdir(parents=True, exist_ok=True)
    table.write_csv(out / "wgd_calls.tsv", separator="\t", null_value="NA")
    per_t = table.unique("tumour_id")
    print(f"tumours: {per_t.height}, clones: {table.height}")
    print(per_t.group_by("class", "status").len().sort("class", "status"))
    print(f"wrote {out / 'wgd_calls.tsv'}")


if __name__ == "__main__":
    main()
