"""Step 06: compare the new table (step 05) with the previous release table.

Usage:
  pixi run compare_tables --old <previous table .tsv/.csv> --out output/<name>

Writes `<out>/06_compare/`: summary.json (counts only), per_tumour.tsv
(which fields agree per tumour), clone_diff.tsv (clones that differ in
anything but status).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import compare  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    new = pl.read_csv(args.out / "05_table" / "wgd_calls.tsv", separator="\t", null_values=["NA"],
                      infer_schema_length=None)
    res = compare.compare(compare.read_table(args.old), new)
    out = args.out / "06_compare"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(res["summary"], indent=2) + "\n")
    res["per_tumour"].write_csv(out / "per_tumour.tsv", separator="\t", null_value="NA")
    res["clone_diff"].write_csv(out / "clone_diff.tsv", separator="\t", null_value="NA")
    print(json.dumps(res["summary"], indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
