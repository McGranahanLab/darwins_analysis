"""Step 07: cohort summary report (quality per WGD pattern with worked examples, cohort overview,
tumours to review by hand).

Usage:
  pixi run report --out output/<name> [--summary <pgdd_cohort_gd_summary.tsv>] [--review-quality 0.3]

Reads `<out>/04_placed/` and `<out>/05_table/`. `--summary` (default
`<out>/03_summary/pgdd_cohort_gd_summary.tsv`, optional) adds the number of samples per tumour.
Writes `<out>/07_report/`: report.html (self-contained), tumour_summary.tsv, examples.tsv,
quality_factors.tsv, review.tsv and the figures as PNG. See wgd/report.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import report  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--summary", type=Path)
    ap.add_argument("--review-quality", type=float, default=0.3, help="review tumours scoring below this")
    ap.add_argument("--prior-review-distance", type=float, default=0.05,
                    help="flag ploidy-prior moves of samples at least this far from 0.5")
    ap.add_argument("--title", default="WGD calls: cohort summary")
    args = ap.parse_args()
    sp = args.summary or args.out / "03_summary" / "pgdd_cohort_gd_summary.tsv"
    summary = pl.read_csv(sp, separator="\t", null_values=["NA", ""], infer_schema_length=None) if sp.exists() else None
    res = report.run(args.out / "07_report", args.out / "04_placed", args.out / "05_table", summary,
                     args.review_quality, args.prior_review_distance, args.title)
    tt, rv = res["tumours"], res["review"]
    print(tt.group_by("pattern").len().sort("pattern"))
    print(f"to review: {rv['tumour_id'].n_unique() if rv.height else 0} tumours")
    print(f"wrote {res['report']}")


if __name__ == "__main__":
    main()
