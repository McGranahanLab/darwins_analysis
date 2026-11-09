"""Step 03: per-sample PGDD cohort summary (pgdd_cohort_gd_summary.tsv).

Python port of 03_generate_summary.R.

Usage:
  pixi run generate_summary --pgdd-input output/<name>/00_pgdd_input/pgdd_input_all_tumours.tsv \\
      --pgdd-regions <pgdd_dir>/pgdd_GDs_per_region.tsv --out output/<name>

Writes `<out>/03_summary/pgdd_cohort_gd_summary.tsv`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import summary  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pgdd-input", type=Path, required=True, help="pgdd_input_all_tumours.tsv (step 00)")
    ap.add_argument("--pgdd-regions", type=Path, required=True, help="pgdd_GDs_per_region.tsv")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    table = summary.cohort_summary(summary.read_tsv(args.pgdd_input), summary.read_tsv(args.pgdd_regions))
    out = args.out / "03_summary"
    out.mkdir(parents=True, exist_ok=True)
    table.write_csv(out / "pgdd_cohort_gd_summary.tsv", separator="\t", null_value="NA")
    print(f"samples: {table.height}")
    print(f"wrote {out / 'pgdd_cohort_gd_summary.tsv'}")


if __name__ == "__main__":
    main()
