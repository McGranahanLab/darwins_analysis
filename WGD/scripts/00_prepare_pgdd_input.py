"""Step 00: build the PGDD input table from CONIPHER outputs and MSAI CN segments.

Python port of 00_prepare_pgdd_input.R. Run before PGDD (detect_par_gd).

Usage:
  pixi run prepare_pgdd_input --conipher-dir <dir> --msai-dir <dir> --out output/<name>
      [--keep-nontree-clusters] [--tumours LTX0001_1 LTX0002_1]

`--conipher-dir` holds one directory per tumour (<tumour>.SCoutput.CLEAN.tsv,
treeTable.tsv, allTrees.txt); `--msai-dir` holds
<tumour>/<tumour>_copynumber_MSAI_adjusted_segments_tree_input.tsv.
Writes `<out>/00_pgdd_input/pgdd_input_all_tumours.tsv` (the file PGDD reads)
and `run.json` (counts and skipped tumours with reasons).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import pgdd_input  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conipher-dir", type=Path, required=True)
    ap.add_argument("--msai-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--keep-nontree-clusters", action="store_true",
                    help="keep mutations whose tree cluster is not an edge of the default tree")
    ap.add_argument("--tumours", nargs="+", help="only these tumour directories")
    args = ap.parse_args()

    res = pgdd_input.prepare_all(args.conipher_dir, args.msai_dir, tumours=args.tumours,
                                 remove_nontree_clusters=not args.keep_nontree_clusters)
    out = args.out / "00_pgdd_input"
    out.mkdir(parents=True, exist_ok=True)
    pgdd_input.write_tsv(res.table, out / "pgdd_input_all_tumours.tsv")

    t = res.table
    run = {
        "conipher_dir": str(args.conipher_dir),
        "msai_dir": str(args.msai_dir),
        "remove_nontree_clusters": not args.keep_nontree_clusters,
        "tumours": t["tumour_id"].n_unique(),
        "samples": t["sample_id"].n_unique(),
        "mutations": t.height,
        "skipped": res.skipped,
    }
    (out / "run.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps({k: v for k, v in run.items() if k != "skipped"} | {"skipped": len(res.skipped)}, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
