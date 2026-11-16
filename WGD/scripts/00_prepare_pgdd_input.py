"""Step 00: build the PGDD input table from CONIPHER outputs and MSAI CN segments.

Python port of 00_prepare_pgdd_input.R. Run before PGDD (detect_par_gd).

Usage:
  pixi run prepare_pgdd_input --conipher-dir <dir> --msai-dir <dir> --out output/<name>
      [--keep-nontree-clusters] [--tumours LTX0001_1 LTX0002_1]

`--conipher-dir` holds one directory per tumour (<tumour>.SCoutput.CLEAN.tsv,
treeTable.tsv, allTrees.txt); `--msai-dir` holds
<tumour>/<tumour>_copynumber_MSAI_adjusted_segments_tree_input.tsv.
A pre-flight pass checks every tumour first. Any missing/empty file, missing
column, unreadable tree, or SAMPLE ids differing between SCoutput, treeTable and
the CN segments stops the run (exit 1) with all problems listed.
`--preflight-only` runs just that check.

Writes to `<out>/00_pgdd_input/`:
  pgdd_input_all_tumours.tsv   the file PGDD reads (tumours with a tree)
  single_clone_samples.tsv     per-sample num_gds of single-clone tumours (not run through PGDD)
  single_clone_events.tsv      per single-clone tumour: n_clonal_gds on its only clone, quality, quality_code
  run.json                     counts
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import polars as pl  # noqa: E402

from wgd import pgdd_input  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conipher-dir", type=Path, required=True)
    ap.add_argument("--msai-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--keep-nontree-clusters", action="store_true",
                    help="keep mutations whose tree cluster is not an edge of the default tree")
    ap.add_argument("--tumours", nargs="+", help="only these tumour directories")
    ap.add_argument("--preflight-only", action="store_true", help="only check inputs, write nothing")
    args = ap.parse_args()

    try:
        if args.preflight_only:
            dirs = sorted(p.name for p in args.conipher_dir.iterdir() if p.is_dir())
            if args.tumours is not None:
                missing = sorted(set(args.tumours) - set(dirs))
                if missing:
                    raise pgdd_input.InputError(f"selected tumours not found: {', '.join(missing)}")
                dirs = [d for d in dirs if d in args.tumours]
            pgdd_input.preflight(args.conipher_dir, args.msai_dir, dirs)
            print(f"pre-flight OK: {len(dirs)} tumours")
            return
        res = pgdd_input.prepare_all(args.conipher_dir, args.msai_dir, tumours=args.tumours,
                                     remove_nontree_clusters=not args.keep_nontree_clusters)
    except pgdd_input.InputError as e:
        sys.exit(f"ERROR: {e}")
    out = args.out / "00_pgdd_input"
    out.mkdir(parents=True, exist_ok=True)
    pgdd_input.write_tsv(res.table, out / "pgdd_input_all_tumours.tsv")
    pgdd_input.write_tsv(res.single_clone, out / "single_clone_samples.tsv")
    events = pgdd_input.single_clone_events(res.single_clone)
    pgdd_input.write_tsv(events, out / "single_clone_events.tsv")

    t = res.table
    run = {
        "conipher_dir": str(args.conipher_dir),
        "msai_dir": str(args.msai_dir),
        "remove_nontree_clusters": not args.keep_nontree_clusters,
        "tumours": t["tumour_id"].n_unique(),
        "samples": t["sample_id"].n_unique(),
        "mutations": t.height,
        "single_clone_tumours": res.single_clone["tumour_id"].n_unique(),
        "single_clone_low_quality": {t: c for t, c in events.filter(pl.col("quality") == "Low")
                                     .select("tumour_id", "quality_code").iter_rows()},
    }
    (out / "run.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps(run, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
