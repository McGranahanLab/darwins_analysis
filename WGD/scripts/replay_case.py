"""Replay steps 04-05 for one tumour from the files Piotr sends for an edge case.

Usage:
  pixi run replay_case [--pgdd templates/pgdd.zip] [--cp-table templates/cp_table.csv] \
      [--ratios templates/ratio_scores.csv] [--tumour LTX0981-Tumour1] [--params params.toml] [--out output/cases]

`--pgdd`: the bootstrapped PGDD run (zip or directory; all tumours). `--cp-table`: ALPACA
cp_table for this tumour (clone x sample). `--ratios`: ALPACA WGD ratio scores. The tumour
defaults to the one in `--ratios`. Clone CN is not given, so the reversion check is skipped.
Writes `<out>/<tumour>/` with events_placed.tsv, tumours_placed.tsv, wgd_events.tsv,
wgd_calls.tsv, modified_tumours.tsv and prints the answer with and without the ploidy prior.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wgd import params as P  # noqa: E402
from wgd.replay import pgdd_tables, read_case, replay  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pgdd", type=Path, default=ROOT / "templates" / "pgdd.zip")
    ap.add_argument("--cp-table", type=Path, default=ROOT / "templates" / "cp_table.csv")
    ap.add_argument("--ratios", type=Path, default=ROOT / "templates" / "ratio_scores.csv")
    ap.add_argument("--tumour")
    ap.add_argument("--params", type=Path, default=ROOT / "params.toml")
    ap.add_argument("--out", type=Path, default=ROOT / "output" / "cases")
    args = ap.parse_args()

    tumour, ratios, cp = read_case(args.cp_table, args.ratios, args.tumour)
    tb = pgdd_tables(args.pgdd, tumour)
    if tb["events"].is_empty() and tb["regions"].is_empty():
        sys.exit(f"{tumour}: not in {args.pgdd}")
    params = P.load(args.params)

    out = args.out / tumour
    out.mkdir(parents=True, exist_ok=True)
    for label, p in (("strict (no prior)", dataclasses.replace(params, prior_enabled=False)), ("final", params)):
        r = replay(tb, ratios, cp, p)
        t = r["tumours"].row(0, named=True)
        print(f"{label:18} GDs {r['gds']:30} status {t['status']:10} quality {t['quality']}  {t['quality_code'] or ''}")
    print(f"ploidy_prior: {t['ploidy_prior']}")
    for k in ("events", "tumours", "wgd_events", "wgd_calls", "modified_tumours"):
        name = {"events": "events_placed", "tumours": "tumours_placed"}.get(k, k)
        r[k].write_csv(out / f"{name}.tsv", separator="\t", null_value="NA")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
