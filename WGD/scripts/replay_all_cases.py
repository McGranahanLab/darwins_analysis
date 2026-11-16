"""Replay every regression case (tests/cases/<tumour>/) and write them as one small cohort.

Usage:
  pixi run replay_all_cases [--out output/cases_cohort] [--params params.toml]

Writes `<out>/04_placed/` (events_placed.tsv, tumours_placed.tsv) and `<out>/05_table/`
(wgd_calls.tsv, wgd_events.tsv, modified_tumours.tsv), the layout step 07 (report) reads.
Useful to try the report before a full cohort run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wgd import params as P  # noqa: E402
from wgd.replay import pgdd_tables, read_case, replay  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=Path, default=ROOT / "tests" / "cases")
    ap.add_argument("--params", type=Path, default=ROOT / "params.toml")
    ap.add_argument("--out", type=Path, default=ROOT / "output" / "cases_cohort")
    args = ap.parse_args()
    params = P.load(args.params)
    parts: dict[str, list[pl.DataFrame]] = {}
    for d in sorted(p for p in args.cases.iterdir() if p.is_dir()):
        t, ratios, cp = read_case(d / "cp_table.csv", d / "ratio_scores.csv", d.name)
        r = replay(pgdd_tables(d / "pgdd", t), ratios, cp, params)
        for k in ("events", "tumours", "wgd_events", "wgd_calls", "modified_tumours"):
            parts.setdefault(k, []).append(r[k])
        print(f"{t}: {r['gds']}")
    for sub, keys in (("04_placed", {"events": "events_placed", "tumours": "tumours_placed"}),
                      ("05_table", {"wgd_events": "wgd_events", "wgd_calls": "wgd_calls",
                                    "modified_tumours": "modified_tumours"})):
        (args.out / sub).mkdir(parents=True, exist_ok=True)
        for k, name in keys.items():
            pl.concat(parts[k], how="diagonal_relaxed").write_csv(args.out / sub / f"{name}.tsv",
                                                                   separator="\t", null_value="NA")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
