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
import tempfile
import zipfile
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from wgd import final_table, pgdd, placement  # noqa: E402
from wgd import params as P  # noqa: E402


def pgdd_tables(src: Path, tumour: str) -> dict[str, pl.DataFrame]:
    if src.suffix == ".zip":
        tmp = Path(tempfile.mkdtemp(prefix="pgdd_"))
        zipfile.ZipFile(src).extractall(tmp)
        src = next(p.parent for p in tmp.rglob(pgdd.FILES["events"]))
    t = pgdd.read_pgdd(src)
    out = {k: (v.filter(pl.col("tumour_id") == tumour) if "tumour_id" in v.columns else v) for k, v in t.items()}
    summary = pl.read_csv(src / "pgdd_cohort_gd_summary.tsv", separator="\t", null_values=["NA"],
                          infer_schema_length=None).filter(pl.col("tumour_id") == tumour)
    return out | {"summary": summary}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pgdd", type=Path, default=ROOT / "templates" / "pgdd.zip")
    ap.add_argument("--cp-table", type=Path, default=ROOT / "templates" / "cp_table.csv")
    ap.add_argument("--ratios", type=Path, default=ROOT / "templates" / "ratio_scores.csv")
    ap.add_argument("--tumour")
    ap.add_argument("--params", type=Path, default=ROOT / "params.toml")
    ap.add_argument("--out", type=Path, default=ROOT / "output" / "cases")
    args = ap.parse_args()

    ratios = pl.read_csv(args.ratios)
    tumour = args.tumour or ratios["tumour_id"].unique().item()
    ratios = ratios.filter(pl.col("tumour_id") == tumour).select("tumour_id", "clone", "parent", "ratio_A", "ratio_B")
    cp = (pl.read_csv(args.cp_table).unpivot(index="clone", variable_name="sample_id", value_name="proportion")
            .with_columns(pl.lit(tumour).alias("tumour_id")))
    tb = pgdd_tables(args.pgdd, tumour)
    if tb["events"].is_empty() and tb["regions"].is_empty():
        sys.exit(f"{tumour}: not in {args.pgdd}")
    params = P.load(args.params)

    out = args.out / tumour
    out.mkdir(parents=True, exist_ok=True)
    for label, p in (("strict (no prior)", dataclasses.replace(params, prior_enabled=False)), ("final", params)):
        res = placement.place_events(tb["events"], tb["tumours"], tb["mut_counts"], ratios, cp,
                                     summary=tb["summary"], ploidy_rule="fill", regions=tb["regions"], params=p)
        ev = final_table.build_events(res.events, res.tumours, tb["tumours"], ratios)
        calls = final_table.build_table(res.events, res.tumours, tb["tumours"], ratios, event_table=ev)
        r = res.tumours.row(0, named=True)
        gds = ",".join(f"{c}={n}" for c, n in zip(calls["clone"], calls["number_of_gds_at_clone"]) if n) or "none"
        print(f"{label:18} GDs {gds:30} status {r['status']:10} quality {r['quality']}  {r['quality_code'] or ''}")
    print(f"ploidy_prior: {r['ploidy_prior']}")
    res.events.write_csv(out / "events_placed.tsv", separator="\t", null_value="NA")
    res.tumours.write_csv(out / "tumours_placed.tsv", separator="\t", null_value="NA")
    ev.write_csv(out / "wgd_events.tsv", separator="\t", null_value="NA")
    calls.write_csv(out / "wgd_calls.tsv", separator="\t", null_value="NA")
    final_table.modifications(res.tumours, res.events, calls).write_csv(out / "modified_tumours.tsv",
                                                                         separator="\t", null_value="NA")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
