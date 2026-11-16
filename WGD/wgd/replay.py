"""Replay steps 04-05 for one tumour (edge cases, `scripts/replay_case.py`, regression tests)."""

from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path

import polars as pl

from wgd import final_table, pgdd, placement
from wgd import params as P


def pgdd_tables(src: Path, tumour: str) -> dict[str, pl.DataFrame]:
    """PGDD tables (zip or directory) restricted to one tumour, plus the cohort summary as `summary`."""
    src = Path(src)
    if src.suffix == ".zip":
        tmp = Path(tempfile.mkdtemp(prefix="pgdd_"))
        zipfile.ZipFile(src).extractall(tmp)
        src = next(p.parent for p in tmp.rglob(pgdd.FILES["events"]))
    t = pgdd.read_pgdd(src)
    out = {k: (v.filter(pl.col("tumour_id") == tumour) if "tumour_id" in v.columns else v) for k, v in t.items()}
    summary = pl.read_csv(src / "pgdd_cohort_gd_summary.tsv", separator="\t", null_values=["NA"],
                          infer_schema_length=None).filter(pl.col("tumour_id") == tumour)
    return out | {"summary": summary}


def read_case(cp_table: Path, ratio_scores: Path, tumour: str | None = None) -> tuple[str, pl.DataFrame, pl.DataFrame]:
    """(tumour, ratios, clone proportions) from Piotr's cp_table.csv (clone x sample) and ratio_scores.csv."""
    ratios = pl.read_csv(ratio_scores)
    tumour = tumour or ratios["tumour_id"].unique().item()
    ratios = ratios.filter(pl.col("tumour_id") == tumour).select("tumour_id", "clone", "parent", "ratio_A", "ratio_B")
    cp = (pl.read_csv(cp_table).unpivot(index="clone", variable_name="sample_id", value_name="proportion")
            .with_columns(pl.lit(tumour).alias("tumour_id")))
    return tumour, ratios, cp


def replay(tb: dict[str, pl.DataFrame], ratios: pl.DataFrame, cp: pl.DataFrame,
           params: P.Params | None = None) -> dict:
    """Steps 04 (fill) and 05 on one tumour. Clone CN is not given, so no reversion check."""
    res = placement.place_events(tb["events"], tb["tumours"], tb["mut_counts"], ratios, cp, summary=tb["summary"],
                                 ploidy_rule="fill", regions=tb["regions"], params=params)
    ev = final_table.build_events(res.events, res.tumours, tb["tumours"], ratios)
    calls = final_table.build_table(res.events, res.tumours, tb["tumours"], ratios, event_table=ev)
    gds = ",".join(f"{c}={n}" for c, n in sorted(zip(calls["clone"], calls["number_of_gds_at_clone"])) if n) or "none"
    return {"events": res.events, "tumours": res.tumours, "wgd_events": ev, "wgd_calls": calls, "gds": gds,
            "modified_tumours": final_table.modifications(res.tumours, res.events, calls)}
