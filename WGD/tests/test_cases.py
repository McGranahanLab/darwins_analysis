"""Regression: every real tumour discussed with Piotr, replayed from its own files.

tests/cases/<tumour>/: cp_table.csv and ratio_scores.csv as Piotr sent them (git history of
templates/), pgdd/ = that tumour's rows of templates/pgdd.zip, snapshot.json = the full output
(GDs, status, scores, codes, per-event scores) when it was last reviewed.

* test_agreed_answer fails when a rule change alters an answer agreed with Piotr
  (tests/cases/agreed.toml): a clash to discuss first.
* test_snapshot fails when anything else in the output changes (a score, a code). Review the
  diff, tell Piotr, then refresh with `UPDATE_CASES=1 pixi run test`.
"""

from __future__ import annotations

import json
import os
import tomllib
from pathlib import Path

import pytest

from wgd import params as P
from wgd.replay import pgdd_tables, read_case, replay

CASES = Path(__file__).parent / "cases"
AGREED = tomllib.loads((CASES / "agreed.toml").read_text())
PARAMS = P.load(Path(__file__).resolve().parents[1] / "params.toml")


def run_case(tumour: str) -> dict:
    d = CASES / tumour
    t, ratios, cp = read_case(d / "cp_table.csv", d / "ratio_scores.csv", tumour)
    r = replay(pgdd_tables(d / "pgdd", t), ratios, cp, PARAMS)
    row = r["tumours"].row(0, named=True)
    return {
        "gds": r["gds"], "status": row["status"], "quality": row["quality"], "quality_code": row["quality_code"],
        "ploidy_prior": row["ploidy_prior"], "solution_strict": row["solution_strict"],
        "events": {e["wgd_event_id"]: [e["source"], e["event_quality"], e["event_quality_code"]]
                   for e in r["wgd_events"].iter_rows(named=True)},
        "modifications": [[m["modification"], m["samples"]] for m in r["modified_tumours"].iter_rows(named=True)],
    }


def test_every_case_has_an_agreed_answer():
    assert sorted(p.name for p in CASES.iterdir() if p.is_dir()) == sorted(AGREED)


@pytest.mark.parametrize("tumour", sorted(AGREED))
def test_agreed_answer(tumour):
    a, got = AGREED[tumour], run_case(tumour)
    assert got["gds"] == a["gds"], f"{tumour}: agreed {a['gds']} ({a['agreed']}); discuss with Piotr"
    assert got["status"] == a["status"]
    if "max_quality" in a:
        assert got["quality"] <= a["max_quality"]
    if "min_quality" in a:
        assert got["quality"] >= a["min_quality"]
    if "ploidy_prior" in a:
        assert (got["ploidy_prior"] or "").startswith(a["ploidy_prior"])


def test_clear_prior_scores_above_borderline_prior():
    """Piotr, 2026-11-12: LTX0981 (one borderline sample against a clear majority) is more
    believable than LTX0535 (every sample borderline); both clonal."""
    assert run_case("LTX0981-Tumour1")["quality"] > run_case("LTX0535-Tumour1")["quality"]


@pytest.mark.parametrize("tumour", sorted(AGREED))
def test_snapshot(tumour):
    f = CASES / tumour / "snapshot.json"
    got = run_case(tumour)
    if os.environ.get("UPDATE_CASES") or not f.exists():
        f.write_text(json.dumps(got, indent=1, sort_keys=True) + "\n")
    assert got == json.loads(f.read_text())


def test_ltx0272_gd_on_clone5_fails_as_reversion_to_diploid():
    """Piotr, 2026-11-12: reversion to diploid is impossible. GD1 forced onto clone5 (cluster 5
    alone) doubles clone12, 11, 10, 6, 8, which all sit in SU_T1.R2 (num_gds 0): failed, quality 0."""
    import polars as pl
    d = CASES / "LTX0272-Tumour1"
    t, ratios, cp = read_case(d / "cp_table.csv", d / "ratio_scores.csv")
    tb = pgdd_tables(d / "pgdd", t)
    tb["events"] = tb["events"].with_columns(
        pl.when(pl.col("GD_event") == "GD1").then(pl.lit("5")).otherwise(pl.col("clusters")).alias("clusters"))
    r = replay(tb, ratios, cp, PARAMS)
    row = r["tumours"].row(0, named=True)
    assert r["gds"] == "clone2=1,clone5=1"
    assert row["status"] == "failed" and row["quality"] == 0.0
    assert "reversion_to_diploid(clone12,clone11,clone8,clone10,clone6:R2(0<1))" in row["quality_code"]
    ev = {e["clone"]: e for e in r["wgd_events"].iter_rows(named=True)}
    assert ev["clone5"]["event_quality"] == 0.0 and ev["clone2"]["event_quality"] > 0
    assert r["modified_tumours"]["modification"].to_list() == ["failed_reversion_to_diploid"]


def test_ltx0552_minor_clone_pass():
    """Without the minor-clone pass only the majority clones are doubled (clone7, clone14, clone21)."""
    import dataclasses
    d = CASES / "LTX0552-Tumour2"
    t, ratios, cp = read_case(d / "cp_table.csv", d / "ratio_scores.csv")
    r = replay(pgdd_tables(d / "pgdd", t), ratios, cp, dataclasses.replace(PARAMS, fill_minor_clones=False))
    assert r["gds"] == "clone14=1,clone21=1,clone7=1"
    r = replay(pgdd_tables(d / "pgdd", t), ratios, cp, PARAMS)
    flags = {e["clone"]: e["flags"] for e in r["wgd_events"].iter_rows(named=True)}
    assert flags["clone20"] == "minor_clone_fill(R3:0.07,R4:0.04,R5:0.12)" and flags["clone14"] is None
