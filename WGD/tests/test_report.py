"""Step 07 report on the regression cases replayed as one cohort."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from wgd import params as P
from wgd import report
from wgd.replay import pgdd_tables, read_case, replay

CASES = Path(__file__).parent / "cases"


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    out = tmp_path_factory.mktemp("cohort")
    params = P.load(Path(__file__).resolve().parents[1] / "params.toml")
    parts: dict[str, list] = {}
    for t in ("LTX0981-Tumour1", "LTX1055-Tumour1", "LTX1348-Tumour1", "LTX0929-Tumour1", "LTX0025-Tumour1"):
        d = CASES / t
        tid, ratios, cp = read_case(d / "cp_table.csv", d / "ratio_scores.csv", t)
        r = replay(pgdd_tables(d / "pgdd", tid), ratios, cp, params)
        for k in ("events", "tumours", "wgd_events", "wgd_calls"):
            parts.setdefault(k, []).append(r[k])
    for sub, keys in (("04_placed", {"events": "events_placed", "tumours": "tumours_placed"}),
                      ("05_table", {"wgd_events": "wgd_events", "wgd_calls": "wgd_calls"})):
        (out / sub).mkdir()
        for k, name in keys.items():
            pl.concat(parts[k], how="diagonal_relaxed").write_csv(out / sub / f"{name}.tsv", separator="\t",
                                                                   null_value="NA")
    return report.run(out / "07_report", out / "04_placed", out / "05_table"), out


def test_patterns_and_files(cohort):
    res, out = cohort
    p = dict(res["tumours"].select("tumour_id", "pattern").iter_rows())
    assert p == {"LTX0981-Tumour1": "clonal only", "LTX1055-Tumour1": "subclonal, parallel",
                 "LTX1348-Tumour1": "subclonal, parallel", "LTX0929-Tumour1": "subclonal, parallel",
                 "LTX0025-Tumour1": "clonal + subclonal"}
    for f in ("report.html", "tumour_summary.tsv", "examples.tsv", "review.tsv", "quality_factors.tsv",
              "quality_by_pattern.png", "patterns.png", "gd_counts.png", "gd_sources.png", "quality_factors.png"):
        assert (out / "07_report" / f).stat().st_size > 0
    assert "data:image/png;base64," in (out / "07_report" / "report.html").read_text()


def test_examples_only_subclonal_with_commentary(cohort):
    ex = cohort[0]["examples"]
    assert "LTX0981-Tumour1" not in ex["tumour_id"].to_list()       # clonal only: not an example
    lo = ex.filter(pl.col("band") == "low").row(0, named=True)
    assert lo["tumour_id"] == "LTX0929-Tumour1" and "minor-clone GD" in lo["commentary"]
    hi = ex.filter(pl.col("tumour_id") == "LTX1348-Tumour1").row(0, named=True)
    assert hi["band"] == "high" and "alpaca_parent" in hi["commentary"]


def test_review_reasons(cohort):
    rv = cohort[0]["review"]
    got = {(t, r) for t, r in rv.select("tumour_id", "reason").iter_rows()}
    assert ("LTX1055-Tumour1", "forced GD (clone proportions vs ploidy)") in got
    assert ("LTX0929-Tumour1", "several GDs on one minor clone") in got
    assert ("LTX1348-Tumour1", "ALPACA CN to correct (parallel GDs)") in got
    assert rv["severity"][0] == "serious"                            # worst first


def test_explain_every_factor_code():
    from wgd import quality as Q
    codes = ["borderline_ploidy", "prior_raised", "prior_lowered", "alpaca_disagrees", "low_mut_support",
             "minority_coverage", "minor_clone", "reversion_to_diploid"] + list(Q.FIXED)
    sample = {"borderline_ploidy": "0.010:R1", "prior_raised": "R1:0.010", "prior_lowered": "R1:0.010",
              "alpaca_disagrees": "1.03:clone5:mut21", "low_mut_support": "8", "minority_coverage": "0.32:R1",
              "minor_clone": "0.15:R3"}
    for c in codes:
        assert not report.explain(c, 0.5, sample.get(c, "x")).startswith(f"x0.50 {c}(")


def test_no_calls_section(cohort, tmp_path):
    """Piotr, 2026-11-16: tumours failing QC or without an ALPACA tree are listed in the report, not counted."""
    _, out = cohort
    s05 = tmp_path / "05_table"
    s05.mkdir()
    (s05 / "wgd_events.tsv").write_bytes((out / "05_table" / "wgd_events.tsv").read_bytes())
    calls = pl.read_csv(out / "05_table" / "wgd_calls.tsv", separator="\t", null_values=["NA"], infer_schema_length=None)
    extra = pl.DataFrame({"tumour_id": ["LTX9998-Tumour1", "LTX9999-Tumour1", "LTX9999-Tumour1"],
                          "clone": [None, "clone1", "clone2"], "parent": [None, "diploid", "clone1"],
                          "status": ["no_alpaca_tree", "qc_failed", "qc_failed"],
                          "upstream_qc": ["pass", "fail", "fail"], "ccf_qc": ["pass", "pass", "pass"],
                          "truncal_min_ccf": [99.0, 97.5, 97.5],
                          "qc_reason": [None, "upstream: not FINALIZED", "upstream: not FINALIZED"]})
    pl.concat([calls, extra], how="diagonal_relaxed").write_csv(s05 / "wgd_calls.tsv", separator="\t", null_value="NA")
    res = report.run(tmp_path / "07_report", out / "04_placed", s05)
    nc = res["no_calls"]
    assert nc.select("tumour_id", "status", "reason").rows() == [
        ("LTX9998-Tumour1", "no_alpaca_tree", "no ALPACA tree"), ("LTX9999-Tumour1", "qc_failed", "upstream: not FINALIZED")]
    html = res["report"].read_text()
    assert "Tumours without WGD calls" in html and "LTX9999-Tumour1" in html and "LTX9998-Tumour1" in html
    assert (tmp_path / "07_report" / "no_calls.tsv").exists()
    assert not {"LTX9998-Tumour1", "LTX9999-Tumour1"} & set(res["tumours"]["tumour_id"])


def test_report_when_no_tumour_has_calls(tmp_path):
    """2026-11-18 (Piotr's run): every tumour qc_failed -> the report still runs and lists them."""
    s04, s05 = tmp_path / "04", tmp_path / "05"
    s04.mkdir()
    s05.mkdir()
    pl.DataFrame(schema={"tumour_id": pl.Utf8, "triage": pl.Utf8, "status": pl.Utf8, "quality": pl.Float64,
                         "quality_code": pl.Utf8, "ploidy_prior": pl.Utf8, "tumour_factors": pl.Utf8}).write_csv(
        s04 / "tumours_placed.tsv", separator="\t")
    pl.DataFrame(schema={"GD_event": pl.Utf8, "tumour_id": pl.Utf8, "flags": pl.Utf8}).write_csv(
        s04 / "events_placed.tsv", separator="\t")
    pl.DataFrame({"tumour_id": ["LTX9999-Tumour1"], "clone": ["clone1"], "parent": ["diploid"],
                  "number_of_gds_at_clone": [None], "status": ["qc_failed"], "class": [None], "upstream_qc": ["fail"],
                  "ccf_qc": ["pass"], "truncal_min_ccf": [99.0], "qc_pass": [False],
                  "qc_reason": ["upstream: not FINALIZED"]},
                 schema_overrides={"number_of_gds_at_clone": pl.Int64, "class": pl.Utf8}).write_csv(
        s05 / "wgd_calls.tsv", separator="\t", null_value="NA")
    pl.DataFrame(schema={"wgd_event_id": pl.Utf8, "tumour_id": pl.Utf8, "clone": pl.Utf8, "source": pl.Utf8,
                         "flags": pl.Utf8, "event_factors": pl.Utf8}).write_csv(s05 / "wgd_events.tsv", separator="\t")
    res = report.run(tmp_path / "07", s04, s05)
    assert res["tumours"].height == 0 and res["no_calls"]["tumour_id"].to_list() == ["LTX9999-Tumour1"]
    assert "LTX9999-Tumour1" in res["report"].read_text()
