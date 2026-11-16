"""The deck's synthetic examples agree with the correction code, and the deck builds clean."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from pptx import Presentation

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import examples as ex  # noqa: E402
import make_deck  # noqa: E402


@pytest.fixture(scope="module")
def cases():
    return ex.all_cases()


def test_outcomes_shown_on_the_slides(cases):
    assert cases["A"].corrected_clones() == ["clone1"]
    assert cases["B"].candidates() == ["clone1"] and cases["B"].corrected_clones() == []
    assert cases["C"].corrected_clones() == ["clone2", "clone1"]          # deepest first
    assert cases["C2"].corrected_clones() == ["clone2"]
    for k in ("D1", "D2", "D3", "E"):
        assert cases[k].candidates() == []
    assert set(cases["D1"].flags()) == {("clone5", "gd_without_alpaca_doubling"),
                                        ("clone6", "gd_without_alpaca_doubling")}
    assert cases["E"].flags() == [("clone2", "alpaca_reversion")]


def test_example_a_numbers(cases):
    a = cases["A"]
    assert a.cn_after("clone1") == [(1, 1), (1, 1), (1, 1), (2, 1)]
    assert [round(a.ratio(c, False), 2) for c in a.parent] == [2.67, 1.08, 1.0]
    assert [round(a.ratio(c, True), 2) for c in a.parent] == [1.33, 2.17, 2.0]


def test_cost_tables_match_best_cn(cases):
    """The cost columns on the slides pick the same x as the correction itself."""
    a, c = cases["A"], cases["C"]
    for allele in (0, 1):
        rows = ex.cost_table(a, "clone1", allele)
        assert [r["new"] for r in rows] == [v[allele] for v in a.cn_after("clone1")]
        for r in rows:
            assert r["costs"].index(min(r["costs"])) == r["new"]
    after = {k: c.cn_after(k) for k in c.parent}
    assert ex.cost_table(c, "clone2", 0)[0]["new"] == after["clone2"][0][0]
    assert ex.cost_table(c, "clone1", 0, child_cn={"clone2": after["clone2"]})[0]["new"] == after["clone1"][0][0]
    tie = ex.cost_table(c, "clone1", 0)[0]["costs"]                       # clone1 first: x = 1 and 2 tie
    assert tie[1] == tie[2] == min(tie)


def test_deck_builds(tmp_path):
    out = make_deck.build(tmp_path / "deck.pptx")
    prs = Presentation(out)
    assert len(prs.slides) == 12
    for s in prs.slides:
        assert s.has_notes_slide and s.notes_slide.notes_text_frame.text.strip()
        txt = " ".join(sh.text_frame.text for sh in s.shapes if sh.has_text_frame)
        assert not re.search(r"\b(LTX|CRUK)\d", txt), "real-looking tumour id on a slide"
