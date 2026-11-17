"""The WGD deck's worked examples match the pipeline code, and the deck builds clean."""

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
    nf = cases["nested_fill"]
    assert nf.gds() == {"clone2": 1, "clone4": 1} and nf.tumour["status"] == "resolved"
    old, new = cases["stacked_notebook"], cases["stacked_fill"]
    assert old.gds() == {"clone1": 1, "clone3": 1, "clone2": 1} and old.levels()["R2"] == 2      # R2 over-doubled
    assert old.tumour["status"] == "failed"                                                       # caught today
    assert new.gds() == {"clone1": 1, "clone3": 2} and new.tumour["status"] == "resolved"
    assert new.levels() == new.num_gds
    ps, pa = cases["prior_strict"], cases["prior_applied"]
    assert ps.tumour["status"] == "unresolved" and pa.gds() == {"clone1": 1}
    assert pa.tumour["ploidy_prior"].startswith("applied:SYN0001_SU_T1-R4 0>1") and pa.tumour["quality"] > ps.tumour["quality"]
    assert cases["relaxed"].gds() == {"clone4": 1} and "minority_coverage(0.40:SYN0001_SU_T1-R1)" in cases["relaxed"].tumour["quality_code"]
    rv = cases["reversion"].tumour
    assert rv["status"] == "failed" and rv["quality"] == 0 and rv["quality_code"].startswith("reversion_to_diploid")
    fw = cases["fill_walk"]
    steps = ex.fill_steps(fw)
    assert [st["added"] for st in steps] == [None, "clone2", "clone4"] and fw.tumour["status"] == "resolved"
    assert steps[-1]["levels"] == fw.num_gds
    for st in steps[1:]:   # the chosen clone pushes no sample above target and raises the most samples
        ok = [c for c in st["candidates"] if not c["over"]]
        best = max(len(c["raised"]) for c in ok)
        assert next(c for c in ok if c["clone"] == st["added"])["raised"].__len__() == best
    fx = cases["fixed"]
    assert fx.gds() == {"clone1": 1, "clone4": 1} and fx.tumour["quality_code"] == "multicluster_parent_child"


def test_deck_builds(tmp_path):
    prs = Presentation(make_deck.build(tmp_path / "deck.pptx"))
    assert len(prs.slides) == 18
    for s in prs.slides:
        assert s.has_notes_slide and s.notes_slide.notes_text_frame.text.strip()
        txt = " ".join(sh.text_frame.text for sh in s.shapes if sh.has_text_frame)
        txt += " ".join(c.text for sh in s.shapes if sh.has_table for r in sh.table.rows for c in r.cells)
        assert not re.search(r"\b(LTX|CRUK)\d", txt), "real-looking tumour id on a slide"
