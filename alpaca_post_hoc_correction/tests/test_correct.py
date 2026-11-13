"""Post-hoc ALPACA correction: parallel GDs promoted to an absent parent (LTX1348-like) and the
ratio port checked against tx_data's ALPACA ratios."""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from apc import correct

T = "LTX9001-Tumour1"
SEGS = ["1_0_100000000", "2_0_50000000", "3_0_50000000", "5_0_100000000"]


def alpaca(cn: dict[str, list[tuple[int, int]]]) -> pl.DataFrame:
    return pl.DataFrame([{"tumour_id": T, "segment": s, "clone": c, "pred_CN_A": a, "pred_CN_B": b}
                         for c, v in cn.items() for s, (a, b) in zip(SEGS, v)])


def calls(parent: dict[str, str], gds: dict[str, int], quality: float = 0.9) -> pl.DataFrame:
    return pl.DataFrame({"tumour_id": [T] * len(parent), "clone": list(parent), "parent": list(parent.values()),
                         "number_of_gds_at_clone": [gds.get(c, 0) for c in parent], "quality": [quality] * len(parent)})


def props(present: dict[str, float]) -> pl.DataFrame:
    return pl.DataFrame({"tumour_id": [T] * len(present), "clone": list(present),
                         "sample_id": ["s"] * len(present), "proportion": list(present.values())})


# LTX1348-like: clone1 (MRCA) absent; clone5 and clone6 doubled in parallel. True CN: clone1 = 1|1
# (2|1 on chr5, a gain before the GDs), children = 2 x clone1 with private changes. ALPACA, blind to
# the GDs, put the shared doubling on clone1.
PARENT = {"clone1": "diploid", "clone5": "clone1", "clone6": "clone1"}
GDS = {"clone5": 1, "clone6": 1}
ALPACA = {"clone1": [(2, 2), (2, 2), (2, 2), (4, 2)],
          "clone5": [(2, 2), (2, 2), (3, 2), (4, 2)],     # chr3: private gain after its GD
          "clone6": [(2, 2), (2, 1), (2, 2), (4, 2)]}     # chr2: private loss after its GD


def test_parallel_gds_absent_parent_halved():
    res = correct.run(alpaca(ALPACA), props({"clone5": 0.6, "clone6": 0.4}), calls(PARENT, GDS))
    cn = {(c, s): (a, b) for c, s, a, b in res["alpaca_corrected"].select("clone", "segment", "pred_CN_A",
                                                                         "pred_CN_B").iter_rows()}
    assert [cn[("clone1", s)] for s in SEGS] == [(1, 1), (1, 1), (1, 1), (2, 1)]
    assert [cn[("clone5", s)] for s in SEGS] == [(2, 2), (2, 2), (3, 2), (4, 2)]       # observed: unchanged
    c = res["corrections"].row(0, named=True)
    assert (c["clone"], c["children_with_gd"], c["parent_looked_doubled"]) == ("clone1", "clone5(1),clone6(1)", True)
    r = {x["clone"]: x for x in res["ratio_comparison"].iter_rows(named=True)}
    # length-weighted (chr1, chr5 100 Mb; chr2, chr3 50 Mb): clone1 looked doubled, clone5 didn't; now the reverse
    assert r["clone1"]["ratio_A"] == pytest.approx(800 / 300) and r["clone1"]["ratio_A_corrected"] == pytest.approx(400 / 300)
    assert r["clone5"]["ratio_A"] == pytest.approx(325 / 300) and r["clone5"]["ratio_A_corrected"] == pytest.approx(650 / 300)
    rows = res["alpaca_corrected"].filter(pl.col("clone") == "clone1")
    assert rows["corrected"].all() and rows["pred_CN_A_original"].to_list() == [2, 2, 2, 4]


def test_present_parent_or_no_gd_child_untouched():
    # parent has cells: observed, never changed
    res = correct.run(alpaca(ALPACA), props({"clone1": 0.2, "clone5": 0.5, "clone6": 0.3}), calls(PARENT, GDS))
    assert res["corrections"].is_empty() and not res["alpaca_corrected"]["corrected"].any()
    # no GD on any child: nothing to correct
    res = correct.run(alpaca(ALPACA), props({"clone5": 0.6, "clone6": 0.4}), calls(PARENT, {}))
    assert res["corrections"].is_empty()
    # the parent itself carries the GD (clonal): keep ALPACA
    res = correct.run(alpaca(ALPACA), props({"clone5": 0.6, "clone6": 0.4}), calls(PARENT, {"clone1": 1}))
    assert res["corrections"].is_empty()


def test_mixed_children_and_low_quality_skip():
    # clone1 absent with one doubled child (clone5) and one undoubled present child (clone7 = 2|2):
    # clone7 says clone1 is 2|2, clone5 (4|4 after its GD) agrees -> ALPACA's 2|2 is kept
    parent = {"clone1": "diploid", "clone5": "clone1", "clone7": "clone1"}
    cn = {"clone1": [(2, 2)] * 4, "clone5": [(4, 4)] * 4, "clone7": [(2, 2)] * 4}
    res = correct.run(alpaca(cn), props({"clone5": 0.5, "clone7": 0.5}), calls(parent, {"clone5": 1}))
    assert not res["alpaca_corrected"]["corrected"].any()
    # quality gate: below min_quality the tumour is skipped
    res = correct.run(alpaca(ALPACA), props({"clone5": 0.6, "clone6": 0.4}), calls(PARENT, GDS, 0.2), min_quality=0.3)
    assert res["skipped"].row(0) == (T, "call quality below 0.3") and res["corrections"].is_empty()


def test_nested_absent_ancestors_bottom_up():
    # clone1 -> clone2 (absent) -> {clone3 GD, clone4 GD}; clone1 present at 1|1 everywhere
    parent = {"clone1": "diploid", "clone2": "clone1", "clone3": "clone2", "clone4": "clone2"}
    cn = {"clone1": [(1, 1)] * 4, "clone2": [(2, 2)] * 4, "clone3": [(2, 2)] * 4, "clone4": [(2, 2)] * 4}
    res = correct.run(alpaca(cn), props({"clone1": 0.2, "clone3": 0.4, "clone4": 0.4}),
                      calls(parent, {"clone3": 1, "clone4": 1}))
    got = res["alpaca_corrected"].filter(pl.col("clone") == "clone2")
    assert got["pred_CN_A"].to_list() == [1] * 4 and got["pred_CN_B"].to_list() == [1] * 4


def test_chain_of_missing_clones_both_corrected():
    """Piotr, 2026-11-13: clone1 absent, its child clone2 absent, clone2's children clone3, clone4
    doubled in parallel: correct clone2 and clone1 (both missing), nothing else."""
    parent = {"clone1": "diploid", "clone2": "clone1", "clone3": "clone2", "clone4": "clone2"}
    cn = {"clone1": [(2, 2)] * 4, "clone2": [(2, 2)] * 4, "clone3": [(2, 2)] * 4, "clone4": [(2, 2)] * 4}
    res = correct.run(alpaca(cn), props({"clone3": 0.5, "clone4": 0.5}), calls(parent, {"clone3": 1, "clone4": 1}))
    got = res["alpaca_corrected"].filter(pl.col("corrected"))
    assert sorted(set(got["clone"])) == ["clone1", "clone2"]
    for c in ("clone1", "clone2"):
        g = res["alpaca_corrected"].filter(pl.col("clone") == c)
        assert g["pred_CN_A"].to_list() == [1] * 4 and g["pred_CN_B"].to_list() == [1] * 4
    assert res["corrections"]["clone"].to_list() == ["clone2", "clone1"]           # deepest first
    assert res["flags"].is_empty()
    # a present clone in the chain stops it: clone1 observed -> only clone2
    res = correct.run(alpaca(cn), props({"clone1": 0.1, "clone3": 0.5, "clone4": 0.4}),
                      calls(parent, {"clone3": 1, "clone4": 1}))
    assert res["corrections"]["clone"].to_list() == ["clone2"]


def test_flags_only():
    # GD called on clone5 but ALPACA (observed) shows no doubling vs clone1, which has cells: flagged, not fixed
    res = correct.run(alpaca(ALPACA), props({"clone1": 0.1, "clone5": 0.5, "clone6": 0.4}), calls(PARENT, GDS))
    f = {(r["clone"], r["flag"]) for r in res["flags"].iter_rows(named=True)}
    assert ("clone5", "gd_without_alpaca_doubling") in f and not res["alpaca_corrected"]["corrected"].any()
    # reversion: clone2 below doubled clone1 back at 1|1
    parent = {"clone1": "diploid", "clone2": "clone1"}
    cn = {"clone1": [(2, 2)] * 4, "clone2": [(1, 1)] * 4}
    res = correct.run(alpaca(cn), props({"clone1": 0.5, "clone2": 0.5}), calls(parent, {"clone1": 1}))
    assert [(r["clone"], r["flag"]) for r in res["flags"].iter_rows(named=True)] == [("clone2", "alpaca_reversion")]
    # after correcting the LTX1348-like case, nothing is left to flag
    res = correct.run(alpaca(ALPACA), props({"clone5": 0.6, "clone6": 0.4}), calls(PARENT, GDS))
    assert res["flags"].is_empty()


TX = Path(os.environ.get("TX_DATA_REPO", Path(__file__).resolve().parents[3] / "tx_data")) / "data"


@pytest.mark.skipif(not (TX / "alpaca.parquet").exists(), reason="tx_data not found")
def test_ratio_port_matches_alpaca_on_tx_data():
    a = pl.read_parquet(TX / "alpaca.parquet")
    ref = pl.read_parquet(TX / "alpaca_wgd_ratios.parquet")
    for (t,), g in a.group_by("tumour_id"):
        r = ref.filter(pl.col("tumour_id") == t)
        got = correct.wgd_ratios(g, dict(zip(r["clone"], r["parent"])))
        m = got.join(r, on=["tumour_id", "clone"], suffix="_ref")
        assert m.height == r.height
        for k in ("ratio_A", "ratio_B", "ratio_total", "weight_sum_A"):
            assert (m[k] - m[f"{k}_ref"]).abs().max() < 1e-9, (t, k)
