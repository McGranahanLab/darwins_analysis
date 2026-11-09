"""Steps 04-06 on a small synthetic tumour, one test per notebook rule.

Tree (ALPACA):        c1 (MRCA) ── c2 ──┬── c4     R1 = 100% c4
                               │        └── c5     R2 = 100% c5
                               └── c3              R3 = 100% c3
"""

from __future__ import annotations

import polars as pl
import pytest

from wgd import compare, final_table, placement

T = "LTX0901-Tumour1"
PARENT = {"clone1": "diploid", "clone2": "clone1", "clone3": "clone1", "clone4": "clone2", "clone5": "clone2"}
CP = {"R1": "clone4", "R2": "clone5", "R3": "clone3"}


def ratios(**r) -> pl.DataFrame:
    base = {c: 1.0 for c in PARENT} | {f"clone{k[1:]}": v for k, v in r.items()}
    return pl.DataFrame({"tumour_id": [T] * 5, "clone": list(PARENT), "parent": list(PARENT.values()),
                         "ratio_A": [base[c] for c in PARENT], "ratio_B": [1.0] * 5})


def clone_props() -> pl.DataFrame:
    rows = [{"tumour_id": T, "clone": c, "sample_id": f"LTX0901_SU_T1-{s}--abc", "proportion": float(CP[s] == c)}
            for s in CP for c in PARENT]
    return pl.DataFrame(rows)


def sample(s: str) -> str:
    return f"LTX0901_SU_T1-{s}--abc"


def events(*evs) -> pl.DataFrame:
    """evs: (is_clonal, mut_supported, clusters, samples)."""
    return pl.DataFrame([{"GD_event": f"GD{i}", "GD_event_id": f"{T}_GD{i}", "tumour_id": T, "is_clonal": c,
                          "is_subclonal_mutation_supported": m, "clusters": cl,
                          "samples": ",".join(sample(s) for s in ss)}
                         for i, (c, m, cl, ss) in enumerate(evs, 1)],
                        schema_overrides={"clusters": pl.Utf8})


def tumour(first_gd="Subclonal", n_sub=1, n_second=0, n_clonal=0) -> pl.DataFrame:
    return pl.DataFrame({"tumour_id": [T], "First_GD": [first_gd], "num_subclonal_gds": [n_sub],
                         "num_second_gd": [n_second], "num_clonal_gds": [n_clonal]})


def mut_counts(**cn2) -> pl.DataFrame:
    rows = [{"tumour_id": T, "cluster_id": k[1:], "num_cn2": v, "discover_pass_all": True, "check_pass_all": False}
            for k, v in cn2.items()]
    return pl.DataFrame(rows, schema={"tumour_id": pl.Utf8, "cluster_id": pl.Utf8, "num_cn2": pl.Int64,
                                      "discover_pass_all": pl.Boolean, "check_pass_all": pl.Boolean})


def run(ev, tu=None, mc=None, rt=None, **kw):
    return placement.place_events(ev, tu if tu is not None else tumour(), mc if mc is not None else mut_counts(),
                                  rt if rt is not None else ratios(), clone_props(), **kw)


def placed(res) -> list[str | None]:
    return res.events["placed_as"].to_list()


def test_triage():
    ev = events((True, False, None, ["R1", "R2", "R3"]))
    assert placement.triage(tumour("No GD", 0), ev)[T] == "no_gd"
    assert placement.triage(tumour("Clonal", 0), ev)[T] == "clonal_gd_only"
    assert placement.triage(tumour("Clonal", 1), ev)[T] == "alpaca_called"
    ev2 = events((False, True, "4", ["R1"]))
    assert placement.triage(tumour(), ev2)[T] == "mut_supported"


def test_C_alpaca_single_event():
    res = run(events((False, False, None, ["R1", "R2"])), rt=ratios(c2=1.95))
    assert placed(res) == ["alpaca_clone2"]
    assert res.tumours["status"].item() == "resolved"


def test_C_greedy_file_order():
    # notebook behaviour kept on purpose (doc section 5, item 1)
    ev = events((False, False, None, ["R1"]), (False, False, None, ["R2"]))
    res = run(ev, tumour(n_sub=2), rt=ratios(c2=1.95, c4=1.6, c5=1.1))
    assert placed(res) == ["alpaca_clone2", "alpaca_clone5"]
    assert res.tumours["status"].item() == "unresolved"  # 1.1 < 1.5
    res = run(ev.reverse(), tumour(n_sub=2), rt=ratios(c2=1.95, c4=1.6, c5=1.1))
    assert placed(res) == ["alpaca_clone2", "alpaca_clone4"]
    assert res.tumours["status"].item() == "resolved"


def test_C_mrca_never_a_candidate_and_max_ab():
    ev = events((False, False, None, ["R3"]))
    res = run(ev, rt=ratios(c1=3.0, c3=1.2))
    assert placed(res) == ["alpaca_clone3"]
    rt = ratios(c3=1.2).with_columns(pl.when(pl.col("clone") == "clone3").then(1.8).otherwise(1.0).alias("ratio_B"))
    res = run(ev, rt=rt, score="max_AB")
    assert res.events["alpaca_score"].item() == 1.8 and res.tumours["status"].item() == "resolved"


def test_E_parent_child_to_child():
    res = run(events((False, True, "2,4", ["R1"]), (False, True, "4,5", ["R2"])))
    assert placed(res) == ["multicluster_clone4", "4,5"]
    assert res.events["placed_clone"].to_list() == ["clone4", None]
    assert res.events["flags"].to_list()[0] == "multicluster_parent_child"


def test_E_old_clone_prefix_format():
    res = run(events((False, True, "clone2,clone4", ["R1"])))
    assert placed(res) == ["multicluster_clone4"]


def test_F_clonal_to_mrca():
    res = run(events((True, False, None, ["R1", "R2", "R3"]), (False, True, "4", ["R1"])), tumour("Clonal", 1))
    assert placed(res) == ["clonal_clone1", "4"]


def test_G_most_ancestral_allowed_clone():
    # mutation-supported GD on clone5 forbids clone5 and its ancestors c2, c1
    res = run(events((False, True, "5", ["R2"]), (False, False, None, ["R1"])), tumour(n_sub=2))
    assert placed(res) == ["5", "ploidy_clone4"]


def test_G_climbs_when_allowed():
    # mutation-supported GD on clone3: forbids c3 and c1; R1's lineage c4 -> c2 -> c1 tops out at c2
    res = run(events((False, True, "3", ["R3"]), (False, False, None, ["R1"])), tumour(n_sub=2))
    assert placed(res) == ["3", "ploidy_clone2"]


def test_G_false_call_and_skips():
    # mutation GD on clone2 disallows its descendants c4, c5 and ancestor c1: nothing left in R1
    res = run(events((False, True, "2", ["R1", "R2"]), (False, False, None, ["R1"])), tumour(n_sub=2))
    assert placed(res) == ["2", "false_call"]
    # notebook quirk kept: a plain cluster label does not disallow its own clone
    res = run(events((False, True, "4", ["R1"]), (False, False, None, ["R1"])), tumour(n_sub=2))
    assert placed(res) == ["4", "ploidy_clone4"]
    res = run(events((False, True, "4", ["R1"]), (False, False, None, ["R2", "R3"])), tumour(n_sub=2))
    assert placed(res) == ["4", None] and res.events["flags"][1] == "G_skipped_multi_sample"
    res = run(events((False, True, "4", ["R1"]), (False, False, None, ["R2"])), tumour(n_sub=2, n_second=1))
    assert res.events["flags"][1] == "G_skipped_second_gd"
    res = run(events((False, True, "4", ["R1"]), (False, False, None, ["R2"]), (False, False, None, ["R3"])),
              tumour(n_sub=3))
    assert res.events["flags"].to_list()[1:] == ["G_skipped_other_ploidy_event"] * 2


def test_H_path_count_vs_ploidy():
    summ = pl.DataFrame({"tumour_id": [T] * 3, "sample_name_hash": ["a", "b", "c"], "num_gds": [2, 1, 0]})
    res = run(events((False, True, "5", ["R2"]), (False, False, None, ["R1"])), tumour(n_sub=2), summary=summ)
    t = res.tumours.row(0, named=True)
    assert t["H_mismatch"] and t["H_summary_max_gds"] == 2 and t["H_tree_path_max_gds"] == 1


def test_I_reversion():
    cn = pl.DataFrame({"tumour_id": [T] * 4, "clone": ["clone2", "clone2", "clone5", "clone5"],
                       "pred_CN_A": [2, 2, 1, 1]})
    res = run(events((False, True, "2", ["R1", "R2"])), alpaca_cn=cn)
    assert res.events["flags"].item() == "reversion"


def test_supporting_mut_count():
    mc = mut_counts(c4=12)
    mc = pl.concat([mc, mc.with_columns(pl.lit(30).cast(pl.Int64).alias("num_cn2"), pl.lit(False).alias("discover_pass_all"))])
    res = run(events((False, True, "4", ["R1"])), mc=mc)
    assert res.events["supporting_mut_count"].item() == 12  # the 30-row passes neither filter


def test_qc_filter():
    res = run(events((False, True, "4", ["R1"])), qc_tumours=set())
    assert res.tumours.is_empty()


def table(res, tu=None, rt=None):
    return final_table.build_table(res.events, res.tumours, tu if tu is not None else tumour(),
                                   rt if rt is not None else ratios())


def test_table_rows_and_counts():
    res = run(events((True, False, None, ["R1", "R2", "R3"]), (False, True, "4", ["R1"]), (False, True, "5", ["R2"])),
              tumour("Clonal", 2), mc=mut_counts(c4=20, c5=20))
    tb = table(res)
    assert tb.columns == list(final_table.COLUMNS)
    assert tb["clone"].to_list() == sorted(PARENT)
    n = dict(zip(tb["clone"], tb["number_of_gds_at_clone"]))
    assert n == {"clone1": 1, "clone2": 0, "clone3": 0, "clone4": 1, "clone5": 1}
    assert tb["total_gds_per_tumour"].unique().to_list() == [3]
    assert tb["class"].unique().to_list() == ["mut_supported"]
    assert tb["status"].unique().to_list() == ["unresolved"]
    assert tb["sub_class"].unique().to_list() == ["parallel only"]
    assert tb["parent"].to_list() == [PARENT[c] for c in sorted(PARENT)]


def test_table_clonal_only_uses_num_clonal_gds():
    res = run(events((True, False, None, ["R1", "R2", "R3"])), tumour("Clonal", 0, n_clonal=2))
    tb = table(res, tumour("Clonal", 0, n_clonal=2))
    assert dict(zip(tb["clone"], tb["number_of_gds_at_clone"]))["clone1"] == 2
    assert tb["class"][0] == "clonal_wgd" and tb["status"][0] == "resolved" and tb["sub_class"][0] is None


def test_sub_class():
    tree = placement.Tree(PARENT)
    assert final_table.sub_class(tree, ["clone4"]) == "single subclonal GD"
    assert final_table.sub_class(tree, ["clone2", "clone4"]) == "sequential only"
    assert final_table.sub_class(tree, ["clone4", "clone5"]) == "parallel only"
    assert final_table.sub_class(tree, ["clone2", "clone4", "clone3"]) == "parallel + sequential"


def test_quality_low_supporting_count():
    t2 = "LTX0902-Tumour1"
    ev = pl.DataFrame({"tumour_id": [T, t2], "GD_event_id": ["a", "b"], "is_clonal": [False, False],
                       "placed_clone": ["clone4", "clone4"], "supporting_mut_count": [3, 40]})
    tu = pl.DataFrame({"tumour_id": [T, t2], "triage": ["mut_supported"] * 2, "status": ["unresolved"] * 2,
                       "n_unplaced": [0, 0]})
    rt = pl.concat([ratios(), ratios().with_columns(pl.lit(t2).alias("tumour_id"))])
    pt = pl.DataFrame({"tumour_id": [T, t2], "num_clonal_gds": [0, 0]})
    tb = final_table.build_table(ev, tu, pt, rt)
    q = dict(tb.unique("tumour_id").select("tumour_id", "quality").iter_rows())
    assert q == {T: "Low", t2: "High"}
    assert tb.filter(pl.col("tumour_id") == T)["quality_code"][0] == "low_supporting_mut_count"


def test_compare():
    res = run(events((True, False, None, ["R1", "R2", "R3"])), tumour("Clonal", 0, n_clonal=1))
    new = table(res, tumour("Clonal", 0, n_clonal=1))
    old = new.select("tumour_id", "clone", "parent", "number_of_gds_at_clone", "total_gds_per_tumour",
                     "quality", "status", "class", "sub_class")
    s = compare.compare(old, new)["summary"]
    assert s["placement_match"] == {"match": 1, "of": 1} and s["class_match"]["match"] == 1
    moved = old.with_columns(pl.when(pl.col("clone") == "clone1").then(0).when(pl.col("clone") == "clone2").then(1)
                             .otherwise(0).alias("number_of_gds_at_clone"))
    r = compare.compare(moved, new)
    assert r["summary"]["placement_match"]["match"] == 0
    assert sorted(r["clone_diff"]["clone"].to_list()) == ["clone1", "clone2"]


def test_sample_key():
    assert placement.sample_key("LTX0001_SU_T1.R1--ab12") == "LTX0001_SU_T1-R1"
    with pytest.raises(ValueError, match="not in clone_proportions"):
        run(events((False, False, None, ["R9"])))
