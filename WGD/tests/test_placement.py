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
    """Notebook rules (C/G) unless ploidy_rule is given."""
    kw.setdefault("ploidy_rule", "notebook")
    return placement.place_events(ev, tu if tu is not None else tumour(), mc if mc is not None else mut_counts(),
                                  rt if rt is not None else ratios(), clone_props(), **kw)


def regions_df(t: str, gds: dict[str, int]) -> pl.DataFrame:
    return pl.DataFrame({"tumour_id": [t] * len(gds), "sample_id": [sample(s) for s in gds],
                         "num_gds": list(gds.values())})


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
    res = run(events((False, False, None, ["R1", "R2"])), rt=ratios(c2=1.95), ploidy_rule="notebook")
    assert placed(res) == ["alpaca_clone2"]
    assert res.tumours["status"].item() == "resolved"


def test_C_greedy_file_order():
    # notebook behaviour kept on purpose (doc section 5, item 1)
    ev = events((False, False, None, ["R1"]), (False, False, None, ["R2"]))
    res = run(ev, tumour(n_sub=2), rt=ratios(c2=1.95, c4=1.6, c5=1.1), ploidy_rule="notebook")
    assert placed(res) == ["alpaca_clone2", "alpaca_clone5"]
    assert res.tumours["status"].item() == "unresolved"  # 1.1 < 1.5
    res = run(ev.reverse(), tumour(n_sub=2), rt=ratios(c2=1.95, c4=1.6, c5=1.1), ploidy_rule="notebook")
    assert placed(res) == ["alpaca_clone2", "alpaca_clone4"]
    assert res.tumours["status"].item() == "resolved"


def test_C_mrca_never_a_candidate_and_max_ab():
    ev = events((False, False, None, ["R3"]))
    res = run(ev, rt=ratios(c1=3.0, c3=1.2), ploidy_rule="notebook")
    assert placed(res) == ["alpaca_clone3"]
    rt = ratios(c3=1.2).with_columns(pl.when(pl.col("clone") == "clone3").then(1.8).otherwise(1.0).alias("ratio_B"))
    res = run(ev, rt=rt, score="max_AB", ploidy_rule="notebook")
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


def test_quality_helpers():
    from wgd import quality as Q
    assert Q.graded(0.0, 0.1) == 0.5 and Q.graded(0.05, 0.1) == 0.75 and Q.graded(0.3, 0.1) == 1.0
    assert Q.graded(1.25, 1.5, 1.0) == 0.75
    fr = {"a": (0.507, 0.095, 0.03), "b": (0.42, 0.07, 0.03), "c": (0.9, 0.2, 0.05)}
    m = Q.ploidy_margins({"a": 1, "b": 0, "c": 1}, fr)
    assert m == pytest.approx({"a": 0.007, "b": 0.08, "c": 0.3})
    assert Q.flips({"a": 1, "b": 0, "c": 1}, fr, m) == [{"a": 0, "b": 0, "c": 1}, {"a": 1, "b": 1, "c": 1},
                                                         {"a": 0, "b": 1, "c": 1}]
    assert Q.score([("x", 0.5, "1"), ("y", 1.0, ""), ("z", 0.8, "")]) == (0.4, "x(1);z")


def summary_df(t: str, fr: dict[str, tuple]) -> pl.DataFrame:
    return pl.DataFrame({"tumour_id": [t] * len(fr), "sample_name_hash": [sample(s) for s in fr],
                         "genome_frac_nMaj2": [v[0] for v in fr.values()], "genome_frac_nMaj3": [v[1] for v in fr.values()],
                         "genome_frac_nMaj4": [v[2] for v in fr.values()],
                         "num_gds": [int(v[0] > 0.5) for v in fr.values()]})


def test_quality_borderline_and_ambiguous():
    """LTX0535-like (Piotr, 2026-11-10): genome fractions all near 0.5. R1, R2 just above
    (one subclonal GD on clone2); R3 at 0.49. Flipping R3 gives a clonal GD instead."""
    reg = regions_df(T, {"R1": 1, "R2": 1, "R3": 0})
    summ = summary_df(T, {"R1": (0.51, 0.1, 0.03), "R2": (0.53, 0.1, 0.03), "R3": (0.49, 0.1, 0.03)})
    res = run(events((False, False, None, ["R1", "R2"])), ploidy_rule="fill", regions=reg, summary=summ,
              rt=ratios(c2=2.0))
    t = res.tumours.row(0, named=True)
    assert placed(res)[-1] == "fill_clone2"
    assert t["quality_code"].startswith("borderline_ploidy(0.010:R")
    # R1 down -> clone5 (R2 only); R2 down -> clone4 (R1 only); R3 up -> one clonal GD
    assert "ambiguous_solution(alt:clone5=1|clone4=1|clone1=1)" in t["quality_code"]
    assert t["quality"] == pytest.approx(0.55 * 0.5, abs=1e-3)
    tb = final_table.build_table(res.events, res.tumours, tumour(), ratios(c2=2.0))
    assert tb["quality"].unique().to_list() == [t["quality"]]


def test_quality_clear_case_and_reversion():
    reg = regions_df(T, {"R1": 1, "R2": 1, "R3": 0})
    summ = summary_df(T, {"R1": (0.95, 0.1, 0.03), "R2": (0.9, 0.1, 0.03), "R3": (0.1, 0.05, 0.0)})
    kw = dict(ploidy_rule="fill", regions=reg, summary=summ, rt=ratios(c2=2.0))
    res = run(events((False, False, None, ["R1", "R2"])), **kw)
    assert res.tumours["quality"].item() == 1.0 and res.tumours["quality_code"].item() is None
    # clone4 (below the doubled clone2) looks diploid again
    cn = pl.DataFrame({"tumour_id": [T] * 4, "clone": ["clone2", "clone2", "clone4", "clone4"],
                       "pred_CN_A": [2, 2, 1, 1], "pred_CN_B": [2, 2, 1, 1]})
    res = run(events((False, False, None, ["R1", "R2"])), alpaca_cn=cn, **kw)
    assert res.tumours["quality"].item() == 0.2 and res.tumours["quality_code"].item() == "reversion(clone4)"


def test_quality_low_mut_support_and_alpaca():
    reg = regions_df(T, {"R1": 1, "R2": 0, "R3": 0})
    summ = summary_df(T, {"R1": (0.95, 0.1, 0.03), "R2": (0.1, 0.05, 0.0), "R3": (0.1, 0.05, 0.0)})
    res = run(events((False, True, "4", ["R1"])), mc=mut_counts(c4=5), ploidy_rule="fill", regions=reg,
              summary=summ, rt=ratios(c4=1.25))
    q = res.tumours.row(0, named=True)
    assert q["quality_code"] == "low_mut_support(5);alpaca_disagrees(1.25:clone4)"
    assert q["quality"] == pytest.approx(0.75 * 0.75, abs=1e-3)


def test_compare():
    res = run(events((True, False, None, ["R1", "R2", "R3"])), tumour("Clonal", 0, n_clonal=1))
    new = table(res, tumour("Clonal", 0, n_clonal=1))
    old = new.select("tumour_id", "clone", "parent", "number_of_gds_at_clone", "total_gds_per_tumour",
                     "quality", "status", "class", "sub_class")
    s = compare.compare(old, new)["summary"]
    assert s["placement_match"] == {"match": 1, "of": 1} and s["class_match"]["match"] == 1
    assert s["placement_match_by_new_status"] == {"resolved": {"match": 1, "of": 1}}
    moved = old.with_columns(pl.when(pl.col("clone") == "clone1").then(0).when(pl.col("clone") == "clone2").then(1)
                             .otherwise(0).alias("number_of_gds_at_clone"))
    r = compare.compare(moved, new)
    assert r["summary"]["placement_match"]["match"] == 0
    assert sorted(r["clone_diff"]["clone"].to_list()) == ["clone1", "clone2"]


def test_sample_key():
    assert placement.sample_key("LTX0001_SU_T1.R1--ab12") == "LTX0001_SU_T1-R1"
    assert placement.sample_key("LTX0025_SU_T1.R4..267e0b2fh2yt") == "LTX0025_SU_T1-R4"  # old PGDD, R-mangled
    assert placement.sample_key("LTX0025_SU_T1-R4--267e0b2fh2yt") == "LTX0025_SU_T1-R4"
    with pytest.raises(ValueError, match="not in clone_proportions"):
        run(events((False, False, None, ["R9"])))


def test_qc_pass_reads_na(tmp_path):
    qc = tmp_path / "qc.tsv"
    qc.write_text("patient\ttumour_id\tselected_solution\tFINALIZED\n"
                  "LTX0901\t1\tTRUE\tTRUE\nLTX0902\t1\tTRUE\tTRUE\nLTX0903\t1\tTRUE\tFALSE\n")
    for p, ccf in (("LTX0901", "103"), ("LTX0902", "85")):
        d = tmp_path / p / "mutation_trees" / "tumour_1"
        d.mkdir(parents=True)
        (d / "clusterInfo.txt").write_text(
            "clusterID\ttruncal\ttreeClust\tSAMPLE\tmeanCCF\n"
            f"1\tTRUE\tTRUE\tR1\t{ccf}\n1\tTRUE\tTRUE\tR2\t95\n2\tNA\tNA\tR1\tNA\n3\tTRUE\tFALSE\tR1\tNA\n"
            "4\tFALSE\tTRUE\tR1\t0\n")
    tx_qc, passed = placement.qc_pass(qc, tmp_path)
    assert tx_qc == {"LTX0901-Tumour1", "LTX0902-Tumour1"}
    assert passed == {"LTX0901-Tumour1"}
    # a truncal tree cluster without meanCCF is an error
    f = tmp_path / "LTX0901" / "mutation_trees" / "tumour_1" / "clusterInfo.txt"
    f.write_text("clusterID\ttruncal\ttreeClust\tSAMPLE\tmeanCCF\n1\tTRUE\tTRUE\tR1\tNA\n")
    with pytest.raises(ValueError, match="missing meanCCF"):
        placement.qc_pass(qc, tmp_path)


def test_ltx0025_like_greedy_case():
    """Real case (Piotr, 2026-11-09): clonal GD + two more GDs in R4 (clones 3, 5 only).
    Correct: 2 GDs on clone3. Notebook rule C picks clone3 (1.48) then clone2 (1.41),
    because a clone already taken can't be reused; both < 1.5, so the tumour is unresolved."""
    parent = {"clone6": "clone2", "clone2": "clone1", "clone3": "clone2", "clone1": "diploid",
              "clone4": "clone6", "clone7": "clone1", "clone5": "clone3"}
    ra = {"clone6": 1.0, "clone2": 1.406, "clone3": 1.485, "clone1": 1.948, "clone4": 0.96, "clone7": 1.0, "clone5": 1.05}
    rt = pl.DataFrame({"tumour_id": [T] * 7, "clone": list(parent), "parent": list(parent.values()),
                       "ratio_A": [ra[c] for c in parent], "ratio_B": [1.0] * 7})
    cp = {"R1": {"clone7": 1.0}, "R2": {"clone4": 0.62, "clone6": 0.38}, "R4": {"clone3": 0.6, "clone5": 0.4}}
    cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s), "proportion": cp[s].get(c, 0.0)}
                         for s in cp for c in parent])
    ev = events((True, False, None, ["R1", "R2", "R4"]), (False, False, None, ["R4"]), (False, False, None, ["R4"]))
    tu = tumour("Clonal", 2, n_second=1, n_clonal=1)
    res = placement.place_events(ev, tu, mut_counts(), rt, cpdf, ploidy_rule="notebook")
    assert placed(res) == ["clonal_clone1", "alpaca_clone3", "alpaca_clone2"]
    assert res.tumours["status"].item() == "unresolved"
    # fill rule: R4's two extra GDs both go on clone3, resolved automatically
    reg = regions_df(T, {"R1": 1, "R2": 1, "R4": 3})
    res = placement.place_events(ev, tu, mut_counts(), rt, cpdf, ploidy_rule="fill", regions=reg)
    assert placed(res) == ["clonal_clone1", "superseded_by_fill", "superseded_by_fill", "fill_clone3", "fill_clone3"]
    assert res.tumours["status"].item() == "resolved" and res.tumours["fill"].item() == "complete"
    tb = final_table.build_table(res.events, res.tumours, tu, rt)
    assert dict(zip(tb["clone"], tb["number_of_gds_at_clone"])) == {
        "clone1": 1, "clone2": 0, "clone3": 2, "clone4": 0, "clone5": 0, "clone6": 0, "clone7": 0}
    assert tb["total_gds_per_tumour"][0] == 3 and tb["sub_class"][0] == "sequential only"
    assert tb["n_unplaced_events"][0] == 0


def test_ltx0073_mutation_plus_ploidy():
    """Real case (Piotr, 2026-11-09): clonal GD, mutation-supported GD on clone3, and PGDD's
    ploidy-only events GD3 (R1+R4), GD4 (R3), GD5 (R1). R4's 2nd GD is already explained by
    clone3 (63% of R4). Expected: clone1 = 1, clone3 = 1, clone7 = 2, clone14 = 1."""
    parent = {"clone15": "clone7", "clone4": "clone3", "clone11": "clone4", "clone5": "clone11", "clone1": "diploid",
              "clone12": "clone5", "clone7": "clone1", "clone14": "clone5", "clone16": "clone14", "clone3": "clone8",
              "clone6": "clone4", "clone8": "clone1"}
    ra = {"clone15": 1.085, "clone4": 0.977, "clone11": 1.063, "clone5": 1.057, "clone1": 3.304, "clone12": 1.003,
          "clone7": 1.138, "clone14": 1.106, "clone16": 1.029, "clone3": 0.984, "clone6": 0.908, "clone8": 1.002}
    rt = pl.DataFrame({"tumour_id": [T] * 12, "clone": list(parent), "parent": list(parent.values()),
                       "ratio_A": [ra[c] for c in parent], "ratio_B": [1.0] * 12})
    cp = {"R1": {"clone7": 0.56, "clone15": 0.44}, "R2": {"clone6": 0.7787, "clone11": 0.2213},
          "R3": {"clone14": 0.52, "clone16": 0.48}, "R4": {"clone3": 0.32, "clone8": 0.37, "clone11": 0.01, "clone12": 0.30}}
    cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s), "proportion": cp[s].get(c, 0.0)}
                         for s in cp for c in parent])
    ev = events((True, False, None, ["R1", "R2", "R3", "R4"]), (False, True, "3", ["R2", "R3"]),
                (False, False, None, ["R1", "R4"]), (False, False, None, ["R3"]), (False, False, None, ["R1"]))
    tu = tumour("Clonal", 4, n_second=2, n_clonal=1)
    reg = regions_df(T, {"R1": 3, "R2": 2, "R3": 3, "R4": 2})
    res = placement.place_events(ev, tu, mut_counts(c3=11), rt, cpdf, ploidy_rule="fill", regions=reg)
    tb = final_table.build_table(res.events, res.tumours, tu, rt)
    n = {c: v for c, v in zip(tb["clone"], tb["number_of_gds_at_clone"]) if v}
    assert n == {"clone1": 1, "clone3": 1, "clone7": 2, "clone14": 1}
    assert res.tumours["status"].item() == "resolved"
    # clone5 would also bring R3 to 3, but it covers clone12 (30% of R4): more spill than clone14
    assert tb["total_gds_per_tumour"][0] == 5 and tb["class"][0] == "mut_supported"


def test_sample_levels_and_fill_limits():
    tree = placement.Tree(PARENT)
    props = {s: {CP[s]: 1.0} for s in CP}
    assert placement.sample_levels(tree, props, {"clone1": 1, "clone2": 1}) == {"R1": 2, "R2": 2, "R3": 1}
    # clone4 is only 40% of R1 and clone2 would over-double R2: the relaxed pass takes clone4
    props2 = {**props, "R1": {"clone4": 0.4, "clone2": 0.6}}
    f = placement.fill_ploidy(tree, props2, {"R1": 2, "R2": 1, "R3": 1}, {"clone1": 1}, {})
    assert [(a["clone"], a.get("relaxed")) for a in f.added] == [("clone4", True)] and f.remaining == {}
    assert f.added[0]["coverage"] == {"R1": 0.4}
    # over-explained: a fixed GD on clone2 pushes R2 above its target
    f = placement.fill_ploidy(tree, props, {"R1": 2, "R2": 1, "R3": 1}, {"clone1": 1, "clone2": 1}, {})
    assert f.over == {"R2": 1} and not f.added


def test_fill_resolution_and_status():
    reg = regions_df(T, {"R1": 2, "R2": 2, "R3": 1})
    # one ploidy-only event in R1 + R2: fill puts it on clone2
    res = run(events((True, False, None, ["R1", "R2", "R3"]), (False, False, None, ["R1", "R2"])),
              tumour("Clonal", 1, n_clonal=1), ploidy_rule="fill", regions=reg)
    assert placed(res)[-1] == "fill_clone2" and res.tumours["status"].item() == "resolved"
    # R1 and R3 doubled but not R2: no single clone, so two parallel GDs (clone4, clone3)
    reg2 = regions_df(T, {"R1": 2, "R2": 1, "R3": 2})
    res = run(events((True, False, None, ["R1", "R2", "R3"]), (False, False, None, ["R1", "R3"])),
              tumour("Clonal", 1, n_clonal=1), ploidy_rule="fill", regions=reg2)
    assert sorted(placed(res)[2:]) == ["fill_clone3", "fill_clone4"]
    assert res.tumours["status"].item() == "resolved"


def test_ltx0501_prefers_ancestral_clone_without_overdoubling():
    """Real case (Piotr, 2026-11-10): mutation-supported GD on clone4 (R1); PGDD's ploidy-only
    GD1 covers R2 + R3. Correct (old, hand-curated): clone4, clone5, clone6. clone9 would also
    raise R2, but clone6 additionally doubles its 8% share of R1, which is at 1 GD already."""
    parent = {"clone8": "clone7", "clone1": "diploid", "clone4": "clone1", "clone7": "clone9", "clone5": "clone1",
              "clone3": "clone4", "clone9": "clone6", "clone2": "clone5", "clone6": "clone1"}
    ra = {"clone8": 0.954, "clone1": 1.947, "clone4": 1.115, "clone7": 1.008, "clone5": 1.011, "clone3": 0.986,
          "clone9": 1.053, "clone2": 0.928, "clone6": 1.030}
    rt = pl.DataFrame({"tumour_id": [T] * 9, "clone": list(parent), "parent": list(parent.values()),
                       "ratio_A": [ra[c] for c in parent], "ratio_B": [1.0] * 9})
    cp = {"R1": {"clone3": 0.2, "clone4": 0.7174, "clone6": 0.0826}, "R2": {"clone7": 0.32, "clone8": 0.34, "clone9": 0.34},
          "R3": {"clone2": 0.39, "clone5": 0.61}}
    cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s), "proportion": cp[s].get(c, 0.0)}
                         for s in cp for c in parent])
    ev = events((False, False, None, ["R2", "R3"]), (False, True, "4", ["R1"]))
    tu = tumour("Subclonal", 2)
    reg = regions_df(T, {"R1": 1, "R2": 1, "R3": 1})
    res = placement.place_events(ev, tu, mut_counts(c4=12), rt, cpdf, ploidy_rule="fill", regions=reg)
    tb = final_table.build_table(res.events, res.tumours, tu, rt)
    assert {c: v for c, v in zip(tb["clone"], tb["number_of_gds_at_clone"]) if v} == {"clone4": 1, "clone5": 1, "clone6": 1}
    assert res.tumours["status"].item() == "resolved"
    assert tb["sub_class"][0] == "parallel only"


def test_ltx1338_goes_as_high_as_the_data_allow():
    """Real case (Piotr, 2026-11-10): mutation-supported GD on clone2 (explains R1, R2 and R3
    through clone11); R4 still needs one. clone9 (extinct) and its only child clone15 cover
    the same cells, so the GD goes on the more ancestral clone9, not clone15."""
    parent = {"clone1": "diploid", "clone14": "clone15", "clone6": "clone7", "clone7": "clone12", "clone5": "clone14",
              "clone15": "clone9", "clone10": "clone7", "clone4": "clone2", "clone11": "clone12", "clone9": "clone1",
              "clone12": "clone2", "clone2": "clone1"}
    ra = {"clone1": 1.957, "clone14": 1.077, "clone6": 0.989, "clone7": 1.068, "clone5": 1.029, "clone15": 1.0,
          "clone10": 0.984, "clone4": 0.996, "clone11": 1.026, "clone9": 0.910, "clone12": 1.012, "clone2": 1.071}
    rt = pl.DataFrame({"tumour_id": [T] * 12, "clone": list(parent), "parent": list(parent.values()),
                       "ratio_A": [ra[c] for c in parent], "ratio_B": [1.0] * 12})
    cp = {"R1": {"clone6": 1.0}, "R2": {"clone10": 1.0}, "R3": {"clone1": 0.35, "clone11": 0.65},
          "R4": {"clone4": 0.2424, "clone5": 0.45, "clone14": 0.22, "clone15": 0.0876}}
    cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s), "proportion": cp[s].get(c, 0.0)}
                         for s in cp for c in parent])
    ev = events((False, True, "2", ["R1", "R2"]), (False, False, None, ["R3", "R4"]))
    tu = tumour("Subclonal", 2)
    reg = regions_df(T, {"R1": 1, "R2": 1, "R3": 1, "R4": 1})
    res = placement.place_events(ev, tu, mut_counts(c2=30), rt, cpdf, ploidy_rule="fill", regions=reg)
    tb = final_table.build_table(res.events, res.tumours, tu, rt)
    assert {c: v for c, v in zip(tb["clone"], tb["number_of_gds_at_clone"]) if v} == {"clone2": 1, "clone9": 1}
    assert res.tumours["status"].item() == "resolved"
    assert "fill_tie" not in (res.events["flags"].drop_nulls().to_list())


def test_ltx1055_relaxed_minority_clone():
    """Real case (Piotr, 2026-11-11): mutation-supported GD on clone4 (FR1, FR2). SU R1 is fully
    doubled (genome_frac_nMaj2 0.9997) but cp_table gives it 68% clone2, which is undoubled in
    R2-R7, and only 32% clone3. No clone is the majority of R1 without doubling R2-R7, so the
    relaxed pass puts the GD on clone3 (ALPACA ratio 3.2) and lowers quality."""
    parent = {"clone1": "diploid", "clone2": "clone1", "clone3": "clone2", "clone6": "clone4", "clone4": "clone1"}
    ra = {"clone1": 1.033, "clone2": 1.0, "clone3": 3.205, "clone6": 0.987, "clone4": 2.36}
    rt = pl.DataFrame({"tumour_id": [T] * 5, "clone": list(parent), "parent": list(parent.values()),
                       "ratio_A": [ra[c] for c in parent], "ratio_B": [1.0] * 5})
    cp = {"F1": {"clone4": 1.0}, "F2": {"clone4": 0.38, "clone6": 0.62}, "R1": {"clone2": 0.68, "clone3": 0.32},
          **{f"R{i}": {"clone2": 1.0} for i in range(2, 8)}}
    cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s), "proportion": cp[s].get(c, 0.0)}
                         for s in cp for c in parent])
    ev = events((False, True, "4", ["F1", "F2"]), (False, False, None, ["R1"]))
    tu = tumour("Subclonal", 2)
    reg = regions_df(T, {"F1": 1, "F2": 1, "R1": 1, **{f"R{i}": 0 for i in range(2, 8)}})
    res = placement.place_events(ev, tu, mut_counts(c4=13), rt, cpdf, ploidy_rule="fill", regions=reg)
    tb = final_table.build_table(res.events, res.tumours, tu, rt)
    assert {c: v for c, v in zip(tb["clone"], tb["number_of_gds_at_clone"]) if v} == {"clone3": 1, "clone4": 1}
    t = res.tumours.row(0, named=True)
    assert t["status"] == "resolved"
    assert "minority_coverage(0.32:R1)" in t["quality_code"]
    assert t["quality"] == pytest.approx(0.2 + 0.3 * 0.32 / 0.5, abs=1e-3)  # forced placement: heavy penalty


def test_event_quality_per_clone():
    """Two GDs in one tumour: a well-supported one (clone4) and a weak one (clone5, 5 mutations)."""
    reg = regions_df(T, {"R1": 1, "R2": 1, "R3": 0})
    summ = summary_df(T, {"R1": (0.95, 0.1, 0.03), "R2": (0.9, 0.1, 0.03), "R3": (0.1, 0.05, 0.0)})
    ev = events((False, True, "4", ["R1"]), (False, True, "5", ["R2"]))
    res = run(ev, tumour(n_sub=2), mc=mut_counts(c4=40, c5=5), ploidy_rule="fill", regions=reg, summary=summ,
              rt=ratios(c4=2.0, c5=2.0))
    e = {r["placed_clone"]: r for r in res.events.iter_rows(named=True)}
    assert e["clone4"]["event_quality"] == 1.0 and e["clone4"]["event_quality_code"] is None
    assert e["clone5"]["event_quality"] == 0.75 and e["clone5"]["event_quality_code"] == "low_mut_support(5)"
    assert res.tumours["quality"].item() == 0.75
    tb = final_table.build_table(res.events, res.tumours, tumour(n_sub=2), ratios(c4=2.0, c5=2.0))
    q = dict(zip(tb["clone"], tb["worst_event_quality"]))
    assert q["clone4"] == 1.0 and q["clone5"] == 0.75 and q["clone2"] is None
    assert tb["quality"].unique().to_list() == [0.75]


def test_ltx1018_forced_third_gd_on_clone2():
    """Real case (Piotr, 2026-11-11): clonal GD + mutation-supported GD on clone3. R2 needs a
    third GD; clone2 (+ clone7) is only 34% of R2, and clone5 would also double R1's clone8 line
    (over-doubling). Place it on clone2, but penalise heavily."""
    parent = {"clone3": "clone6", "clone15": "clone13", "clone12": "clone3", "clone11": "clone12", "clone5": "clone1",
              "clone8": "clone5", "clone7": "clone2", "clone6": "clone1", "clone4": "clone11", "clone13": "clone8",
              "clone2": "clone12", "clone9": "clone4", "clone1": "diploid"}
    ra = {"clone3": 1.227, "clone15": 1.0, "clone12": 1.223, "clone11": 0.999, "clone5": 1.019, "clone8": 0.998,
          "clone7": 1.001, "clone6": 1.352, "clone4": 1.002, "clone13": 1.0, "clone2": 1.043, "clone9": 1.0, "clone1": 2.177}
    rt = pl.DataFrame({"tumour_id": [T] * 13, "clone": list(parent), "parent": list(parent.values()),
                       "ratio_A": [ra[c] for c in parent], "ratio_B": [1.0] * 13})
    cp = {"F1": {"clone3": 0.36, "clone4": 0.09, "clone9": 0.2, "clone11": 0.35},
          "R1": {"clone1": 0.07, "clone5": 0.52, "clone6": 0.07, "clone8": 0.15, "clone13": 0.08, "clone15": 0.11},
          "R2": {"clone2": 0.15, "clone3": 0.1937, "clone4": 0.18, "clone5": 0.2063, "clone7": 0.19, "clone11": 0.03,
                 "clone12": 0.05}}
    cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s), "proportion": cp[s].get(c, 0.0)}
                         for s in cp for c in parent])
    ev = events((True, False, None, ["F1", "R1", "R2"]), (False, True, "3", ["F1", "R2"]), (False, False, None, ["R2"]))
    tu = tumour("Clonal", 2, n_second=1, n_clonal=1)
    reg = regions_df(T, {"F1": 2, "R1": 1, "R2": 3})
    res = placement.place_events(ev, tu, mut_counts(c3=8), rt, cpdf, ploidy_rule="fill", regions=reg)
    tb = final_table.build_table(res.events, res.tumours, tu, rt)
    assert {c: v for c, v in zip(tb["clone"], tb["number_of_gds_at_clone"]) if v} == {"clone1": 1, "clone3": 1, "clone2": 1}
    row = tb.filter(pl.col("clone") == "clone2").row(0, named=True)
    assert "minority_coverage(0.34:R2)" in row["worst_event_quality_code"]
    assert row["worst_event_quality"] < 0.45


def test_event_table_ids_and_rejection():
    """Per-event table with unique ids; rejecting the weak GD recomputes counts and scores."""
    reg = regions_df(T, {"R1": 1, "R2": 1, "R3": 0})
    summ = summary_df(T, {"R1": (0.95, 0.1, 0.03), "R2": (0.9, 0.1, 0.03), "R3": (0.1, 0.05, 0.0)})
    ev = events((False, True, "4", ["R1"]), (False, True, "5", ["R2"]))
    tu = tumour(n_sub=2)
    rt = ratios(c4=2.0, c5=2.0)
    res = run(ev, tu, mc=mut_counts(c4=40, c5=2), ploidy_rule="fill", regions=reg, summary=summ, rt=rt)
    et = final_table.build_events(res.events, res.tumours, tu, rt)
    assert et["wgd_event_id"].to_list() == [f"{T}:clone4:1", f"{T}:clone5:1"]
    assert et["wgd_event_id"].n_unique() == et.height
    assert et["source"].to_list() == ["mutation", "mutation"] and not et["rejected"].any()
    full = final_table.build_table(res.events, res.tumours, tu, rt)
    assert full["total_gds_per_tumour"][0] == 2 and full["quality"][0] == pytest.approx(0.6)  # 2 muts: 0.5+0.5*0.2
    # reject GDs with event_quality < 0.7: clone5 goes, the tumour score recovers to 1
    et_r = final_table.build_events(res.events, res.tumours, tu, rt, reject_below=0.7)
    assert et_r.filter("rejected")["wgd_event_id"].to_list() == [f"{T}:clone5:1"]
    filt = final_table.build_table(res.events, res.tumours, tu, rt, reject_below=0.7, event_table=et_r)
    n = {c: v for c, v in zip(filt["clone"], filt["number_of_gds_at_clone"]) if v}
    assert n == {"clone4": 1}
    row = filt.row(0, named=True)
    assert row["total_gds_per_tumour"] == 1 and row["n_rejected_events"] == 1
    assert row["quality"] == 1.0 and row["quality_code"] is None
    assert row["sub_class"] == "single subclonal GD"
