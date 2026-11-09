"""Step 00 (PGDD input) and step 03 (cohort summary) on small synthetic CONIPHER/MSAI inputs.

`test_parity_with_r` runs Piotr's R scripts on the same inputs and compares the
outputs. It needs Rscript with data.table: set WGD_RSCRIPT (e.g.
"micromamba run -n pgdd Rscript"), otherwise it is skipped.
"""

from __future__ import annotations

import math
import os
import shlex
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from wgd import pgdd_input, summary

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
TREE = "### 2 trees\n# tree 1\n1\t2\n1\t3\n# tree 2\n1\t2\n2\t3\n"


def _sample(case: str, r: int) -> str:
    return f"{case.split('_')[0]}_SU_T1-R{r}--{case[-1]}{r}abc"


def _write_tumour(conipher: Path, msai: Path, case: str, *, tree: str = TREE,
                  cn_samples: list[str] | None = None, seed: int = 0) -> None:
    """Clusters 1 (clonal), 2, 3 on the tree; 5 not in treeTable; 6 in treeTable but not an edge."""
    rng = np.random.default_rng(seed)
    samples = [_sample(case, 1), _sample(case, 2)]
    rows, tree_rows = [], []
    pos = 1000
    for cluster, n in ((1, 6), (2, 4), (3, 4), (5, 3), (6, 3)):
        for _ in range(n):
            pos += 137
            mid = f":{(pos % 3) + 1}:{pos}:C:T"
            for s in samples:
                a, b = (2, 2) if s.endswith("1abc") else (1, 1)
                rows.append(dict(CASE_ID=case, SAMPLE=s, CHR=(pos % 3) + 1, POS=pos, REF="C", ALT="T",
                                 CLUSTER=cluster, MUT_COPY=round(float(rng.uniform(0.5, 2.2)), 4),
                                 COPY_NUMBER_A=a, COPY_NUMBER_B=b, PLOIDY=3.1 if a == 2 else 2.05,
                                 mutation_id=mid))
                if cluster != 5:
                    tree_rows.append(dict(SAMPLE=s, mutation_id=mid, treeCLUSTER=cluster))
    t = conipher / case
    t.mkdir(parents=True)
    pl.DataFrame(rows).write_csv(t / f"{case}.SCoutput.CLEAN.tsv", separator="\t")
    pl.DataFrame(tree_rows).write_csv(t / "treeTable.tsv", separator="\t")
    (t / "allTrees.txt").write_text(tree)

    segs = []
    for s in cn_samples or samples:
        doubled = s.endswith("1abc")
        # (start, end, A, B): 60% at 2|2 (or 1|1), then a 2.5 (rounds to 2) and a 3.5 (rounds to 4) segment, one NA
        for start, end, a, b in ((1, 600, 2, 2) if doubled else (1, 600, 1, 1),
                                 (601, 800, 2.5, 1), (801, 900, 3.5, 1), (901, 1000, None, None)):
            segs.append(dict(SAMPLE=s, chr=1, STARTPOS=start * 1000, ENDPOS=end * 1000,
                             COPY_NUMBER_A=a, COPY_NUMBER_B=b))
    m = msai / case
    m.mkdir(parents=True)
    pl.DataFrame(segs).write_csv(m / f"{case}_copynumber_MSAI_adjusted_segments_tree_input.tsv",
                                 separator="\t", null_value="NA")


@pytest.fixture
def inputs(tmp_path):
    conipher, msai = tmp_path / "conipher", tmp_path / "msai"
    _write_tumour(conipher, msai, "LTX0901_1", seed=1)
    _write_tumour(conipher, msai, "LTX0902_1", seed=2)
    _write_tumour(conipher, msai, "LTX0903_1", tree="### 1 trees\n# tree 1\n1\t1\n")  # single clone (CONIPHER format)
    return conipher, msai


def test_prepare(inputs):
    res = pgdd_input.prepare_all(*inputs, log=lambda *_: None)
    t = res.table
    assert t.columns == list(pgdd_input.OUTPUT_COLS)
    assert sorted(t["tumour_id"].unique()) == ["LTX0901-Tumour1", "LTX0902-Tumour1"]
    assert res.skipped == {}
    assert res.single_clone["tumour_id"].unique().to_list() == ["LTX0903-Tumour1"]
    # cluster 5 is not in treeTable; cluster 6 is not an edge of the tree
    assert set(t["cluster_id"].unique()) == {1, 2, 3}
    assert t.filter(pl.col("is_clonal_cluster"))["cluster_id"].unique().to_list() == [1]
    assert t["chromosome"].str.starts_with("chr").all()
    assert t.height == 2 * 2 * (6 + 4 + 4)


def test_num_gds_and_fractions(inputs):
    t = pgdd_input.prepare_all(*inputs, log=lambda *_: None).table
    per = t.select("sample_id", "num_gds", "frac_maj_ge2", "frac_maj_ge3", "frac_maj_ge4", "ploidy").unique()
    r1 = per.filter(pl.col("sample_id") == _sample("LTX0901_1", 1)).row(0, named=True)
    r2 = per.filter(pl.col("sample_id") == _sample("LTX0901_1", 2)).row(0, named=True)
    # R1: 600 kb at 2|2 + 200 kb at 2.5 (rounds to 2) + 100 kb at 3.5 (rounds to 4); NA segment only in the total
    assert r1["frac_maj_ge2"] == pytest.approx(900 / 1000, abs=1e-3)
    assert r1["frac_maj_ge3"] == pytest.approx(100 / 1000, abs=1e-3)
    assert r1["frac_maj_ge4"] == pytest.approx(100 / 1000, abs=1e-3)
    assert r1["num_gds"] == 1 and r1["ploidy"] == 3.1
    # R2: 1|1 except the 2.5 and 3.5 segments -> 30% doubled
    assert r2["frac_maj_ge2"] == pytest.approx(300 / 1000, abs=1e-3)
    assert r2["num_gds"] == 0


def test_level_confidence():
    a, b, w = np.array([2.0, 1.0]), np.array([2.0, 1.0]), np.array([1.0, 1.0])
    assert pgdd_input._level_confidence(a, b, w, 2.0) == pytest.approx(0.5)
    assert pgdd_input._level_confidence(np.array([5.0]), np.array([4.0]), np.array([1.0]), 4.0) == 1.0
    assert math.isnan(pgdd_input._level_confidence(np.array([np.nan]), np.array([1.0]), np.array([1.0]), 2.0))


def test_default_tree_only(tmp_path):
    f = tmp_path / "allTrees.txt"
    f.write_text("### 2 trees\n# tree 1\n1\t2\n1\t3\n# tree 2\n1\t2\n2\t4\n")
    edges = pgdd_input.default_tree_edges(f)
    assert pgdd_input.tree_clusters(edges) == {1, 2, 3}
    assert pgdd_input.clonal_cluster(edges) == 1


def test_single_clone(inputs):
    res = pgdd_input.prepare_all(*inputs, log=lambda *_: None)
    sc = res.single_clone
    assert sc.columns == list(pgdd_input.SINGLE_CLONE_SCHEMA)
    assert sc["clone"].unique().to_list() == ["clone1"]
    assert dict(zip(sc["sample_id"], sc["num_gds"])) == {_sample("LTX0903_1", 1): 1, _sample("LTX0903_1", 2): 0}
    ev = pgdd_input.single_clone_events(sc).row(0, named=True)
    assert ev["tumour_id"] == "LTX0903-Tumour1" and ev["clone"] == "clone1"
    assert ev["n_clonal_gds"] == 1 and ev["num_gds_discordant"]  # synthetic R2 is not doubled


def test_one_edge_two_clone_tree_goes_to_pgdd(tmp_path):
    conipher, msai = tmp_path / "conipher", tmp_path / "msai"
    _write_tumour(conipher, msai, "LTX0905_1", tree="### 1 trees\n# tree 1\n1\t2\n")
    res = pgdd_input.prepare_all(conipher, msai, log=lambda *_: None)
    assert res.single_clone.is_empty()
    assert set(res.table["cluster_id"].unique()) == {1, 2}


def test_sample_mismatch_fails_before_processing(inputs, monkeypatch):
    conipher, msai = inputs
    _write_tumour(conipher, msai, "LTX0904_1", cn_samples=[_sample("LTX0904_1", 1)])  # CN missing R2
    _write_tumour(conipher, msai, "LTX0906_1")
    (conipher / "LTX0906_1" / "treeTable.tsv").write_text(
        "SAMPLE\tmutation_id\ttreeCLUSTER\nLTX0906_SU_T1-R9--x\t:1:1:C:T\t1\n")
    called = []
    monkeypatch.setattr(pgdd_input, "prepare_tumour", lambda *a, **k: called.append(a))
    with pytest.raises(pgdd_input.SampleMismatchError) as e:
        pgdd_input.prepare_all(conipher, msai, log=lambda *_: None)
    assert "LTX0904_1" in str(e.value) and "LTX0906_1" in str(e.value) and "2 tumour(s)" in str(e.value)
    assert called == []


def test_summary(inputs, tmp_path):
    t = pgdd_input.prepare_all(*inputs, log=lambda *_: None).table
    samples = t.select("sample_id", "tumour_id").unique().sort("sample_id")
    regions = samples.head(3).with_columns(pl.Series("num_gds", [2, 1, 0]))
    s = summary.cohort_summary(t, regions)
    assert s.columns == list(summary.COLUMNS)
    assert s.height == samples.height
    assert s["source"].unique().to_list() == ["pgdd"]
    assert s["tumour_id"].null_count() == 0  # 4th sample is not in regions: tumour_id from the input
    got = dict(zip(s["sample_name_hash"], s["Second GD predicted"]))
    assert [got[x] for x in regions["sample_id"]] == [1, 0, 0]
    r1 = s.filter(pl.col("sample_name_hash") == _sample("LTX0901_1", 1)).row(0, named=True)
    assert r1["genome_frac_nMaj2"] == pytest.approx(0.9, abs=1e-3) and r1["Ploidy"] == 3.1


def test_summary_single_clone(inputs):
    res = pgdd_input.prepare_all(*inputs, log=lambda *_: None)
    regions = res.table.select("sample_id", "tumour_id").unique().with_columns(pl.lit(0).alias("num_gds"))
    s = summary.cohort_summary(res.table, regions, res.single_clone)
    sc = s.filter(pl.col("source") == "single_clone")
    assert sc["tumour_id"].unique().to_list() == ["LTX0903-Tumour1"]
    assert dict(zip(sc["sample_name_hash"], sc["GD predicted"])) == {_sample("LTX0903_1", 1): 1,
                                                                     _sample("LTX0903_1", 2): 0}
    assert sc["genome_frac_nMaj2"].null_count() == 0


def _assert_same(r: pl.DataFrame, py: pl.DataFrame) -> None:
    assert r.columns == py.columns
    assert r.height == py.height
    for c in r.columns:
        if r.schema[c].is_numeric() and py.schema[c].is_numeric():
            assert np.allclose(r[c].cast(pl.Float64).to_numpy(), py[c].cast(pl.Float64).to_numpy(), equal_nan=True), c
        else:
            assert r[c].cast(pl.Utf8).to_list() == py[c].cast(pl.Utf8).to_list(), c


def _rscript() -> list[str] | None:
    cmd = os.environ.get("WGD_RSCRIPT")
    return shlex.split(cmd) if cmd else None


@pytest.mark.skipif(_rscript() is None, reason="set WGD_RSCRIPT to compare with the R scripts")
def test_parity_with_r(inputs, tmp_path):
    conipher, msai = inputs
    r_out, py_out = tmp_path / "r", tmp_path / "py"
    subprocess.run([*_rscript(), str(SCRIPTS / "00_prepare_pgdd_input.R"), str(conipher), str(msai), str(r_out)],
                   check=True, capture_output=True)
    res = pgdd_input.prepare_all(conipher, msai, log=lambda *_: None)
    py_out.mkdir()
    pgdd_input.write_tsv(res.table, py_out / "pgdd_input_all_tumours.tsv")

    r = summary.read_tsv(r_out / "pgdd_input_all_tumours.tsv")
    py = summary.read_tsv(py_out / "pgdd_input_all_tumours.tsv")
    _assert_same(r, py)

    regions = py.select("sample_id", "tumour_id").unique().sort("sample_id").with_columns(
        pl.Series("num_gds", [2, 1, 0, 1]))
    regions.write_csv(tmp_path / "regions.tsv", separator="\t")
    subprocess.run([*_rscript(), str(SCRIPTS / "03_generate_summary.R"), str(r_out / "pgdd_input_all_tumours.tsv"),
                    str(tmp_path / "regions.tsv"), str(r_out)], check=True, capture_output=True)
    rs = summary.read_tsv(r_out / "pgdd_cohort_gd_summary.tsv").sort("sample_name_hash")
    ps = summary.cohort_summary(py, regions).sort("sample_name_hash").drop("source")  # R has no source column
    assert rs.columns == ps.columns
    _assert_same(rs, ps)
