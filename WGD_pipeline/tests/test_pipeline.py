"""WGD_pipeline: config handling, staging, and (if pixi + the PGDD R env exist) a full run on a
synthetic store, for both configurations (config.wgs.toml, config.wes.toml). Run with the WGD environment:
  pixi run --manifest-path ../WGD/pixi.toml pytest tests/ -q
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import polars as pl
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "pipeline"))
sys.path.insert(0, str(HERE / "tests"))
sys.path.insert(0, str(HERE.parent / "WGD"))

import make_fixture  # noqa: E402
import prep  # noqa: E402
from wgd import params as P  # noqa: E402
from wgd.pgdd import EXPECTED_SETTINGS  # noqa: E402


PGDD_SRC = os.environ.get("PGDD_SRC", str(Path("~/scratch/PGDD_amf71").expanduser()))
# The synthetic store has 30 mutations per cluster, far below the WGS thresholds in config.wgs.toml
# (Piotr, 2026-11-19: 100 / 50 / 30), so the tests use the WES-scale ones for both cohorts; the pipeline is the same.
FIXTURE_PGDD = {"discover_num_muts_threshold": 10, "discover_num_2_cpn_muts_threshold": 5,
                "check_num_2_cpn_muts_threshold": 3}
COHORTS = ["wgs", "wes"]
# What may differ between config.wgs.toml and config.wes.toml: everything else must be equal.
COHORT_KEYS = {("run", "cohort"), ("run", "out_dir"), ("tx_data", "tables")} | {
    ("pgdd", "settings", k) for k in FIXTURE_PGDD}


def tables(cohort: str) -> dict[str, str]:
    return tomllib.loads((HERE / f"config.{cohort}.toml").read_text())["tx_data"]["tables"]


def fixture(tmp_path: Path, cohort: str = "wgs") -> Path:
    return make_fixture.build(tmp_path / "tx", tables(cohort))


def config(tmp_path: Path, cohort: str = "wgs", **override) -> Path:
    text = (HERE / f"config.{cohort}.toml").read_text()
    text = re.sub(r'(?m)^src = ".*?"', f'src = "{PGDD_SRC}"', text)
    for k, v in FIXTURE_PGDD.items():
        text = re.sub(rf"(?m)^{k} = \S+", f"{k} = {v}", text)
    text = re.sub(r'(?m)^out_dir = ".*?"', f'out_dir = "{tmp_path / "run"}"', text)
    text = text.replace('dir = "../../tx_data/data"', f'dir = "{tmp_path / "tx"}"')
    text = re.sub(r'manual_qc_table = ".*"', f'manual_qc_table = "{tmp_path / "tx" / "manual_qc.tsv"}"', text)
    for k, v in override.items():
        text = text.replace(k, v)
    p = tmp_path / "config.toml"
    p.write_text(text)
    return p


def _drop(d: dict, keys: set, prefix: tuple = ()) -> dict:
    return {k: (_drop(v, keys, prefix + (k,)) if isinstance(v, dict) else v)
            for k, v in d.items() if prefix + (k,) not in keys}


def test_wes_and_wgs_configs_differ_only_in_cohort_settings():
    """One pipeline, two configs: they may differ only in cohort, output dir, tables and PGDD counts."""
    wgs_cfg, wes_cfg = (tomllib.loads((HERE / f"config.{c}.toml").read_text()) for c in COHORTS)
    assert _drop(wgs_cfg, COHORT_KEYS) == _drop(wes_cfg, COHORT_KEYS)
    assert (wgs_cfg["run"]["cohort"], wes_cfg["run"]["cohort"]) == ("wgs", "wes")
    assert wgs_cfg["run"]["out_dir"] != wes_cfg["run"]["out_dir"]
    # WES: Piotr's 2026-11-09 settings (wgd.pgdd.EXPECTED_SETTINGS); WGS: 100 / 50 / 30 (2026-11-19)
    assert {k: wes_cfg["pgdd"]["settings"][k] for k in EXPECTED_SETTINGS} == EXPECTED_SETTINGS
    assert [wgs_cfg["pgdd"]["settings"][k] for k in FIXTURE_PGDD] == [100, 50, 30]
    assert set(tables("wgs")) == set(tables("wes")) == set(prep.INPUTS) | set(prep.QC_INPUTS)


@pytest.mark.parametrize("cohort", COHORTS)
def test_check_rejects_tables_of_the_other_cohort(tmp_path, cohort):
    fixture(tmp_path, cohort)
    assert not [p for p in prep.check(prep.load(config(tmp_path, cohort))) if "cohort" in p]
    other = "wes" if cohort == "wgs" else "wgs"
    flipped = config(tmp_path, cohort, **{f'cohort = "{cohort}"': f'cohort = "{other}"'})
    assert any("tables of the other cohort" in p for p in prep.check(prep.load(flipped)))
    t = tables(cohort)["alpaca"]
    mixed = config(tmp_path, cohort, **{f'alpaca = "{t}"': 'alpaca = "wgs_alpaca"' if cohort == "wes"
                                        else 'alpaca = "alpaca"'})
    assert any("tables of the other cohort" in p for p in prep.check(prep.load(mixed)))


def test_old_pgdd_checkout_is_caught(tmp_path):
    fixture(tmp_path)
    old = tmp_path / "pgdd_old" / "R"
    old.mkdir(parents=True)
    (old / "parrallel_gd_caller.R").write_text("detect_par_gd <- function(input, mut_cpn_2_threshold = 1.5) {}\n")
    cfg = config(tmp_path)
    cfg.write_text(re.sub(r'(?m)^src = ".*?"', f'src = "{old.parent}"', cfg.read_text()))
    assert any("does not know settings" in p for p in prep.check(prep.load(cfg)))


@pytest.mark.parametrize("cohort", COHORTS)
def test_wgd_params_match_wgd_defaults(tmp_path, cohort):
    """[wgd.params] in each config is WGD/params.toml; the written file loads to the same Params."""
    cfg = prep.load(HERE / f"config.{cohort}.toml")
    prep.write(cfg, tmp_path)
    assert P.load(tmp_path / "config" / "wgd_params.toml") == P.load(HERE.parent / "WGD" / "params.toml")
    s = (tmp_path / "config" / "pgdd_settings.json").read_text()
    assert '"check_mut_cpn_2_threshold": 1.25' in s and '"n_boot": 500' in s


@pytest.mark.parametrize("cohort", COHORTS)
def test_check_lists_missing_tables(tmp_path, cohort):
    fixture(tmp_path, cohort)
    t = tables(cohort)
    (tmp_path / "tx" / f"{t['alpaca']}.parquet").unlink()
    probs = prep.check(prep.load(config(tmp_path, cohort)))
    assert [p for p in probs if "missing" in p] == [f"tx_data table missing: {t['alpaca']}: steps 04, correction: "
                                                     "ALPACA clone CN per segment"]
    (tmp_path / "tx" / f"{t['conipher_cluster_info']}.parquet").unlink()
    (tmp_path / "tx" / "manual_qc.tsv").unlink()
    probs = prep.check(prep.load(config(tmp_path, cohort)))
    assert any(t["conipher_cluster_info"] in p for p in probs) and any("manual_qc_table" in p for p in probs)
    off = config(tmp_path, cohort, **{"enabled = true ": "enabled = false"})
    assert not any("cluster_info" in p or "manual_qc" in p for p in prep.check(prep.load(off)))


def test_derived_calls_file_and_case_id(tmp_path):
    fixture(tmp_path)
    assert prep.get(prep.load(config(tmp_path)), "derived.calls_file") == "wgd_calls_rejected_below_0.3.tsv"
    off = config(tmp_path, **{"quality_cutoff = 0.3": "quality_cutoff = -1"})
    assert prep.get(prep.load(off), "derived.calls_file") == "wgd_calls.tsv"
    assert prep.case_id("LTX0002-Tumour1") == "LTX0002_1"


@pytest.mark.parametrize("cohort", COHORTS)
def test_stage_layout(tmp_path, cohort):
    fixture(tmp_path, cohort)
    res = prep.stage(prep.load(config(tmp_path, cohort)), tmp_path / "run")
    assert res["tumours"] == 3 and res["cohort"] == cohort
    c = tmp_path / "run" / "inputs" / "conipher" / "SYN0001_1"
    assert (c / "SYN0001_1.SCoutput.CLEAN.tsv").exists() and (c / "allTrees.txt").read_text().startswith("# tree 1\n1\t2")
    tt = pl.read_csv(c / "treeTable.tsv", separator="\t")
    assert tt.columns == ["mutation_id", "SAMPLE", "treeCLUSTER"]
    sc = pl.read_csv(c / "SYN0001_1.SCoutput.CLEAN.tsv", separator="\t")
    assert {"CASE_ID", "SAMPLE", "CHR", "POS", "CLUSTER", "MUT_COPY", "COPY_NUMBER_A", "mutation_id"} <= set(sc.columns)
    ci = tmp_path / "run" / "inputs" / "qc" / "release" / "SYN0001" / "mutation_trees" / "tumour_1" / "clusterInfo.txt"
    assert pl.read_csv(ci, separator="\t").columns == ["clusterID", "SAMPLE", "truncal", "treeClust", "meanCCF"]
    assert (tmp_path / "run" / "inputs" / "qc" / "manual_qc.tsv").exists()
    seg = pl.read_csv(tmp_path / "run" / "inputs" / "msai" / "SYN0001_1" /
                      "SYN0001_1_copynumber_MSAI_adjusted_segments_tree_input.tsv", separator="\t")
    assert {"SAMPLE", "STARTPOS", "ENDPOS", "COPY_NUMBER_A", "COPY_NUMBER_B"} <= set(seg.columns)
    cp = pl.read_parquet(tmp_path / "run" / "inputs" / "tx_data" / "clone_proportions.parquet")
    assert cp.columns == ["tumour_id", "clone", "sample_id", "proportion"]


def _tools() -> bool:
    pixi = os.environ.get("PIXI", "pixi")
    return shutil.which(pixi) is not None and shutil.which("micromamba") is not None


@pytest.mark.skipif(not _tools(), reason="needs pixi and micromamba (pgdd env) on PATH")
@pytest.mark.parametrize("cohort", COHORTS)
def test_full_run_on_synthetic_store(tmp_path, cohort):
    fixture(tmp_path, cohort)
    cfg = config(tmp_path, cohort, **{'rscript = "Rscript"': 'rscript = "micromamba run -n pgdd Rscript"',
                              "n_boot = 500": "n_boot = 50"})
    r = subprocess.run([str(HERE / "run.sh"), str(cfg)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    run = tmp_path / "run"
    calls = pl.read_csv(run / "wgd" / "05_table" / "wgd_calls_rejected_below_0.3.tsv", separator="\t",
                        null_values=["NA"])
    gds = {(t, c): n for t, c, n in calls.select("tumour_id", "clone", "number_of_gds_at_clone").iter_rows() if n}
    # SYN0003 is single-clone and doubled: placed by ploidy on its only clone (2026-11-17)
    assert gds == {("SYN0001-Tumour1", "clone1"): 1, ("SYN0002-Tumour1", "clone2"): 1, ("SYN0003-Tumour1", "clone1"): 1}
    assert (run / "wgd" / "07_report" / "report.html").exists()
    assert (run / "correction" / "alpaca_corrected.parquet").exists()
    used = tomllib.loads((run / "logs" / "config.used.toml").read_text())
    assert used["wgd"]["quality_cutoff"] == 0.3 and used["run"]["cohort"] == cohort
    assert f"TRACERx {cohort.upper()}" in (run / "wgd" / "07_report" / "report.html").read_text()


def test_stage_stops_when_no_tumours(tmp_path):
    """2026-11-18 (Piotr's run): an empty wgs_conipher_mutations must stop `stage`, not step 00."""
    fixture(tmp_path)
    m = pl.read_parquet(tmp_path / "tx" / "wgs_conipher_mutations.parquet")
    m.clear().write_parquet(tmp_path / "tx" / "wgs_conipher_mutations.parquet")
    with pytest.raises(SystemExit, match="no tumours to stage"):
        prep.stage(prep.load(config(tmp_path)), tmp_path / "run")


@pytest.mark.parametrize("cohort", COHORTS)
def test_stage_skips_tumour_without_segments(tmp_path, cohort):
    """2026-11-26 (Piotr): a tumour without an input is left out and recorded; "stop" keeps the old behaviour."""
    fixture(tmp_path, cohort)
    seg = tmp_path / "tx" / f"{tables(cohort)['segments']}.parquet"
    pl.read_parquet(seg).filter(pl.col("tumour_id") != "SYN0002-Tumour1").write_parquet(seg)
    res = prep.stage(prep.load(config(tmp_path, cohort)), tmp_path / "run")
    assert res["tumours"] == 2 and res["failed_tumours"] == 1
    f = pl.read_csv(tmp_path / "run" / "inputs" / "failed_tumours.tsv", separator="\t")
    assert f.rows() == [("SYN0002-Tumour1", "stage", f"no CN segments ({tables(cohort)['segments']})")]
    stop = config(tmp_path, cohort, **{'on_tumour_failure = "skip"': 'on_tumour_failure = "stop"'})
    with pytest.raises(SystemExit, match="SYN0002-Tumour1: no CN segments"):
        prep.stage(prep.load(stop), tmp_path / "run2")


@pytest.mark.skipif(not _tools(), reason="needs pixi and micromamba (pgdd env) on PATH")
def test_full_run_completes_with_failing_tumours(tmp_path):
    """One tumour fails at stage (no segments), one at step 00 (SAMPLE ids differ); the run completes and
    failed_tumours.tsv lists both with step and reason."""
    fixture(tmp_path)
    seg = tmp_path / "tx" / "wgs_consensus_segments.parquet"
    s = pl.read_parquet(seg).filter(pl.col("tumour_id") != "SYN0003-Tumour1")
    s.with_columns(pl.when(pl.col("tumour_id") == "SYN0002-Tumour1").then(pl.col("sample") + "_x")
                   .otherwise(pl.col("sample")).alias("sample")).write_parquet(seg)
    cfg = config(tmp_path, **{'rscript = "Rscript"': 'rscript = "micromamba run -n pgdd Rscript"',
                              "n_boot = 500": "n_boot = 50"})
    r = subprocess.run([str(HERE / "run.sh"), str(cfg)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    run = tmp_path / "run"
    f = pl.read_csv(run / "failed_tumours.tsv", separator="\t")
    assert f.select("tumour_id", "step").rows() == [("SYN0002-Tumour1", "wgd00"), ("SYN0003-Tumour1", "stage")]
    assert "SAMPLE ids differ" in f["reason"][0]
    calls = pl.read_csv(run / "wgd" / "05_table" / "wgd_calls.tsv", separator="\t", null_values=["NA"])
    assert set(calls["tumour_id"]) == {"SYN0001-Tumour1"}
    assert "2 tumour(s) left out of the run" in r.stdout
