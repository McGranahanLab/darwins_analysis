"""WGD_WGS_pipeline: config handling, staging, and (if pixi + the PGDD R env exist) a full run on a
synthetic WGS store. Run with the WGD environment:
  pixi run --manifest-path ../WGD/pixi.toml pytest tests/ -q
"""

from __future__ import annotations

import os
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
import wgs  # noqa: E402
from wgd import params as P  # noqa: E402


def config(tmp_path: Path, **override) -> Path:
    text = (HERE / "config.toml").read_text()
    text = text.replace('out_dir = "output/run"', f'out_dir = "{tmp_path / "run"}"')
    text = text.replace('dir = "../../tx_data/data"', f'dir = "{tmp_path / "tx"}"')
    for k, v in override.items():
        text = text.replace(k, v)
    p = tmp_path / "config.toml"
    p.write_text(text)
    return p


def test_wgd_params_match_wgd_defaults(tmp_path):
    """[wgd.params] in config.toml is WGD/params.toml; the written file loads to the same Params."""
    cfg = wgs.load(HERE / "config.toml")
    wgs.write(cfg, tmp_path)
    assert P.load(tmp_path / "config" / "wgd_params.toml") == P.load(HERE.parent / "WGD" / "params.toml")
    s = (tmp_path / "config" / "pgdd_settings.json").read_text()
    assert '"check_mut_cpn_2_threshold": 1.25' in s and '"n_boot": 500' in s


def test_check_lists_missing_tables(tmp_path):
    make_fixture.build(tmp_path / "tx")
    (tmp_path / "tx" / "wgs_alpaca.parquet").unlink()
    probs = wgs.check(wgs.load(config(tmp_path)))
    assert [p for p in probs if "missing" in p] == ["tx_data table missing: wgs_alpaca: steps 04, correction: "
                                                     "ALPACA clone CN per segment"]
    on = config(tmp_path, **{"enabled = false": "enabled = true"})
    assert any("wgs_manual_qc" in p for p in wgs.check(wgs.load(on)))


def test_derived_calls_file_and_case_id(tmp_path):
    make_fixture.build(tmp_path / "tx")
    assert wgs.get(wgs.load(config(tmp_path)), "derived.calls_file") == "wgd_calls_rejected_below_0.3.tsv"
    off = config(tmp_path, **{"quality_cutoff = 0.3": "quality_cutoff = -1"})
    assert wgs.get(wgs.load(off), "derived.calls_file") == "wgd_calls.tsv"
    assert wgs.case_id("LTX0002-Tumour1") == "LTX0002_1"


def test_stage_layout(tmp_path):
    make_fixture.build(tmp_path / "tx")
    res = wgs.stage(wgs.load(config(tmp_path)), tmp_path / "run")
    assert res["tumours"] == 3
    c = tmp_path / "run" / "inputs" / "conipher" / "SYN0001_1"
    assert (c / "SYN0001_1.SCoutput.CLEAN.tsv").exists() and (c / "allTrees.txt").read_text().startswith("# tree 1\n1\t2")
    tt = pl.read_csv(c / "treeTable.tsv", separator="\t")
    assert tt.columns == ["mutation_id", "SAMPLE", "treeCLUSTER"]
    seg = pl.read_csv(tmp_path / "run" / "inputs" / "msai" / "SYN0001_1" /
                      "SYN0001_1_copynumber_MSAI_adjusted_segments_tree_input.tsv", separator="\t")
    assert {"SAMPLE", "STARTPOS", "ENDPOS", "COPY_NUMBER_A", "COPY_NUMBER_B"} <= set(seg.columns)
    cp = pl.read_parquet(tmp_path / "run" / "inputs" / "tx_data" / "clone_proportions.parquet")
    assert cp.columns == ["tumour_id", "clone", "sample_id", "proportion"]


def _tools() -> bool:
    pixi = os.environ.get("PIXI", "pixi")
    return shutil.which(pixi) is not None and shutil.which("micromamba") is not None


@pytest.mark.skipif(not _tools(), reason="needs pixi and micromamba (pgdd env) on PATH")
def test_full_run_on_synthetic_store(tmp_path):
    make_fixture.build(tmp_path / "tx")
    cfg = config(tmp_path, **{'rscript = "Rscript"': 'rscript = "micromamba run -n pgdd Rscript"',
                              "n_boot = 500": "n_boot = 50"})
    r = subprocess.run([str(HERE / "run.sh"), str(cfg)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    run = tmp_path / "run"
    calls = pl.read_csv(run / "wgd" / "05_table" / "wgd_calls_rejected_below_0.3.tsv", separator="\t",
                        null_values=["NA"])
    gds = {(t, c): n for t, c, n in calls.select("tumour_id", "clone", "number_of_gds_at_clone").iter_rows() if n}
    assert gds == {("SYN0001-Tumour1", "clone1"): 1, ("SYN0002-Tumour1", "clone2"): 1}
    assert (run / "wgd" / "07_report" / "report.html").exists()
    assert (run / "correction" / "alpaca_corrected.parquet").exists()
    assert tomllib.loads((run / "logs" / "config.used.toml").read_text())["wgd"]["quality_cutoff"] == 0.3
