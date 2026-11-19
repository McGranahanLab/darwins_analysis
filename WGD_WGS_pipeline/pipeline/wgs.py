"""Helper for run.sh: read config.toml and stage tx_data WGS tables into the layouts the steps read.

  python wgs.py check  <config>            validate the config; list missing tx_data tables (exit 1 if any)
  python wgs.py get    <config> <key>      print one value (dotted key; paths resolved; lists space-joined)
  python wgs.py write  <config> <run_dir>  write pgdd_settings.json and wgd_params.toml into <run_dir>/config/
  python wgs.py stage  <config> <run_dir>  tx_data (+ the manual QC file) -> <run_dir>/inputs/{conipher,msai,tx_data,qc}

Staged layout (what each step expects):
  inputs/conipher/<CASE>/<CASE>.SCoutput.CLEAN.tsv, treeTable.tsv, allTrees.txt   step 00
  inputs/msai/<CASE>/<CASE>_copynumber_MSAI_adjusted_segments_tree_input.tsv      step 00
  inputs/tx_data/{clone_proportions,alpaca,alpaca_wgd_ratios}.parquet              steps 02, 04, 05, correction
  inputs/qc/manual_qc.tsv, inputs/qc/release/<patient>/mutation_trees/tumour_<n>/clusterInfo.txt   step 04 QC
CASE is the CONIPHER case id (`LTX0001_1` for `LTX0001-Tumour1`).
"""

from __future__ import annotations

import json
import shutil
import sys
import tomllib
from pathlib import Path

import polars as pl

# tx_data WGS tables the pipeline reads (tx_data PR #11, 2026-11-17).
TABLES = {
    "wgs_consensus_segments": ("step 00: MSAI-adjusted CN segments per sample", True),
    "wgs_clone_proportions": ("steps 02, 04, correction: clone composition of each sample", True),
    "wgs_conipher_mutations": ("step 00: CONIPHER SCoutput.CLEAN + treeTable per mutation", True),
    "wgs_conipher_tree_edges": ("step 00: CONIPHER default tree (allTrees.txt)", True),
    "wgs_alpaca": ("steps 04, correction: ALPACA clone CN per segment", True),
    "wgs_alpaca_wgd_ratios": ("steps 02, 04, 05, correction: ALPACA tree + WGD ratios", True),
}
QC_TABLES = {
    "wgs_conipher_cluster_info": "CCF QC: CONIPHER clusterInfo (truncal, treeClust, meanCCF)",
}
# Columns each staged file must carry (the step 00 readers check the rest).
MUT_OUT = {  # SCoutput.CLEAN column step 00 reads <- tx_data wgs_conipher_mutations column
    "SAMPLE": "sample", "CHR": "chr", "POS": "pos", "REF": "ref", "ALT": "alt", "CLUSTER": "cluster",
    "MUT_COPY": "mut_copy", "COPY_NUMBER_A": "cn_a", "COPY_NUMBER_B": "cn_b", "mutation_id": "mutation_id",
    "PLOIDY": "ploidy"}


def load(config: Path) -> dict:
    cfg = tomllib.loads(Path(config).read_text())
    cfg["_base"] = Path(config).resolve().parent
    return cfg


def path(cfg: dict, value: str) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else (cfg["_base"] / p).resolve()


def get(cfg: dict, key: str):
    if key == "derived.calls_file":      # step 05 table the correction reads
        q = float(cfg["wgd"]["quality_cutoff"])
        return "wgd_calls.tsv" if q < 0 else f"wgd_calls_rejected_below_{q:g}.tsv"
    v = cfg
    for k in key.split("."):
        v = v[k]
    if key in ("run.out_dir", "tx_data.dir", "pgdd.src", "qc.manual_qc_table"):
        return str(path(cfg, v))
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, list):
        return " ".join(map(str, v))
    return v


def case_id(tumour_id: str) -> str:
    """`LTX0001-Tumour1` -> CONIPHER case `LTX0001_1` (inverse of pgdd_input.canonical_tumour_id)."""
    return tumour_id.replace("-Tumour", "_")


def missing_tables(cfg: dict) -> list[str]:
    tx = path(cfg, cfg["tx_data"]["dir"])
    need = dict(TABLES)
    if cfg.get("qc", {}).get("enabled"):
        need |= {k: (v, False) for k, v in QC_TABLES.items()}
    return [f"{t}: {need[t][0]}" for t in need if not (tx / f"{t}.parquet").exists()]


def check(cfg: dict) -> list[str]:
    problems = []
    for key in ("run.out_dir", "tx_data.dir", "pgdd.src", "pgdd.rscript", "wgd.quality_cutoff", "correction.min_quality"):
        try:
            get(cfg, key)
        except KeyError:
            problems.append(f"config: missing {key}")
    if not problems:
        if not path(cfg, cfg["tx_data"]["dir"]).is_dir():
            problems.append(f"tx_data.dir not found: {path(cfg, cfg['tx_data']['dir'])}")
        caller = path(cfg, cfg["pgdd"]["src"]) / "R" / "parrallel_gd_caller.R"
        if not caller.exists():
            problems.append(f"pgdd.src has no R/parrallel_gd_caller.R: {path(cfg, cfg['pgdd']['src'])}")
        else:  # an older PGDD (e.g. the wlippa fork before bootstrap) rejects the settings at run time
            src = caller.read_text()
            old = [k for k in cfg["pgdd"].get("settings", {}) if k not in src]
            if old:
                problems.append(f"pgdd.src {caller.parent.parent} does not know settings {old}: use "
                                "amf71/ParallelGDDetect at cb0437b or later (bootstrap/stability arguments)")
        problems += [f"tx_data table missing: {m}" for m in missing_tables(cfg)]
        if cfg.get("qc", {}).get("enabled"):
            q = cfg["qc"].get("manual_qc_table")
            if not q or not path(cfg, q).is_file():
                problems.append(f"qc.manual_qc_table not found: {q}")
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "WGD"))
    from wgd import params as P
    try:
        P.Params(**{**_flat_params(cfg)})
    except TypeError as e:
        problems.append(f"[wgd.params]: {e}")
    return problems


def _flat_params(cfg: dict) -> dict:
    p = dict(cfg.get("wgd", {}).get("params", {}))
    fixed = p.pop("fixed", None)
    for k in ("ambiguity_majorities", "prior_raised_factor", "prior_lowered_factor", "minor_clone_factor"):
        if k in p:
            p[k] = tuple(p[k])
    if fixed is not None:
        from wgd import params as P
        p["fixed"] = {**P.Params().fixed, **fixed}
    return p


def _toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    return repr(v)


def write(cfg: dict, run_dir: Path) -> None:
    out = run_dir / "config"
    out.mkdir(parents=True, exist_ok=True)
    (out / "pgdd_settings.json").write_text(json.dumps(cfg["pgdd"]["settings"], indent=2) + "\n")
    p = dict(cfg["wgd"].get("params", {}))
    fixed = p.pop("fixed", {})
    lines = [f"{k} = {_toml_value(v)}" for k, v in p.items()]
    if fixed:
        lines += ["", "[fixed]"] + [f"{k} = {_toml_value(v)}" for k, v in fixed.items()]
    (out / "wgd_params.toml").write_text("# written by WGD_WGS_pipeline from [wgd.params]\n" + "\n".join(lines) + "\n")


def _tumours(cfg: dict, df: pl.DataFrame) -> pl.DataFrame:
    sel = cfg["run"].get("tumours") or []
    return df.filter(pl.col("tumour_id").is_in(sel)) if sel else df


def stage(cfg: dict, run_dir: Path) -> dict:
    tx = path(cfg, cfg["tx_data"]["dir"])
    miss = missing_tables(cfg)
    if miss:
        raise SystemExit("tx_data tables missing:\n  " + "\n  ".join(miss))
    read = lambda t: _tumours(cfg, pl.read_parquet(tx / f"{t}.parquet"))  # noqa: E731
    inp = run_dir / "inputs"
    seg, mut, edges = read("wgs_consensus_segments"), read("wgs_conipher_mutations"), read("wgs_conipher_tree_edges")
    tumours = sorted(set(mut["tumour_id"]))
    if not tumours:
        sel = cfg["run"].get("tumours") or []
        counts = {t: pl.read_parquet(tx / f"{t}.parquet").height for t in TABLES}
        raise SystemExit(
            f"no tumours to stage: wgs_conipher_mutations in {tx} has no rows"
            + (f" for [run].tumours {sel}" if sel else "") + ".\n  rows per table: "
            + ", ".join(f"{t}={n}" for t, n in counts.items())
            + "\n  If it is empty, check the tx_data build log for 'wgs_conipher_mutations: read N/M WGS tumours'"
              " and 'no file for ... tumour(s)' (CONIPHER root conipher_wgs in catalog/_sources.yml).")
    no_seg = sorted(set(tumours) - set(seg["tumour_id"]))
    if no_seg:
        raise SystemExit(f"tumours with CONIPHER mutations but no CN segments: {no_seg[:10]}")
    for t in tumours:
        case = case_id(t)
        cdir, mdir = inp / "conipher" / case, inp / "msai" / case
        cdir.mkdir(parents=True, exist_ok=True)
        mdir.mkdir(parents=True, exist_ok=True)
        m = mut.filter(pl.col("tumour_id") == t)
        (m.select(pl.lit(case).alias("CASE_ID"), *[pl.col(src).alias(dst) for dst, src in MUT_OUT.items()])
          .write_csv(cdir / f"{case}.SCoutput.CLEAN.tsv", separator="\t", null_value="NA"))
        (m.filter(pl.col("tree_cluster").is_not_null())
          .select("mutation_id", pl.col("sample").alias("SAMPLE"), pl.col("tree_cluster").alias("treeCLUSTER"))
          .write_csv(cdir / "treeTable.tsv", separator="\t"))
        e = edges.filter(pl.col("tumour_id") == t).sort("parent_cluster", "child_cluster")
        (cdir / "allTrees.txt").write_text("# tree 1\n" + "".join(
            f"{p}\t{c}\n" for p, c in e.select("parent_cluster", "child_cluster").iter_rows()))
        (seg.filter(pl.col("tumour_id") == t)
            .select(pl.col("sample").alias("SAMPLE"), pl.col("chr").alias("CHR"), pl.col("start").alias("STARTPOS"),
                    pl.col("end").alias("ENDPOS"), pl.col("cn_major").alias("COPY_NUMBER_A"),
                    pl.col("cn_minor").alias("COPY_NUMBER_B"))
            .write_csv(mdir / f"{case}_copynumber_MSAI_adjusted_segments_tree_input.tsv", separator="\t"))
    txo = inp / "tx_data"
    txo.mkdir(parents=True, exist_ok=True)
    read("wgs_clone_proportions").select("tumour_id", "clone", "sample_id", "proportion").write_parquet(
        txo / "clone_proportions.parquet")
    read("wgs_alpaca").write_parquet(txo / "alpaca.parquet")
    read("wgs_alpaca_wgd_ratios").write_parquet(txo / "alpaca_wgd_ratios.parquet")
    qc = None
    if cfg.get("qc", {}).get("enabled"):
        qdir = inp / "qc"
        qdir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path(cfg, cfg["qc"]["manual_qc_table"]), qdir / "manual_qc.tsv")
        ci = read("wgs_conipher_cluster_info")
        rb = lambda c: pl.when(pl.col(c)).then(pl.lit("TRUE")).when(~pl.col(c)).then(pl.lit("FALSE"))  # noqa: E731
        for (t,), g in ci.group_by("tumour_id"):
            patient, n = t.split("-Tumour")
            d = qdir / "release" / patient / "mutation_trees" / f"tumour_{n}"
            d.mkdir(parents=True, exist_ok=True)
            (g.select(pl.col("cluster_id").alias("clusterID"), pl.col("sample").alias("SAMPLE"),
                      rb("truncal").alias("truncal"), rb("tree_clust").alias("treeClust"),
                      pl.col("mean_ccf").alias("meanCCF"))
              .write_csv(d / "clusterInfo.txt", separator="\t", null_value="NA"))
        qc = str(qdir)
    summary = {"tumours": len(tumours), "mutations": mut.height, "segments": seg.height, "qc": qc}
    (inp / "staged.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv: list[str]) -> None:
    if len(argv) < 2 or argv[0] not in ("check", "get", "write", "stage"):
        raise SystemExit(__doc__)
    cmd, cfg = argv[0], load(Path(argv[1]))
    if cmd == "check":
        problems = check(cfg)
        print("\n".join(problems) if problems else "config OK")
        raise SystemExit(1 if problems else 0)
    if cmd == "get":
        print(get(cfg, argv[2]))
    elif cmd == "write":
        write(cfg, Path(argv[2]))
    elif cmd == "stage":
        print(json.dumps(stage(cfg, Path(argv[2]))))


if __name__ == "__main__":
    main(sys.argv[1:])
