"""Write a synthetic tx_data store for the WGS pipeline (every table it reads, including the proposed
ones tx_data doesn't have yet), for the end-to-end test. Synthetic ids only.

  SYN0001-Tumour1  2 samples, tree 1 -> {2, 3}; both samples doubled; truncal mutations at 2 copies (clonal GD)
  SYN0002-Tumour1  2 samples, tree 1 -> {2, 3}; R1 doubled, R2 not; cluster 2 (R1) at 2 copies (subclonal GD)
  SYN0003-Tumour1  1 sample, single clone (tree 1 1); doubled
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

GENOME = [("1", 0, 100_000_000), ("2", 0, 60_000_000), ("3", 0, 40_000_000)]   # 200 Mb


def build(out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    tumours = {
        "SYN0001-Tumour1": dict(tree=[(1, 2), (1, 3)], samples={"R1": 2, "R2": 3}, doubled={"R1": True, "R2": True},
                                mut_copy={1: {"R1": 2, "R2": 2}, 2: {"R1": 1, "R2": 0}, 3: {"R1": 0, "R2": 1}}),
        "SYN0002-Tumour1": dict(tree=[(1, 2), (1, 3)], samples={"R1": 2, "R2": 3}, doubled={"R1": True, "R2": False},
                                mut_copy={1: {"R1": 1, "R2": 1}, 2: {"R1": 2, "R2": 0}, 3: {"R1": 0, "R2": 1}}),
        "SYN0003-Tumour1": dict(tree=[(1, 1)], samples={"R1": 1}, doubled={"R1": True},
                                mut_copy={1: {"R1": 2}}),
    }
    seg, mut, edges, cp, alp, rat, info = [], [], [], [], [], [], []
    qc = ["patient\ttumour_id\tselected_solution\tFINALIZED"]
    for t, d in tumours.items():
        pat = t.split("-")[0]
        name = lambda s: f"{pat}_SU_T1.{s}"  # noqa: E731
        hashed = lambda s: f"{pat}_SU_T1-{s}--abc123"  # noqa: E731
        for s, dbl in d["doubled"].items():
            for chrom, start, end in GENOME:   # chr1 + chr2 doubled (80%) when doubled, else all 1|1
                a = 2.0 if dbl and chrom in ("1", "2") else 1.0
                seg.append({"tumour_id": t, "sample": name(s), "sample_name_hash": hashed(s), "chr": chrom,
                            "start": start + 1, "end": end, "cn_major": a, "cn_minor": a})
        for p, c in d["tree"]:
            edges.append({"tumour_id": t, "parent_cluster": p, "child_cluster": c, "n_trees": 1})
        for cl in sorted({x for e in d["tree"] for x in e}):
            for s in d["samples"]:
                info.append({"tumour_id": t, "cluster_id": cl, "sample": name(s), "truncal": cl == 1,
                             "tree_clust": True, "mean_ccf": 99.0 if cl == 1 else 40.0})
        qc.append(f"{pat}\t{t.split('Tumour')[1]}\tTRUE\tTRUE")
        clusters = sorted({x for e in d["tree"] for x in e})
        for cl in clusters:
            for i in range(30):
                mid = f"{pat}:1:{cl * 1_000_000 + i * 1000}:C:T"
                for s in d["samples"]:
                    mc = d["mut_copy"][cl][s]
                    a = 2.0 if d["doubled"][s] else 1.0
                    mut.append({"tumour_id": t, "sample": name(s), "chr": "1", "pos": cl * 1_000_000 + i * 1000,
                                "ref": "C", "alt": "T", "cluster": cl, "mut_copy": float(mc), "cn_a": a,
                                "cn_b": a, "ploidy": 4.0 if d["doubled"][s] else 2.0, "mutation_id": mid,
                                "tree_cluster": cl})
        for s, leaf in d["samples"].items():
            for cl in clusters:
                cp.append({"tumour_id": t, "clone": f"clone{cl}", "sample_id": hashed(s),
                           "proportion": 1.0 if cl == leaf else 0.0})
        parent = {f"clone{c}": f"clone{p}" for p, c in d["tree"] if p != c}
        parent[f"clone{clusters[0]}"] = "diploid"
        for c, p in parent.items():
            dbl_clone = (t == "SYN0001-Tumour1" and c == "clone1") or (t == "SYN0002-Tumour1" and c == "clone2") \
                or t == "SYN0003-Tumour1"
            rat.append({"tumour_id": t, "clone": c, "parent": p, "ratio_A": 2.0 if dbl_clone else 1.0,
                        "ratio_B": 2.0 if dbl_clone else 1.0, "ratio_total": 2.0 if dbl_clone else 1.0})
            for chrom, start, end in GENOME:
                inherited = t == "SYN0001-Tumour1" or c == "clone2" or t == "SYN0003-Tumour1"
                cn = 2 if inherited and chrom in ("1", "2") and (t != "SYN0002-Tumour1" or c == "clone2") else 1
                alp.append({"tumour_id": t, "segment": f"{chrom}_{start}_{end}", "clone": c, "pred_CN_A": cn,
                            "pred_CN_B": cn})
    pl.DataFrame(seg).write_parquet(out / "wgs_consensus_segments.parquet")
    pl.DataFrame(mut).write_parquet(out / "wgs_conipher_mutations.parquet")
    pl.DataFrame(edges).write_parquet(out / "wgs_conipher_tree_edges.parquet")
    pl.DataFrame(cp).write_parquet(out / "wgs_clone_proportions.parquet")
    pl.DataFrame(alp).write_parquet(out / "wgs_alpaca.parquet")
    pl.DataFrame(rat).write_parquet(out / "wgs_alpaca_wgd_ratios.parquet")
    pl.DataFrame(info).write_parquet(out / "wgs_conipher_cluster_info.parquet")
    (out / "manual_qc.tsv").write_text("\n".join(qc) + "\n")
    return out


if __name__ == "__main__":
    print(build(Path(sys.argv[1])))
