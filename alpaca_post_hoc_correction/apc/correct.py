"""Correct ALPACA clone copy number with agreed WGD calls (Piotr, 2026-11-12/13).

ALPACA knows nothing about genome doubling. It infers each clone's allele-specific integer copy
number by parsimony over the tree. A clone that is absent from every sample (an extinct
ancestor) is constrained only by its relatives. If two or more of its children doubled
*independently* (parallel GDs), the cheapest explanation for ALPACA is one shared gain on the
absent parent. The parent then looks doubled (WGD ratio ~2 vs its own parent) and the doubled
children look unchanged (~1). Once the WGD pipeline has called the GDs on the children, the
parent's copy number can be re-inferred knowing that it did *not* double.

Rule (one clone at a time, deepest first so corrected values feed the next level up):
  a clone is corrected when it is absent from every sample, carries no GD in the calls, and at
  least one child carries a GD. For every segment and allele, its new copy number x minimises
      sum over children c  |CN_c - 2^g_c * x|   +   |x - CN_parent|
  where g_c is the number of GDs on child c (each GD doubles what the child inherits) and
  CN_parent is the parent's copy number (1 for `diploid`). Ties keep ALPACA's value if it is
  among the minimisers, else the smallest. Clones with cells in a sample are observed and are
  never changed.

Ratios are then recomputed exactly as ALPACA's `analysis.calculate_wgd_ratios`
(length-weighted mean of child / parent CN per allele and in total; segments whose parent
CN is 0 are skipped; the MRCA's parent is diploid 1|1).
"""

from __future__ import annotations

import numpy as np
import polars as pl

ALLELES = ("pred_CN_A", "pred_CN_B")


# ---- tree helpers -----------------------------------------------------------------------

def depth(parent: dict[str, str], c: str) -> int:
    d = 0
    while parent.get(c, "diploid") != "diploid" and d < 10_000:
        c, d = parent[c], d + 1
    return d


def children(parent: dict[str, str], c: str) -> list[str]:
    return sorted(k for k, p in parent.items() if p == c)


# ---- the correction ---------------------------------------------------------------------

def clones_to_correct(parent: dict[str, str], gds: dict[str, int], present: set[str]) -> list[str]:
    """Absent clones without a GD that have at least one child with a GD, deepest first."""
    out = [c for c in parent if c not in present and not gds.get(c, 0)
           and any(gds.get(ch, 0) for ch in children(parent, c))]
    return sorted(out, key=lambda c: (-depth(parent, c), c))


def best_cn(child_cn: list[int], child_mult: list[int], parent_cn: int, original: int) -> int:
    """argmin_x sum_c |CN_c - m_c x| + |x - CN_parent| over x = 0..max; ties: ALPACA's value, then smallest."""
    hi = max([parent_cn, original] + list(child_cn)) + 1
    xs = np.arange(hi + 1)
    cost = np.abs(xs - parent_cn).astype(float)
    for cn, m in zip(child_cn, child_mult):
        cost += np.abs(cn - m * xs)
    best = np.flatnonzero(cost == cost.min())
    return int(original) if original in best else int(best[0])


def correct_tumour(cn: pl.DataFrame, parent: dict[str, str], gds: dict[str, int],
                   present: set[str]) -> tuple[pl.DataFrame, list[dict]]:
    """cn: one tumour's ALPACA rows (segment, clone, pred_CN_A, pred_CN_B, ...).
    Returns (corrected rows with `pred_CN_*_original` and `corrected`, one record per corrected clone)."""
    wide = {(s, c): [a, b] for s, c, a, b in cn.select("segment", "clone", *ALLELES).iter_rows()}
    original = {k: list(v) for k, v in wide.items()}
    segments = sorted({s for s, _ in wide})
    log = []
    for p in clones_to_correct(parent, gds, present):
        kids = [ch for ch in children(parent, p) if any((s, ch) in wide for s in segments)]
        if not any((s, p) in wide for s in segments) or not kids:
            continue
        changed = 0
        for s in segments:
            if (s, p) not in wide:
                continue
            for i in range(2):
                par = 1 if parent[p] == "diploid" else wide.get((s, parent[p]), [1, 1])[i]
                new = best_cn([wide[(s, ch)][i] for ch in kids if (s, ch) in wide],
                              [2 ** gds.get(ch, 0) for ch in kids if (s, ch) in wide], par, wide[(s, p)][i])
                if new != wide[(s, p)][i]:
                    wide[(s, p)][i] = new
                    changed += 1
        log.append({"clone": p, "parent": parent[p], "children": ",".join(kids),
                    "children_with_gd": ",".join(f"{ch}({gds[ch]})" for ch in kids if gds.get(ch, 0)),
                    "allele_segments_changed": changed})
    rows = cn.with_columns(
        pl.Series("pred_CN_A", [wide[(s, c)][0] for s, c in cn.select("segment", "clone").iter_rows()], dtype=pl.Int64),
        pl.Series("pred_CN_B", [wide[(s, c)][1] for s, c in cn.select("segment", "clone").iter_rows()], dtype=pl.Int64),
        pl.Series("pred_CN_A_original", [original[(s, c)][0] for s, c in cn.select("segment", "clone").iter_rows()],
                  dtype=pl.Int64),
        pl.Series("pred_CN_B_original", [original[(s, c)][1] for s, c in cn.select("segment", "clone").iter_rows()],
                  dtype=pl.Int64),
    )
    rows = rows.with_columns(((pl.col("pred_CN_A") != pl.col("pred_CN_A_original"))
                              | (pl.col("pred_CN_B") != pl.col("pred_CN_B_original"))).alias("corrected"))
    # keep tx_data's derived columns consistent
    derived = {"total_cn": pl.col("pred_CN_A") + pl.col("pred_CN_B"),
               "minor_cn": pl.min_horizontal(*ALLELES), "major_cn": pl.max_horizontal(*ALLELES),
               "is_homdel": (pl.col("pred_CN_A") + pl.col("pred_CN_B")) == 0,
               "is_loh": (pl.min_horizontal(*ALLELES) == 0) & (pl.max_horizontal(*ALLELES) > 0)}
    rows = rows.with_columns([e.alias(k) for k, e in derived.items() if k in rows.columns])
    return rows, log


# ---- WGD ratios (port of ALPACA analysis.calculate_wgd_ratios) -----------------------------

def segment_length(segment: str) -> int:
    _, start, end = segment.split("_")
    return int(end) - int(start)


def wgd_ratios(cn: pl.DataFrame, parent: dict[str, str]) -> pl.DataFrame:
    """Per clone: length-weighted mean of child/parent CN for allele A, B and total (one tumour)."""
    t = cn["tumour_id"][0]
    wide = {(s, c): (a, b) for s, c, a, b in cn.select("segment", "clone", *ALLELES).iter_rows()}
    segments = sorted({s for s, _ in wide})
    rows = []
    for c in dict.fromkeys(cn["clone"].to_list()):
        p = parent[c]
        acc = {k: [0.0, 0, 0.0] for k in ("A", "B", "total")}
        for s in segments:
            if (s, c) not in wide:
                continue
            ca, cb = wide[(s, c)]
            pa, pb = (1, 1) if p == "diploid" else wide[(s, p)]
            w = float(segment_length(s))
            for k, ch, par in (("A", ca, pa), ("B", cb, pb), ("total", ca + cb, pa + pb)):
                if par > 0:
                    acc[k][0] += ch / par * w
                    acc[k][1] += 1
                    acc[k][2] += w
        rows.append({"tumour_id": t, "clone": c, "parent": p,
                     **{f"ratio_{k}": (v[0] / v[2] if v[2] else float("nan")) for k, v in acc.items()},
                     **{f"segments_used_{k}": v[1] for k, v in acc.items()},
                     **{f"weight_sum_{k}": v[2] for k, v in acc.items()}})
    return pl.DataFrame(rows)


# ---- cohort driver ---------------------------------------------------------------------------

def run(alpaca: pl.DataFrame, clone_proportions: pl.DataFrame, wgd_calls: pl.DataFrame,
        ratios: pl.DataFrame | None = None, min_quality: float | None = None,
        threshold: float = 1.5) -> dict[str, pl.DataFrame]:
    """Correct every tumour present in both ALPACA and the WGD calls.

    `ratios` (ALPACA `wgd_ratio_scores`): only for the tree; otherwise the parents come from the
    WGD calls. `min_quality`: use only tumours whose call quality is at least this (the "agreed"
    calls); tumours below it are left as ALPACA gave them and listed as skipped.
    """
    present = {(t, c) for t, c, v in clone_proportions.select("tumour_id", "clone", "proportion").iter_rows()
               if v is not None and v > 0}
    tree_src = ratios if ratios is not None else wgd_calls
    trees = {t: dict(zip(g["clone"], g["parent"])) for (t,), g in tree_src.group_by("tumour_id")}
    calls = {t: g for (t,), g in wgd_calls.group_by("tumour_id")}
    out_cn, out_before, out_after, log, skipped = [], [], [], [], []
    for (t,), cn in alpaca.group_by("tumour_id", maintain_order=True):
        g = calls.get(t)
        parent = trees.get(t)
        reason = None
        if g is None:
            reason = "no WGD calls"
        elif parent is None:
            reason = "no tree"
        elif min_quality is not None and "quality" in g.columns and (g["quality"].max() or 0) < min_quality:
            reason = f"call quality below {min_quality:g}"
        if reason:
            skipped.append({"tumour_id": t, "reason": reason})
            out_cn.append(cn.with_columns(pl.col("pred_CN_A").alias("pred_CN_A_original"),
                                          pl.col("pred_CN_B").alias("pred_CN_B_original"),
                                          pl.lit(False).alias("corrected")))
            continue
        gds = {c: int(n) for c, n in zip(g["clone"], g["number_of_gds_at_clone"]) if n}
        pres = {c for (tt, c) in present if tt == t}
        before = wgd_ratios(cn, parent)
        fixed, entries = correct_tumour(cn, parent, gds, pres)
        after = wgd_ratios(fixed, parent)
        out_cn.append(fixed)
        out_before.append(before)
        out_after.append(after)
        rb = {r["clone"]: r for r in before.iter_rows(named=True)}
        ra = {r["clone"]: r for r in after.iter_rows(named=True)}
        for e in entries:
            seg = fixed.filter(pl.col("clone") == e["clone"])
            bp = seg["segment"].map_elements(segment_length, return_dtype=pl.Int64)
            frac = float(bp.filter(seg["corrected"]).sum() / max(bp.sum(), 1))
            p, kids = e["clone"], e["children"].split(",")
            log.append({"tumour_id": t, **e, "bp_fraction_changed": round(frac, 4),
                        "ratio_A_before": rb[p]["ratio_A"], "ratio_A_after": ra[p]["ratio_A"],
                        "ratio_B_before": rb[p]["ratio_B"], "ratio_B_after": ra[p]["ratio_B"],
                        "parent_looked_doubled": max(rb[p]["ratio_A"], rb[p]["ratio_B"]) >= threshold,
                        "children_ratio_before": ",".join(f"{k}:{max(rb[k]['ratio_A'], rb[k]['ratio_B']):.2f}" for k in kids),
                        "children_ratio_after": ",".join(f"{k}:{max(ra[k]['ratio_A'], ra[k]['ratio_B']):.2f}" for k in kids)})
    rb = pl.concat(out_before) if out_before else pl.DataFrame()
    ra = pl.concat(out_after) if out_after else pl.DataFrame()
    comp = (rb.select("tumour_id", "clone", "parent", "ratio_A", "ratio_B", "ratio_total")
              .join(ra.select("tumour_id", "clone", pl.col("ratio_A").alias("ratio_A_corrected"),
                              pl.col("ratio_B").alias("ratio_B_corrected"),
                              pl.col("ratio_total").alias("ratio_total_corrected")), on=["tumour_id", "clone"])
            if rb.height else pl.DataFrame())
    return {"alpaca_corrected": pl.concat(out_cn, how="diagonal_relaxed"),
            "wgd_ratios_corrected": ra, "ratio_comparison": comp,
            "corrections": pl.DataFrame(log, infer_schema_length=None),
            "skipped": pl.DataFrame(skipped, schema={"tumour_id": pl.Utf8, "reason": pl.Utf8})}
