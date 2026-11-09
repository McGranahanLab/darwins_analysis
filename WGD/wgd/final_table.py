"""Step 05: per-clone WGD table in the format of the previous release table.

One row per ALPACA clone of every QC-passed PGDD tumour:

  tumour_id, clone, parent, number_of_gds_at_clone, total_gds_per_tumour,
  quality, status, class, sub_class   (+ quality_code, n_unplaced_events)

* number_of_gds_at_clone: GDs placed on this clone (not cumulative). For
  clonal_gd_only tumours it is `num_clonal_gds` on the MRCA (notebook [40]);
  otherwise the number of events placed on the clone by step 04.
* status: from step 04 (the notebook's automatic queue status); `unresolved`
  wherever the notebook left the call to manual review.
* class: no_gd -> no_wgd, clonal_gd_only -> clonal_wgd, alpaca_called -> ploidy_only,
  mut_supported -> mut_supported.
* sub_class (not in the notebook; rule from the tx_data wgd_calls definitions):
  for ploidy_only / mut_supported tumours with every event placed:
  `single subclonal GD`, `parallel only`, `sequential only` or `parallel + sequential`
  over the subclonal GD clones. Null otherwise.
* quality (notebook [53]-[57]): start at 1 (the viewer's manual quality is not
  available); a tumour with any clone whose supporting mutation count is
  <= the cohort 5th percentile gets 0.5; tumour minimum; < 1 -> Low.
  quality_code says why (`low_supporting_mut_count`).
* n_unplaced_events: events with no clone, or placed on a PGDD cluster that is
  not an ALPACA clone (they add no GD to any row).
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import polars as pl

from wgd.placement import Tree

CLASS = {"no_gd": "no_wgd", "clonal_gd_only": "clonal_wgd",
         "alpaca_called": "ploidy_only", "mut_supported": "mut_supported"}
COLUMNS = ("tumour_id", "clone", "parent", "number_of_gds_at_clone", "total_gds_per_tumour",
           "quality", "status", "class", "sub_class", "quality_code", "n_unplaced_events")


def sub_class(tree: Tree, subclonal_clones: list[str]) -> str | None:
    if not subclonal_clones:
        return None
    if len(subclonal_clones) == 1:
        return "single subclonal GD"
    related = [a == b or a in tree.ancestors(b) or b in tree.ancestors(a)
               for a, b in combinations(subclonal_clones, 2)]
    if all(related):
        return "sequential only"
    if not any(related):
        return "parallel only"
    return "parallel + sequential"


def build_table(events: pl.DataFrame, tumours: pl.DataFrame, pgdd_tumours: pl.DataFrame,
                ratios: pl.DataFrame) -> pl.DataFrame:
    """`events`/`tumours` from step 04; `pgdd_tumours` = pgdd_GDs_per_tumour; `ratios` = alpaca_wgd_ratios."""
    num_clonal = dict(pgdd_tumours.select("tumour_id", "num_clonal_gds").iter_rows())
    ev = {t: g for (t,), g in events.group_by("tumour_id")} if events.height else {}

    # quality, notebook [54]-[56]: 5th percentile of per-clone supporting counts across the cohort
    counts = (events.filter(pl.col("supporting_mut_count").is_not_null() & pl.col("placed_clone").is_not_null())
                    .select("tumour_id", pl.col("placed_clone").alias("clone"), "supporting_mut_count")
              if events.height and "supporting_mut_count" in events.columns else pl.DataFrame())
    low_tumours = set()
    if counts.height:
        q05 = float(np.quantile(counts["supporting_mut_count"].to_numpy(), 0.05))
        low_tumours = set(counts.filter(pl.col("supporting_mut_count") <= q05)["tumour_id"])

    rows = []
    for r in tumours.iter_rows(named=True):
        t = r["tumour_id"]
        g = ratios.filter(pl.col("tumour_id") == t)
        if g.is_empty():
            continue  # no ALPACA tree: reported by step 04
        tree = Tree.from_ratios(g)
        n_at = {c: 0 for c in tree.parent}
        te = ev.get(t)
        subclonal, unplaced = [], 0
        if r["triage"] == "clonal_gd_only":
            n_at[tree.mrca] = int(num_clonal.get(t, 0))
        elif r["triage"] != "no_gd" and te is not None:
            for e in te.iter_rows(named=True):
                c = e["placed_clone"]
                if c is None or c not in n_at:  # unplaced, or a PGDD cluster that is not an ALPACA clone
                    unplaced += 1
                    continue
                n_at[c] += 1
                if not e["is_clonal"]:
                    subclonal.append(c)
        sc = (sub_class(tree, subclonal) if r["triage"] in ("alpaca_called", "mut_supported")
              and unplaced == 0 else None)
        low = t in low_tumours
        for clone, parent in tree.parent.items():
            rows.append({
                "tumour_id": t, "clone": clone, "parent": parent,
                "number_of_gds_at_clone": n_at[clone], "total_gds_per_tumour": sum(n_at.values()),
                "quality": "Low" if low else "High", "status": r["status"], "class": CLASS[r["triage"]],
                "sub_class": sc, "quality_code": "low_supporting_mut_count" if low else None,
                "n_unplaced_events": unplaced,
            })
    schema = {"tumour_id": pl.Utf8, "clone": pl.Utf8, "parent": pl.Utf8, "number_of_gds_at_clone": pl.Int64,
              "total_gds_per_tumour": pl.Int64, "quality": pl.Utf8, "status": pl.Utf8, "class": pl.Utf8,
              "sub_class": pl.Utf8, "quality_code": pl.Utf8, "n_unplaced_events": pl.Int64}
    return pl.DataFrame(rows, schema=schema).sort("tumour_id", "clone")
