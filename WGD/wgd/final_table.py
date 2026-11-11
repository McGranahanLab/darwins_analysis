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
* quality: per-tumour score in [0, 1] from step 04 (wgd.quality); quality_code
  lists the factors that lowered it.
* worst_event_quality / worst_event_quality_code: score of the GD(s) on this clone, from the
  factors that belong to one event (its samples' ploidy margin, mutation support,
  ALPACA ratio, minority coverage, tie, ambiguity of this clone, reversion below it).
  The worst event if the clone carries several; null on clones without a GD.
* n_unplaced_events: events with no clone, or placed on a PGDD cluster that is
  not an ALPACA clone (they add no GD to any row).
* n_rejected_events: GDs dropped by `--reject-below` (0 without it).

The per-event table (`build_events`, `wgd_events.tsv`) has one row per GD with a
unique `wgd_event_id` (tumour:clone:k), its source, PGDD event, quality and factors.
"""

from __future__ import annotations

import json
from itertools import combinations

import polars as pl

from wgd import quality as Q
from wgd.placement import Tree

CLASS = {"no_gd": "no_wgd", "clonal_gd_only": "clonal_wgd",
         "alpaca_called": "ploidy_only", "mut_supported": "mut_supported"}
COLUMNS = ("tumour_id", "clone", "parent", "number_of_gds_at_clone", "total_gds_per_tumour",
           "quality", "status", "class", "sub_class", "quality_code", "worst_event_quality", "worst_event_quality_code",
           "n_unplaced_events", "n_rejected_events")


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


EVENT_COLUMNS = ("wgd_event_id", "tumour_id", "clone", "parent", "gd_index_at_clone", "is_clonal", "source",
                 "pgdd_event_id", "samples", "supporting_mut_count", "alpaca_score",
                 "event_quality", "event_quality_code", "event_factors", "flags", "rejected")


def _source(e: dict) -> str:
    lb = e["placed_as"] or ""
    if e["is_clonal"]:
        return "clonal"
    if lb.startswith("fill_"):
        return "fill_forced" if (e.get("flags") or "").find("minority_coverage") >= 0 else "fill"
    if lb.startswith("multicluster_"):
        return "mutation_multicluster"
    if lb.isdigit():
        return "mutation"
    return lb.split("_")[0] or "unplaced"  # notebook rules: alpaca / ploidy / false_call


MUTATION_SOURCES = ("mutation", "mutation_multicluster")


def build_events(events: pl.DataFrame, tumours: pl.DataFrame, pgdd_tumours: pl.DataFrame,
                 ratios: pl.DataFrame, reject_below: float | None = None,
                 reject_mutation_supported: bool = False) -> pl.DataFrame:
    """One row per GD (Piotr, 2026-11-11), with a unique `wgd_event_id` = tumour:clone:k
    (k-th GD on that clone; `tumour:unplaced:k` for events not on an ALPACA clone).
    PGDD ploidy-only events superseded by the fill are not GDs and are left out.
    `rejected` = event_quality below `reject_below`. Mutation-supported GDs are never rejected
    unless `reject_mutation_supported` (Piotr, 2026-11-11: already triaged by PGDD's bootstrap)."""
    num_clonal = dict(pgdd_tumours.select("tumour_id", "num_clonal_gds").iter_rows())
    ev = {t: g for (t,), g in events.group_by("tumour_id")} if events.height else {}
    rows = []
    for r in tumours.iter_rows(named=True):
        t = r["tumour_id"]
        g = ratios.filter(pl.col("tumour_id") == t)
        te = ev.get(t)
        if g.is_empty() or r["triage"] == "no_gd" or te is None:
            continue
        tree = Tree.from_ratios(g)
        placed: dict[str, list[dict]] = {}
        unplaced: list[dict] = []
        for e in te.iter_rows(named=True):
            if e["placed_as"] == "superseded_by_fill":
                continue
            c = e.get("placed_clone")
            if c in tree.parent:
                placed.setdefault(c, []).append(e)
            else:
                unplaced.append(e)
        if r["triage"] == "clonal_gd_only":  # the count comes from num_clonal_gds (notebook [40])
            n = int(num_clonal.get(t, 0))
            base = placed.get(tree.mrca, [])
            placed = {tree.mrca: (base + [base[-1]] * n)[:n]} if base else {}
            unplaced = []
        for c, es in list(placed.items()) + [(None, unplaced)]:
            for k, e in enumerate(es, 1):
                q = e.get("event_quality")
                src = _source(e)
                rejectable = c is not None and (reject_mutation_supported or src not in MUTATION_SOURCES)
                rows.append({
                    "wgd_event_id": f"{t}:{c or 'unplaced'}:{k}", "tumour_id": t, "clone": c,
                    "parent": tree.parent.get(c) if c else None, "gd_index_at_clone": k if c else None,
                    "is_clonal": e["is_clonal"], "source": src, "pgdd_event_id": e["GD_event_id"],
                    "samples": e.get("samples"), "supporting_mut_count": e.get("supporting_mut_count"),
                    "alpaca_score": e.get("alpaca_score"), "event_quality": q,
                    "event_quality_code": e.get("event_quality_code"), "event_factors": e.get("event_factors"),
                    "flags": e.get("flags"),
                    "rejected": bool(reject_below is not None and rejectable and q is not None and q < reject_below),
                })
    schema = {"wgd_event_id": pl.Utf8, "tumour_id": pl.Utf8, "clone": pl.Utf8, "parent": pl.Utf8,
              "gd_index_at_clone": pl.Int64, "is_clonal": pl.Boolean, "source": pl.Utf8, "pgdd_event_id": pl.Utf8,
              "samples": pl.Utf8, "supporting_mut_count": pl.Int64, "alpaca_score": pl.Float64,
              "event_quality": pl.Float64, "event_quality_code": pl.Utf8, "event_factors": pl.Utf8,
              "flags": pl.Utf8, "rejected": pl.Boolean}
    return pl.DataFrame(rows, schema=schema)


def build_table(events: pl.DataFrame, tumours: pl.DataFrame, pgdd_tumours: pl.DataFrame,
                ratios: pl.DataFrame, reject_below: float | None = None,
                event_table: pl.DataFrame | None = None) -> pl.DataFrame:
    """Per-clone table from the per-event table. With `reject_below`, events whose quality is
    below it are dropped: GD counts, sub_class, worst_event_quality and the tumour score
    (Q.aggregate of the remaining events' factors + tumour-level factors) are recomputed."""
    et = event_table if event_table is not None else build_events(events, tumours, pgdd_tumours, ratios, reject_below)
    by_t = {t: g for (t,), g in et.group_by("tumour_id")} if et.height else {}
    rows = []
    for r in tumours.iter_rows(named=True):
        t = r["tumour_id"]
        g = ratios.filter(pl.col("tumour_id") == t)
        if g.is_empty():
            continue  # no ALPACA tree: reported by step 04
        tree = Tree.from_ratios(g)
        n_at = {c: 0 for c in tree.parent}
        subclonal, unplaced, n_rej = [], 0, 0
        eq: dict[str, tuple[float, str | None]] = {}
        kept_factors = []
        te = by_t.get(t)
        for e in (te.iter_rows(named=True) if te is not None else []):
            if e["clone"] is None:
                unplaced += 1
                continue
            if e["rejected"]:
                n_rej += 1
                continue
            c = e["clone"]
            n_at[c] += 1
            if not e["is_clonal"]:
                subclonal.append(c)
            q = e["event_quality"]
            if q is not None and (c not in eq or q < eq[c][0]):
                eq[c] = (q, e["event_quality_code"])
            kept_factors.append(json.loads(e["event_factors"]) if e["event_factors"] else [])
        sc = (sub_class(tree, subclonal) if r["triage"] in ("alpaca_called", "mut_supported")
              and unplaced == 0 else None)
        quality, code = r.get("quality"), r.get("quality_code")
        if reject_below is not None and r.get("tumour_factors") is not None:
            quality, code = Q.aggregate(kept_factors, json.loads(r["tumour_factors"]))
        for clone, parent in tree.parent.items():
            rows.append({
                "tumour_id": t, "clone": clone, "parent": parent,
                "number_of_gds_at_clone": n_at[clone], "total_gds_per_tumour": sum(n_at.values()),
                "quality": quality, "status": r["status"], "class": CLASS[r["triage"]],
                "sub_class": sc, "quality_code": code,
                "worst_event_quality": eq.get(clone, (None, None))[0],
                "worst_event_quality_code": eq.get(clone, (None, None))[1],
                "n_unplaced_events": unplaced, "n_rejected_events": n_rej,
            })
    schema = {"tumour_id": pl.Utf8, "clone": pl.Utf8, "parent": pl.Utf8, "number_of_gds_at_clone": pl.Int64,
              "total_gds_per_tumour": pl.Int64, "quality": pl.Float64, "status": pl.Utf8, "class": pl.Utf8,
              "sub_class": pl.Utf8, "quality_code": pl.Utf8, "worst_event_quality": pl.Float64,
              "worst_event_quality_code": pl.Utf8, "n_unplaced_events": pl.Int64, "n_rejected_events": pl.Int64}
    return pl.DataFrame(rows, schema=schema).sort("tumour_id", "clone")


def rejection_report(event_table: pl.DataFrame, full: pl.DataFrame, filtered: pl.DataFrame,
                     reject_below: float) -> pl.DataFrame:
    """One row per rejected GD, with the tumour's score and GD total before and after rejection."""
    per_t = lambda tb, sfx: tb.unique("tumour_id").select(  # noqa: E731
        "tumour_id", pl.col("quality").alias(f"tumour_quality_{sfx}"),
        pl.col("total_gds_per_tumour").alias(f"total_gds_{sfx}"), pl.col("quality_code").alias(f"tumour_quality_code_{sfx}"))
    return (event_table.filter("rejected")
            .select("wgd_event_id", "tumour_id", "clone", "gd_index_at_clone", "is_clonal", "source", "pgdd_event_id",
                    "samples", "event_quality", "event_quality_code")
            .with_columns(pl.lit(reject_below).alias("reject_below"))
            .join(per_t(full, "before"), on="tumour_id", how="left")
            .join(per_t(filtered, "after"), on="tumour_id", how="left")
            .sort("tumour_id", "wgd_event_id"))
