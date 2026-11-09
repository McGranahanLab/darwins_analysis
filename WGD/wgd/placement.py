"""Step 04: place PGDD GD events on ALPACA clones, porting `templates/previous_analysis.ipynb` as is.

Stages follow `docs/previous_analysis_logic.md` section 3:

  A  tumour QC (manual QC table + truncal CCF), optional
  B  triage from GDs_per_tumour: no_gd, clonal_gd_only, alpaca_called (no mutation support), mut_supported
  C  alpaca_called tumours: each subclonal event -> highest-ratio candidate clone (greedy, file order)
  E  2-cluster events whose clusters are a parent-child pair -> the child
  F  clonal events -> MRCA
  G  remaining single-sample ploidy-only events -> most ancestral allowed clone
  H  QC flag: max GDs on a root-leaf path vs max sample num_gds (summary)
  I  QC flag: reversion (a descendant clone with median pred_CN_A == 1)
  K1 supporting mutation count per mutation-supported event (max num_cn2 over passing rows)

`placed_as` keeps the notebook's string labels (`alpaca_clone3`, `clonal_clone1`,
`multicluster_clone4`, `ploidy_clone2`, plain `3` for a single supporting
cluster, `false_call`, `alpaca_None`), so events can be checked against the
notebook cell by cell. `placed_clone` is the clone, or null if unplaced.

Deliberate deviations (the notebook would crash or misbehave; see the doc):
* G: when a tumour has no other events, the disallowed ancestors are empty
  (the notebook reused a stale variable from the previous tumour, section 5 item 4).
* G, kept as is: other events labelled with a prefix (`clonal_clone1`) add no
  ancestors or descendants to the disallowed set (the notebook looked up the
  label, not the clone), and other events with a plain cluster label (`4`) do
  not disallow their own clone (it is listed as `4`, not `clone4`), so a
  ploidy-only event can land on the same clone as a mutation-supported one.
* H, K1: events the notebook could not parse (`alpaca_None`, unresolved
  multi-cluster `3,6`, a cluster missing from mut_counts) are skipped or get a
  null count instead of raising.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

ALPACA_WGD_THRESHOLD = 1.5
MAJORITY = 0.5  # a clone's subtree "carries" a sample when it makes up more than this fraction of it
TRIAGE = ("no_gd", "clonal_gd_only", "alpaca_called", "mut_supported")


# ---- ids --------------------------------------------------------------------

def sample_key(sample: str) -> str:
    """Region stem used to match samples across sources: '.' -> '-', then drop the `--hash` suffix.

    The old PGDD run wrote R-mangled names (`LTX0001_SU_T1.R1..hash`); the notebook
    matched them by turning '-' into '.' in cp_table. Replace first, then split.
    """
    s = sample.strip().replace(".", "-")
    return s.split("--")[0]


def strip_clone(clusters: str | None) -> str | None:
    """Old PGDD runs wrote `clone3,clone6`; the current one writes `3,6`."""
    if clusters is None:
        return None
    return ",".join(c.strip().removeprefix("clone") for c in clusters.split(","))


def clone_of(label: str | None) -> str | None:
    """Clone a notebook label points to, or None if the event is not placed."""
    if label is None or label in ("false_call", "alpaca_None", "no_feasible_clone"):
        return None
    if "_" in label:
        return label.split("_", 1)[1]
    if re.fullmatch(r"\d+", label):
        return f"clone{label}"
    return None  # unresolved multi-cluster, e.g. "3,6"


# ---- tree -------------------------------------------------------------------

@dataclass
class Tree:
    parent: dict[str, str]  # clone -> parent ("diploid" for the MRCA); insertion order = file order

    @classmethod
    def from_ratios(cls, ratios: pl.DataFrame) -> "Tree":
        return cls(dict(ratios.select("clone", "parent").iter_rows()))

    @property
    def mrca(self) -> str | None:
        return next((c for c, p in self.parent.items() if p == "diploid"), None)

    def ancestors(self, clone: str) -> list[str]:
        """Parent chain up to and including 'diploid' (notebook get_ancestors)."""
        out = []
        while clone in self.parent:
            clone = self.parent[clone]
            out.append(clone)
        return out

    def children(self, clone: str) -> list[str]:
        return [c for c, p in self.parent.items() if p == clone]

    def descendants(self, clone: str) -> list[str]:
        out, queue = [], self.children(clone)
        while queue:
            c = queue.pop(0)
            out.append(c)
            queue.extend(self.children(c))
        return out

    def subtree(self, clone: str) -> list[str]:
        return [clone] + self.descendants(clone)

    def leaves(self) -> list[str]:
        parents = set(self.parent.values())
        return [c for c in self.parent if c not in parents]

    def lineages(self, clones: list[str]) -> list[list[str]]:
        """notebook order_clones_by_phylogeny: per terminal node of the sub-tree, walk up within `clones`."""
        sub = {c: p for c, p in self.parent.items() if c in clones}
        terminals = [c for c in sub if c not in set(sub.values())]
        out = []
        for t in terminals:
            lin = [t]
            while len(lin) < len(clones):
                p = sub.get(lin[-1])
                if p in clones and p not in lin:
                    lin.append(p)
                else:
                    break
            out.append(lin)
        return out


# ---- A: QC ------------------------------------------------------------------

def qc_pass(qc_table: Path, release_dir: Path, ccf_threshold: float = 90) -> tuple[set[str], set[str]]:
    """(tx_QC tumours, tx_QC & ccf_QC tumours). notebook get_qc_df + additional_ccf_qc.

    `meanCCF` in clusterInfo.txt is in percent (0 to slightly above 100), so the
    cut-off is 90 (Piotr, 2026-11-09; the notebook's 0.9 assumed fractions).
    Missing meanCCF occurs only on clusters off the tree; on a truncal tree
    cluster it is an error.
    """
    qc = pl.read_csv(qc_table, separator="\t", infer_schema_length=None)
    qc = qc.filter(pl.col("selected_solution").cast(pl.Utf8).str.to_lowercase() == "true").with_columns(
        (pl.col("patient") + "-Tumour" + pl.col("tumour_id").cast(pl.Float64).cast(pl.Int64).cast(pl.Utf8))
        .alias("patient_tumour"))
    finalized = (qc.group_by("patient_tumour")
                   .agg((pl.col("FINALIZED").cast(pl.Utf8).str.to_lowercase() == "true").all().alias("ok"))
                   .filter("ok")["patient_tumour"].to_list())
    passed = set()
    for pt in finalized:
        patient, tumour = pt.split("-")
        f = Path(release_dir) / patient / "mutation_trees" / f"tumour_{tumour.replace('Tumour', '')}" / "clusterInfo.txt"
        if not f.exists():
            continue  # as in the notebook: no clusterInfo -> not passed
        # pandas (notebook) reads "NA" as missing; polars needs it spelled out, or meanCCF becomes a string
        ci = (pl.read_csv(f, separator="\t", infer_schema_length=None, null_values=["NA", "NaN", ""])
                .drop_nulls(["treeClust", "truncal"])
                .with_columns(pl.col("meanCCF").cast(pl.Float64, strict=False)))
        as_bool = lambda c: pl.col(c).cast(pl.Utf8).str.to_lowercase().is_in(["true", "1"])  # noqa: E731
        trunk = ci.filter(as_bool("treeClust") & as_bool("truncal"))
        if trunk["meanCCF"].null_count():
            raise ValueError(f"{f}: missing meanCCF on a truncal tree cluster")
        if trunk.height and trunk["meanCCF"].min() > ccf_threshold:
            passed.add(pt)
    return set(finalized), passed


# ---- B: triage --------------------------------------------------------------

def triage(tumours: pl.DataFrame, events: pl.DataFrame) -> dict[str, str]:
    """tumour_id -> one of TRIAGE (notebook cells [11], [14])."""
    mut = set(events.group_by("tumour_id").agg(pl.col("is_subclonal_mutation_supported").any().alias("m"))
                    .filter("m")["tumour_id"])
    out = {}
    for r in tumours.iter_rows(named=True):
        t = r["tumour_id"]
        if r["First_GD"] == "No GD":
            out[t] = "no_gd"
        elif r["First_GD"] == "Clonal" and r["num_subclonal_gds"] == 0:
            out[t] = "clonal_gd_only"
        else:
            out[t] = "mut_supported" if t in mut else "alpaca_called"
    return out


# ---- C, feasible-set rule ------------------------------------------------------------

def feasible_clones(tree: Tree, props: dict[str, dict[str, float]], event_samples: list[str],
                    all_samples: list[str], majority: float = MAJORITY) -> list[str]:
    """Clones that can carry a subclonal event from tree + clone proportions alone.

    `props[sample][clone]` = proportion. A clone qualifies when its subtree makes up
    more than `majority` of every event sample and at most `majority` of every other
    sample of the tumour (the doubling must explain the doubled samples and must not
    reach samples without it). Returned in tree (file) order.
    """
    others = [s for s in all_samples if s not in event_samples]
    out = []
    for c in tree.parent:
        sub = set(tree.subtree(c))
        frac = {s: sum(v for k, v in props.get(s, {}).items() if k in sub) for s in all_samples}
        if all(frac.get(s, 0) > majority for s in event_samples) and all(frac.get(s, 0) <= majority for s in others):
            out.append(c)
    return out


# ---- the pipeline ---------------------------------------------------------------

@dataclass
class Placement:
    events: pl.DataFrame                    # one row per event with placed_as / placed_clone / flags
    tumours: pl.DataFrame                   # one row per tumour: triage, status, flags
    log: list[str] = field(default_factory=list)


def place_events(events: pl.DataFrame, tumours: pl.DataFrame, mut_counts: pl.DataFrame,
                 ratios: pl.DataFrame, clone_proportions: pl.DataFrame,
                 summary: pl.DataFrame | None = None, alpaca_cn: pl.DataFrame | None = None,
                 qc_tumours: set[str] | None = None, score: str = "ratio_A",
                 ploidy_rule: str = "feasible", regions: pl.DataFrame | None = None) -> Placement:
    """Run stages B–I and K1 on QC-passed tumours.

    `score`: ALPACA column used in C. `ploidy_rule`: `feasible` (default, 2026-11-09)
    or `notebook` (the original greedy rule C, for comparison). `regions`
    (pgdd_GDs_per_region) defines a tumour's samples for the feasible rule;
    without it the clone_proportions samples are used.
    """
    log: list[str] = []
    if score == "max_AB":
        ratios = ratios.with_columns(pl.max_horizontal("ratio_A", "ratio_B").alias("max_AB"))
    if qc_tumours is not None:
        tumours = tumours.filter(pl.col("tumour_id").is_in(list(qc_tumours)))
        events = events.filter(pl.col("tumour_id").is_in(list(qc_tumours)))
        mut_counts = mut_counts.filter(pl.col("tumour_id").is_in(list(qc_tumours)))

    events = events.with_columns(pl.col("clusters").map_elements(strip_clone, return_dtype=pl.Utf8))
    tri = triage(tumours, events)
    trees = {t: Tree.from_ratios(g) for (t,), g in ratios.group_by("tumour_id", maintain_order=True)}
    ratio_rows = {(r["tumour_id"], r["clone"]): r[score] for r in ratios.iter_rows(named=True)}
    cp = clone_proportions.with_columns(pl.col("sample_id").map_elements(sample_key, return_dtype=pl.Utf8).alias("key"))
    cp_by_tumour = {t: g for (t,), g in cp.group_by("tumour_id")}
    second_gd = dict(tumours.select("tumour_id", "num_second_gd").iter_rows())

    ev = [dict(r) for r in events.iter_rows(named=True)]
    for e in ev:
        e["placed_as"] = e["clusters"]
        e["alpaca_score"] = None
        e["placement_method"] = None
        e["feasible_clones"] = None
        e["flags"] = []

    def clones_in_samples(tumour: str, samples: list[str]) -> list[str]:
        g = cp_by_tumour.get(tumour)
        if g is None:
            return []
        keys = [sample_key(s) for s in samples]
        missing = sorted(set(keys) - set(g["key"]))
        if missing:
            raise ValueError(f"{tumour}: event samples not in clone_proportions: {missing}")
        tot = g.filter(pl.col("key").is_in(keys)).group_by("clone", maintain_order=True).agg(pl.col("proportion").sum())
        return tot.filter(pl.col("proportion") > 0)["clone"].to_list()

    # C: alpaca_called tumours (no mutation support), subclonal events
    tumour_samples = {}
    if regions is not None:
        for (t,), g in regions.group_by("tumour_id"):
            tumour_samples[t] = sorted({sample_key(x) for x in g["sample_id"]})
    for t in [t for t, c in tri.items() if c == "alpaca_called"]:
        tree = trees.get(t)
        if tree is None:
            log.append(f"C {t}: no ALPACA tree, skipped")
            continue
        sub_events = [e for e in ev if e["tumour_id"] == t and not e["is_clonal"]]
        if ploidy_rule == "feasible":
            g = cp_by_tumour.get(t)
            props: dict[str, dict[str, float]] = {}
            for k, c, v in (g.select("key", "clone", "proportion").iter_rows() if g is not None else []):
                props.setdefault(k, {})[c] = v
            samples_all = tumour_samples.get(t, sorted(props))
            for e in sub_events:
                ev_keys = [sample_key(x) for x in e["samples"].split(",")]
                clones_in_samples(t, e["samples"].split(","))  # fail loudly on unknown samples
                feas = feasible_clones(tree, props, ev_keys, samples_all)
                e["feasible_clones"] = ",".join(feas) or None
                if not feas:
                    e["placed_as"], e["placement_method"] = "no_feasible_clone", "none"
                    e["flags"].append("no_feasible_clone")
                    continue
                if len(feas) == 1:
                    pick, method = feas[0], "unique"
                else:  # a chain of clones fits: the ratio breaks the tie (reuse allowed)
                    pick = max(feas, key=lambda c: -np.inf if ratio_rows.get((t, c)) is None
                               or np.isnan(ratio_rows[(t, c)]) else ratio_rows[(t, c)])
                    method = "ratio"
                e["placed_as"], e["placement_method"] = f"feasible_{pick}", method
                e["alpaca_score"] = ratio_rows.get((t, pick))
            continue
        # notebook rule C: greedy, file order, a clone is used at most once, MRCA excluded
        cand_all = sorted(((c, ratio_rows[(t, c)]) for c, p in tree.parent.items() if p != "diploid"),
                          key=lambda x: -x[1] if x[1] is not None and not np.isnan(x[1]) else np.inf)
        taken: list[str | None] = []
        for e in sub_events:
            present = clones_in_samples(t, e["samples"].split(","))
            present = set(present + sum((tree.ancestors(c) for c in present), []))
            cands = [(c, s) for c, s in cand_all if c in present and c not in taken]
            pick = cands[0] if cands else (None, None)
            taken.append(pick[0])
            e["placed_as"] = f"alpaca_{pick[0]}"
            e["alpaca_score"] = pick[1]
            e["placement_method"] = "notebook_greedy"

    # E: two-cluster parent-child events -> child
    for e in ev:
        if e["placed_as"] is None or "," not in str(e["placed_as"]):
            continue
        clusters = [f"clone{c}" for c in e["placed_as"].split(",")]
        if len(clusters) != 2:
            log.append(f"E {e['GD_event_id']}: {len(clusters)} clusters, skipped")
            continue
        tree = trees.get(e["tumour_id"])
        if tree is None:
            continue
        for child in clusters:
            if tree.parent.get(child) in clusters:
                e["placed_as"] = f"multicluster_{child}"
                e["flags"].append("multicluster_parent_child")
                break

    # F: clonal events -> MRCA
    for e in ev:
        if e["is_clonal"]:
            tree = trees.get(e["tumour_id"])
            if tree is None or tree.mrca is None:
                log.append(f"F {e['GD_event_id']}: no ALPACA tree, clonal event not placed")
                continue
            e["placed_as"] = f"clonal_{tree.mrca}"

    # K1 (cell [46], computed before G as in the notebook)
    mc = mut_counts.with_columns(pl.col("cluster_id").cast(pl.Utf8).map_elements(strip_clone, return_dtype=pl.Utf8))
    has_pass = {"discover_pass_all", "check_pass_all", "num_cn2"} <= set(mc.columns)
    if not has_pass:
        log.append("K1: mut_counts lacks discover_pass_all/check_pass_all/num_cn2; supporting counts not computed")
    for e in ev:
        e["supporting_mut_count"] = None
        if not has_pass or not e["is_subclonal_mutation_supported"]:
            continue
        c = str(e["placed_as"])
        c = c.split("_clone")[1] if "_" in c else c
        if not re.fullmatch(r"\d+", c):
            continue
        rows = mc.filter((pl.col("tumour_id") == e["tumour_id"]) & (pl.col("cluster_id") == c)
                         & (pl.col("discover_pass_all") | pl.col("check_pass_all")))
        if rows.height:
            e["supporting_mut_count"] = int(rows["num_cn2"].max())

    # G: subclonal events still without a label (snapshot of labels after C/E/F, as in the notebook)
    snapshot = {e["GD_event_id"]: e["placed_as"] for e in ev}
    for e in ev:
        if e["is_clonal"] or snapshot[e["GD_event_id"]] is not None:
            continue
        t, eid = e["tumour_id"], e["GD_event_id"]
        if second_gd.get(t, 0) > 0:
            e["flags"].append("G_skipped_second_gd")
            continue
        tree = trees.get(t)
        if tree is None:
            continue
        others = [snapshot[o["GD_event_id"]] for o in ev if o["tumour_id"] == t and o["GD_event_id"] != eid]
        f_desc, f_anc = [], []
        if others:
            if any(o is None for o in others):
                e["flags"].append("G_skipped_other_ploidy_event")
                continue
            if any("," in o for o in others):
                e["flags"].append("G_skipped_other_multicluster")
                continue
            looked_up = [f"clone{o}" if "clone" not in o else o for o in others]
            f_desc = sum((tree.descendants(c) for c in looked_up), [])
            f_anc = sum((tree.ancestors(c) for c in looked_up), [])
        forbidden = set(f_desc + [o.split("_")[-1] if "_" in o else o for o in others] + f_anc)
        if "," in e["samples"]:
            e["flags"].append("G_skipped_multi_sample")
            continue
        present = clones_in_samples(t, [e["samples"]])
        with_anc = [c for c in set(present + sum((tree.ancestors(c) for c in present), [])) if c != "diploid"]
        lineages = tree.lineages(with_anc)
        if not lineages:
            continue
        cands = []
        for lin in lineages:
            cs = []
            for c in lin:
                if c == lin[0] and c not in forbidden:
                    cs.append(c)
                    continue
                if any(d in forbidden for d in tree.descendants(c)) or c in forbidden:
                    continue
                cs.append(c)
            if cs:
                cands.append(cs)
        if not cands:
            e["placed_as"] = "false_call"
        elif len(cands) == 1:
            e["placed_as"] = f"ploidy_{cands[0][-1]}"
        elif all(c[-1] == cands[0][-1] for c in cands):
            e["placed_as"] = f"ploidy_{cands[0][-1]}"
        else:
            e["flags"].append("G_multiple_lineage_candidates")

    for e in ev:
        e["placed_clone"] = clone_of(e["placed_as"])

    # H: path GD count vs max sample num_gds
    h_mismatch = {}
    if summary is not None:
        max_gds = dict(summary.group_by("tumour_id").agg(pl.col("num_gds").max()).iter_rows())
        for t in {e["tumour_id"] for e in ev}:
            te = [e for e in ev if e["tumour_id"] == t]
            if all(e["is_clonal"] for e in te) or t not in trees:
                continue
            labels = {e["placed_as"] for e in te}
            if None in labels or "false_call" in labels:
                continue
            clones = {clone_of(lb) for lb in labels}
            if None in clones:
                continue  # notebook would raise here (alpaca_None, unresolved multi-cluster)
            tree = trees[t]
            # notebook: distinct GD clones on the path; feasible rule: events (a clone can carry several)
            per_clone = ({c: sum(clone_of(e["placed_as"]) == c for e in te) for c in clones}
                         if ploidy_rule == "feasible" else {c: 1 for c in clones})
            path_max = max(sum(per_clone.get(c, 0) for c in [leaf] + tree.ancestors(leaf)) for leaf in tree.leaves())
            if max_gds.get(t) != path_max:
                h_mismatch[t] = (max_gds.get(t), path_max)

    # I: reversion, mutation-supported tumours only
    if alpaca_cn is not None:
        med = (alpaca_cn.group_by("tumour_id", "clone").agg(pl.col("pred_CN_A").median())
                        .filter(pl.col("pred_CN_A") == 1))
        reverted = {(t, c) for t, c, _ in med.iter_rows()}
        for e in ev:
            t = e["tumour_id"]
            if tri.get(t) != "mut_supported" or e["is_clonal"] or e["placed_as"] is None or t not in trees:
                continue
            lb = e["placed_as"]
            clone = lb.split("_")[1] if "_" in lb else f"clone{lb}"
            if any((t, d) in reverted for d in trees[t].descendants(clone)):
                e["flags"].append("reversion")

    # per-tumour status (notebook queue files): no_gd / clonal_gd_only resolved,
    # alpaca_called resolved iff every event's ratio > 1.5, mut_supported unresolved
    rows = []
    for t, cls in tri.items():
        te = [e for e in ev if e["tumour_id"] == t]
        if cls in ("no_gd", "clonal_gd_only"):
            status = "resolved"
        elif cls == "alpaca_called" and ploidy_rule == "feasible":
            # resolved when every event fits exactly one clone by tree + ploidy, or the
            # ratio that picked among several fitting clones is above the threshold
            sub = [e for e in te if not e["is_clonal"]]
            ok = bool(sub) and all(e["placement_method"] == "unique" or (
                e["placement_method"] == "ratio" and e["alpaca_score"] is not None
                and e["alpaca_score"] > ALPACA_WGD_THRESHOLD) for e in sub)
            status = "resolved" if ok else "unresolved"
        elif cls == "alpaca_called":
            scores = [e["alpaca_score"] for e in te if not e["is_clonal"]]
            ok = bool(scores) and all(s is not None and s > ALPACA_WGD_THRESHOLD for s in scores)
            status = "resolved" if ok else "unresolved"
        else:
            status = "unresolved"
        rows.append({"tumour_id": t, "triage": cls, "status": status, "has_alpaca_tree": t in trees,
                     "n_events": len(te), "n_unplaced": sum(e["placed_clone"] is None for e in te),
                     "H_mismatch": t in h_mismatch,
                     "H_summary_max_gds": h_mismatch.get(t, (None, None))[0],
                     "H_tree_path_max_gds": h_mismatch.get(t, (None, None))[1]})

    ev_df = pl.DataFrame([{**{k: v for k, v in e.items() if k != "flags"}, "flags": ";".join(e["flags"]) or None}
                          for e in ev], infer_schema_length=None)
    tumour_schema = {"tumour_id": pl.Utf8, "triage": pl.Utf8, "status": pl.Utf8, "has_alpaca_tree": pl.Boolean,
                     "n_events": pl.Int64, "n_unplaced": pl.Int64, "H_mismatch": pl.Boolean,
                     "H_summary_max_gds": pl.Int64, "H_tree_path_max_gds": pl.Int64}
    return Placement(ev_df, pl.DataFrame(rows, schema=tumour_schema).sort("tumour_id"), log)
