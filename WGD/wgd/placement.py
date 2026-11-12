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

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from wgd import params as P
from wgd import quality as Q

ALPACA_WGD_THRESHOLD = 1.5
MAJORITY = 0.5  # a clone's subtree "carries" a sample when it makes up more than this fraction of it
MIN_COVERAGE = 0.1  # relaxed fill pass: smallest share of a sample a clone must cover
TRIAGE = ("no_gd", "clonal_gd_only", "alpaca_called", "mut_supported")


# ---- ids --------------------------------------------------------------------

def sample_key(sample: str) -> str:
    """Region stem used to match samples across sources: '.' -> '-', then drop the `--hash` suffix.

    The old PGDD run wrote R-mangled names (`LTX0001_SU_T1.R1..hash`); the notebook
    matched them by turning '-' into '.' in cp_table. Replace first, then split.
    """
    s = sample.strip().replace(".", "-")
    return s.split("--")[0]


def region_label(key: str) -> str:
    """Short region name for codes: LTX1018_SU_T1-R2 -> R2, LTX1018_SU_FLN01 -> FLN01."""
    return key.split("_")[-1].split("-")[-1]


def strip_clone(clusters: str | None) -> str | None:
    """Old PGDD runs wrote `clone3,clone6`; the current one writes `3,6`."""
    if clusters is None:
        return None
    return ",".join(c.strip().removeprefix("clone") for c in clusters.split(","))


def clone_of(label: str | None) -> str | None:
    """Clone a notebook label points to, or None if the event is not placed."""
    if label is None or label in ("false_call", "alpaca_None", "superseded_by_fill"):
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


# ---- ploidy fill (tumour-level rule replacing C and G) ---------------------------

def sample_levels(tree: Tree, props: dict[str, dict[str, float]], gds: dict[str, int],
                  majority: float | None = None) -> dict[str, int]:
    """GD level of each sample under a placement `gds` (clone -> number of GDs on it).

    A clone's count is the GDs on its root path. A sample's level is the largest k
    such that clones with count >= k make up more than `majority` of the sample
    (the sample-level analogue of PGDD's >= 50%-of-genome rule).
    """
    majority = MAJORITY if majority is None else majority
    count = {c: sum(gds.get(a, 0) for a in [c] + tree.ancestors(c)) for c in tree.parent}
    out = {}
    for s, p in props.items():
        k = 0
        while sum(v for c, v in p.items() if count.get(c, 0) >= k + 1) > majority:
            k += 1
        out[s] = k
    return out


def _cell_fit(tree: Tree, props: dict[str, dict[str, float]], targets: dict[str, int],
              gds: dict[str, int]) -> tuple[float, float]:
    """(cells with more GDs than their sample's target, cells exactly at it), summed over samples."""
    count = {c: sum(gds.get(a, 0) for a in [c] + tree.ancestors(c)) for c in tree.parent}
    over = sum(v for s, p in props.items() for c, v in p.items() if count.get(c, 0) > targets[s])
    at = sum(v for s, p in props.items() for c, v in p.items() if count.get(c, 0) == targets[s])
    return over, at


@dataclass
class Fill:
    added: list[dict]            # one per added GD: clone, samples raised, spill, tie info
    remaining: dict[str, int]    # sample -> GDs still missing
    over: dict[str, int]         # sample -> GDs too many from the fixed (clonal + mutation) placement


def fill_ploidy(tree: Tree, props: dict[str, dict[str, float]], targets: dict[str, int],
                fixed: dict[str, int], score: dict[str, float], majority: float | None = None) -> Fill:
    """Add GDs until every sample's level reaches its PGDD `num_gds`, never exceeding it.

    Each step adds one GD on the clone that
      (1) raises the most samples that still need a GD, without raising any sample
          past its target;
      (2) leaves the fewest cells *above* their sample's target (`over_cells`): a GD
          must not over-double clones in samples that already have enough
          (LTX0073: clone14, not clone5, which would put clone12 in R4 at 3 GDs vs 2);
      (3) leaves the most cells exactly *at* their sample's target (`at_target`): with
          (2) equal, prefer the more ancestral clone that also doubles minor clones of
          already-doubled samples (LTX0501: clone6, not clone9, as clone6 is 8% of R1,
          which is at 1 GD through clone4);
      (4) is the most ancestral (fewest ancestors): when the data can't tell clones apart,
          e.g. a chain through an extinct clone, the GD goes as high as the data allow
          (LTX1338: clone9, not its only child clone15);
      (5) has the highest ALPACA score; (6) comes first in the tree.
    Cell fractions are summed over samples. A clone can take several GDs.
    """
    majority = MAJORITY if majority is None else majority
    gds = dict(fixed)
    lv = sample_levels(tree, props, gds, majority)
    over = {s: lv[s] - targets[s] for s in props if lv[s] > targets[s]}
    added = []
    for _ in range(sum(max(0, targets[s] - lv[s]) for s in props) + 1):
        need = {s for s in props if lv[s] < targets[s]}
        if not need:
            break
        cands = []
        for c in tree.parent:
            trial = sample_levels(tree, props, {**gds, c: gds.get(c, 0) + 1}, majority)
            raised = {s for s in props if trial[s] > lv[s]}
            if not raised or any(trial[s] > max(targets[s], lv[s]) for s in raised):
                continue
            over_cells, at_target = _cell_fit(tree, props, targets, {**gds, c: gds.get(c, 0) + 1})
            sc = score.get(c)
            sc = -np.inf if sc is None or np.isnan(sc) else sc
            depth = len(tree.ancestors(c))
            cands.append((-len(raised), round(over_cells, 9), -round(at_target, 9), depth, -sc, c, raised))
        if not cands:
            break
        cands.sort(key=lambda x: x[:5])  # stable: ties keep tree order
        best = cands[0]
        tied = [x[5] for x in cands if x[:4] == best[:4]]  # only the ALPACA score (or tree order) separated these
        gds[best[5]] = gds.get(best[5], 0) + 1
        added.append({"clone": best[5], "samples": sorted(best[6]), "over_cells": best[1], "at_target": -best[2],
                      "alternatives": [c for c in tied if c != best[5]]})
        lv = sample_levels(tree, props, gds, majority)
    need = {s: targets[s] - lv[s] for s in props if lv[s] < targets[s]}

    # Relaxed pass (2026-11-10, LTX1055): a sample is doubled but no clone's subtree is the
    # majority of it without over-doubling other samples. Clone proportions are then likely
    # off, so take the clone covering most of the still-undoubled samples' cells, provided no
    # cell anywhere ends above its sample's target. Each such GD is marked `relaxed` with its
    # coverage, which lowers quality (`minority_coverage`).
    for _ in range(sum(need.values())):
        if not need:
            break
        over_now, _ = _cell_fit(tree, props, targets, gds)
        cands = []
        for c in tree.parent:
            sub = set(tree.subtree(c))
            cov = {s: round(sum(v for k, v in props[s].items() if k in sub), 9) for s in need}
            cov = {s: v for s, v in cov.items() if v >= MIN_COVERAGE}
            if not cov:
                continue
            trial = {**gds, c: gds.get(c, 0) + 1}
            if _cell_fit(tree, props, targets, trial)[0] > over_now + 1e-9:
                continue
            if any(lv_t > targets[s] for s, lv_t in sample_levels(tree, props, trial, majority).items()):
                continue
            sc = score.get(c)
            sc = -np.inf if sc is None or np.isnan(sc) else sc
            cands.append((-sum(cov.values()), len(tree.ancestors(c)), -sc, c, cov))
        if not cands:
            break
        cands.sort(key=lambda x: x[:3])
        best = cands[0]
        gds[best[3]] = gds.get(best[3], 0) + 1
        added.append({"clone": best[3], "samples": sorted(best[4]), "over_cells": over_now, "at_target": None,
                      "alternatives": [x[3] for x in cands if x[:2] == best[:2] and x[3] != best[3]],
                      "relaxed": True, "coverage": best[4]})
        for s_ in best[4]:
            need[s_] -= 1
        need = {s_: n for s_, n in need.items() if n > 0}
    return Fill(added, need, over)


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
                 ploidy_rule: str = "fill", regions: pl.DataFrame | None = None,
                 params: P.Params | None = None) -> Placement:
    """Run stages B–I and K1 on QC-passed tumours.

    `score`: ALPACA column used as a tie-breaker. `ploidy_rule`: `fill` (default,
    2026-11-09) places ploidy-only GDs per tumour with `fill_ploidy`, replacing
    rules C and G; `notebook` runs the original rules C and G. `regions`
    (pgdd_GDs_per_region) gives each sample's target `num_gds`; required for `fill`.
    `params`: tunable values (`wgd.params`); defaults when None. With `fill` and a `summary`,
    the majority ploidy prior (`Q.majority_prior`) may move borderline samples' targets first.
    """
    params = params or P.Params()
    P.apply(params)
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

    # C (notebook rule only; `fill` places these events after F, together with G's)
    for t in ([t for t, c in tri.items() if c == "alpaca_called"] if ploidy_rule == "notebook" else []):
        tree = trees.get(t)
        if tree is None:
            log.append(f"C {t}: no ALPACA tree, skipped")
            continue
        sub_events = [e for e in ev if e["tumour_id"] == t and not e["is_clonal"]]
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
    for e in (ev if ploidy_rule == "notebook" else []):
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

    # fill: per tumour, fix clonal + mutation-supported GDs, then add the fewest GDs that
    # bring every sample to its PGDD num_gds. PGDD's ploidy-only events are superseded.
    fill_info: dict[str, Fill | str] = {}
    prior_info: dict[str, dict] = {}   # tumour -> {"applied": bool, "moves": {...}, "targets": {...}}
    fracs_all = _fracs_by_tumour(summary)
    if ploidy_rule == "fill":
        if regions is None:
            raise ValueError("ploidy_rule='fill' needs `regions` (pgdd_GDs_per_region)")
        targets_all = {}
        for t, sid, n in regions.select("tumour_id", "sample_id", "num_gds").iter_rows():
            targets_all.setdefault(t, {})[sample_key(sid)] = int(n)
        template = {k: None for k in ev[0]} if ev else {}
        for t in [t for t, c in tri.items() if c in ("alpaca_called", "mut_supported")]:
            tree = trees.get(t)
            te = [e for e in ev if e["tumour_id"] == t]
            ploidy_ev = [e for e in te if not e["is_clonal"] and e["placed_as"] is None]
            fixed_ev = [e for e in te if e not in ploidy_ev]
            if tree is None:
                fill_info[t] = "no_alpaca_tree"
                continue
            if any(clone_of(e["placed_as"]) not in tree.parent for e in fixed_ev):
                fill_info[t] = "fixed_event_not_on_tree"  # unresolved multi-cluster or off-tree cluster
                for e in fixed_ev:
                    if clone_of(e["placed_as"]) not in tree.parent:
                        e["flags"].append("fill_blocked")
                continue
            targets = targets_all.get(t, {})
            g = cp_by_tumour.get(t)
            props: dict[str, dict[str, float]] = {k: {} for k in targets}
            for k, c, v in (g.select("key", "clone", "proportion").iter_rows() if g is not None else []):
                if k in props and v > 0:
                    props[k][c] = v
            missing = sorted(k for k in targets if not props[k])
            if missing:
                raise ValueError(f"{t}: samples without clone proportions: {missing}")
            fixed: dict[str, int] = {}
            for e in fixed_ev:
                c = clone_of(e["placed_as"])
                fixed[c] = fixed.get(c, 0) + 1
            sc = {c: ratio_rows.get((t, c)) for c in tree.parent}
            res = fill_ploidy(tree, props, targets, fixed, sc)
            if params.prior_enabled and fracs_all.get(t):
                new_targets, moves = Q.majority_prior(targets, fracs_all[t], params.band,
                                                      params.prior_min_clear_fraction,
                                                      params.prior_borderline_majority)
                if moves:
                    res_prior = fill_ploidy(tree, props, new_targets, fixed, sc)
                    applied = not params.prior_require_simpler or _complexity(res_prior) < _complexity(res)
                    prior_info[t] = {"applied": applied, "moves": moves, "targets": new_targets,
                                     "strict": _solution(fixed, res), "prior": _solution(fixed, res_prior)}
                    if applied:
                        res = res_prior
            fill_info[t] = res
            for e in ploidy_ev:
                e["placed_as"], e["placement_method"] = "superseded_by_fill", "fill"
            raised = ({s_: m for s_, m in prior_info[t]["moves"].items() if m[1] > m[0]}
                      if prior_info.get(t, {}).get("applied") else {})
            for i, a in enumerate(res.added, 1):
                mine = {s_: m for s_, m in raised.items() if s_ in a["samples"]}
                ev.append({**template, "GD_event": f"fill{i}", "GD_event_id": f"{t}_fill{i}", "tumour_id": t,
                           "is_clonal": False, "is_subclonal_mutation_supported": False,
                           "samples": ",".join(a["samples"]), "placed_as": f"fill_{a['clone']}",
                           "placement_method": "fill", "alpaca_score": ratio_rows.get((t, a["clone"])),
                           "fill_alternatives": ",".join(a["alternatives"]) or None,
                           "flags": (["fill_tie"] if a["alternatives"] else [])
                                    + ([f"minority_coverage({','.join(f'{region_label(k)}:{v:.2f}' for k, v in a['coverage'].items())})"]
                                       if a.get("relaxed") else [])
                                    + ([f"prior_raised({','.join(f'{region_label(k)}:{m[2]:.3f}' for k, m in mine.items())})"]
                                       if mine else [])})
            if len(res.added) != len(ploidy_ev):
                log.append(f"fill {t}: {len(res.added)} GDs added for {len(ploidy_ev)} PGDD ploidy-only events")

    for e in ev:
        e["placed_clone"] = clone_of(e["placed_as"])
        e.setdefault("fill_alternatives", None)

    # H: path GD count vs max sample num_gds
    h_mismatch = {}
    if summary is not None:
        max_gds = dict(summary.group_by("tumour_id").agg(pl.col("num_gds").max()).iter_rows())
        max_gds.update({t: max(p["targets"].values()) for t, p in prior_info.items() if p["applied"]})
        for t in {e["tumour_id"] for e in ev}:
            te = [e for e in ev if e["tumour_id"] == t]
            if all(e["is_clonal"] for e in te) or t not in trees:
                continue
            labels = {e["placed_as"] for e in te if e["placed_as"] != "superseded_by_fill"}
            if None in labels or "false_call" in labels:
                continue
            clones = {clone_of(lb) for lb in labels}
            if None in clones:
                continue  # notebook would raise here (alpaca_None, unresolved multi-cluster)
            tree = trees[t]
            # notebook: distinct GD clones on the path; fill: events (a clone can carry several)
            per_clone = ({c: sum(clone_of(e["placed_as"]) == c for e in te) for c in clones}
                         if ploidy_rule == "fill" else {c: 1 for c in clones})
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

    # quality (wgd.quality): per tumour, product of evidence factors
    qual = _tumour_quality(ev, tri, trees, cp_by_tumour, regions, summary, ratios, alpaca_cn,
                           tumours, h_mismatch, fill_info, ploidy_rule, prior_info, params)

    # per-tumour status (notebook queue files): no_gd / clonal_gd_only resolved,
    # alpaca_called resolved iff every event's ratio > 1.5, mut_supported unresolved
    rows = []
    for t, cls in tri.items():
        te = [e for e in ev if e["tumour_id"] == t]
        if cls in ("no_gd", "clonal_gd_only"):
            status = "resolved"
        elif ploidy_rule == "fill":
            # resolved when the fixed + added GDs bring every sample exactly to its num_gds
            f = fill_info.get(t)
            status = "resolved" if isinstance(f, Fill) and not f.remaining and not f.over else "unresolved"
        elif cls == "alpaca_called":
            scores = [e["alpaca_score"] for e in te if not e["is_clonal"]]
            ok = bool(scores) and all(s is not None and s > ALPACA_WGD_THRESHOLD for s in scores)
            status = "resolved" if ok else "unresolved"
        else:
            status = "unresolved"
        f = fill_info.get(t)
        rows.append({"tumour_id": t, "triage": cls, "status": status, "has_alpaca_tree": t in trees,
                     "n_events": len(te),
                     "n_unplaced": sum(e["placed_clone"] is None and e["placed_as"] != "superseded_by_fill" for e in te),
                     "fill": (None if f is None else f if isinstance(f, str) else
                              "complete" if not f.remaining and not f.over else
                              ";".join([f"missing:{s}={n}" for s, n in f.remaining.items()]
                                       + [f"over:{s}={n}" for s, n in f.over.items()])),
                     "H_mismatch": t in h_mismatch,
                     "H_summary_max_gds": h_mismatch.get(t, (None, None))[0],
                     "H_tree_path_max_gds": h_mismatch.get(t, (None, None))[1],
                     "quality": qual.get(t, (None, None, None))[0], "quality_code": qual.get(t, (None, None, None))[1],
                     "tumour_factors": qual.get(t, (None, None, None))[2],
                     "ploidy_prior": _prior_text(prior_info.get(t)),
                     "solution_strict": (prior_info.get(t) or {}).get("strict"),
                     "solution_prior": (prior_info.get(t) or {}).get("prior")})

    ev_df = pl.DataFrame([{**{k: v for k, v in e.items() if k != "flags"}, "flags": ";".join(e["flags"]) or None}
                          for e in ev], infer_schema_length=None)
    tumour_schema = {"tumour_id": pl.Utf8, "triage": pl.Utf8, "status": pl.Utf8, "has_alpaca_tree": pl.Boolean,
                     "n_events": pl.Int64, "n_unplaced": pl.Int64, "fill": pl.Utf8, "H_mismatch": pl.Boolean,
                     "H_summary_max_gds": pl.Int64, "H_tree_path_max_gds": pl.Int64,
                     "quality": pl.Float64, "quality_code": pl.Utf8, "tumour_factors": pl.Utf8,
                     "ploidy_prior": pl.Utf8, "solution_strict": pl.Utf8, "solution_prior": pl.Utf8}
    return Placement(ev_df, pl.DataFrame(rows, schema=tumour_schema).sort("tumour_id"), log)


def _fracs_by_tumour(summary: pl.DataFrame | None) -> dict[str, dict[str, tuple]]:
    """tumour -> sample key -> (genome_frac_nMaj2, nMaj3, nMaj4), from the step 03 summary."""
    out: dict[str, dict[str, tuple]] = {}
    if summary is not None and {"genome_frac_nMaj2", "genome_frac_nMaj3", "genome_frac_nMaj4"} <= set(summary.columns):
        for r in summary.iter_rows(named=True):
            out.setdefault(r["tumour_id"], {})[sample_key(r["sample_name_hash"])] = (
                r["genome_frac_nMaj2"], r["genome_frac_nMaj3"], r["genome_frac_nMaj4"])
    return out


def _complexity(f: Fill) -> tuple[int, int, int]:
    """Smaller is simpler: (GDs still missing or in excess, forced relaxed-pass GDs, GDs added)."""
    return (sum(f.remaining.values()) + sum(f.over.values()),
            sum(bool(a.get("relaxed")) for a in f.added), len(f.added))


def _solution(fixed: dict[str, int], f: Fill) -> str:
    """GDs per clone after a fill, e.g. `clone3=1,clone4=1,clone5=1*` (* = forced), `none`."""
    gds = dict(fixed)
    forced = set()
    for a in f.added:
        gds[a["clone"]] = gds.get(a["clone"], 0) + 1
        if a.get("relaxed"):
            forced.add(a["clone"])
    txt = ",".join(f"{c}={n}{'*' if c in forced else ''}" for c, n in sorted(gds.items()) if n) or "none"
    return txt + ("" if not f.remaining and not f.over else " (incomplete)")


def _prior_text(p: dict | None) -> str | None:
    """e.g. `applied:FR1 0>1 (0.010; clear majority)` or `not_simpler:...` (strict targets kept)."""
    if not p:
        return None
    moves = ",".join(f"{region_label(s)} {o}>{n} ({d:.3f}; {k} majority)" for s, (o, n, d, k) in p["moves"].items())
    return f"{'applied' if p['applied'] else 'not_simpler'}:{moves}"


def _gd_vector(te: list[dict], cls: str, tree: Tree, num_clonal: int) -> dict[str, int]:
    """GDs per clone as step 05 will write them."""
    if cls == "no_gd":
        return {}
    if cls == "clonal_gd_only":
        return {tree.mrca: num_clonal} if num_clonal else {}
    out: dict[str, int] = {}
    for e in te:
        c = e["placed_clone"]
        if e["placed_as"] != "superseded_by_fill" and c in tree.parent:
            out[c] = out.get(c, 0) + 1
    return out


def _tumour_quality(ev, tri, trees, cp_by_tumour, regions, summary, ratios, alpaca_cn,
                    tumours, h_mismatch, fill_info, ploidy_rule, prior_info=None,
                    params: P.Params | None = None) -> dict[str, tuple[float, str | None]]:
    params = params or P.Params()
    prior_info = prior_info or {}
    targets_all: dict[str, dict[str, int]] = {}
    if regions is not None:
        for t, sid, n in regions.select("tumour_id", "sample_id", "num_gds").iter_rows():
            targets_all.setdefault(t, {})[sample_key(sid)] = int(n)
    fracs_all = _fracs_by_tumour(summary)
    maxab = {(r["tumour_id"], r["clone"]): max(r["ratio_A"], r["ratio_B"])
             for r in ratios.iter_rows(named=True) if r["ratio_A"] is not None and r["ratio_B"] is not None}
    med_total = {}
    if alpaca_cn is not None and {"pred_CN_A", "pred_CN_B"} <= set(alpaca_cn.columns):
        med_total = {(t, c): v for t, c, v in alpaca_cn.group_by("tumour_id", "clone")
                     .agg((pl.col("pred_CN_A") + pl.col("pred_CN_B")).median()).iter_rows()}
    num_clonal = dict(tumours.select("tumour_id", "num_clonal_gds").iter_rows()) if "num_clonal_gds" in tumours.columns else {}

    out = {}
    for t, cls in tri.items():
        tree = trees.get(t)
        te = [e for e in ev if e["tumour_id"] == t]
        factors: list[tuple[str, float, str]] = []
        targets = targets_all.get(t, {})
        fracs = fracs_all.get(t, {})
        prior = prior_info.get(t) if prior_info.get(t, {}).get("applied") else None
        if prior:
            # moved samples count once, through prior_raised / prior_lowered, not as borderline
            # or as flips that undo the prior
            targets = prior["targets"]
            fracs = {s: f for s, f in fracs.items() if s not in prior["moves"]}
            for s, (o, n, d, _) in prior["moves"].items():
                up = n > o
                factors.append(("prior_raised" if up else "prior_lowered",
                                Q.prior_factor(d, params.band, params.prior_raised_factor if up
                                               else params.prior_lowered_factor),
                                f"{region_label(s)}:{d:.3f}"))
            worst = {}
            for x in factors:
                if x[0] not in worst or x[1] < worst[x[0]][1]:
                    worst[x[0]] = x
            factors = list(worst.values())

        margins = Q.ploidy_margins(targets, fracs)
        # borderline-majority prior: the other samples at that boundary are borderline too; their
        # uncertainty is borderline_ploidy, so flipping them is not counted again as ambiguity
        flip_margins = dict(margins)
        if prior:
            bj = {max(o, n) for o, n, _, k in prior["moves"].values() if k == "borderline"}
            flip_margins = {s: m for s, m in margins.items()
                            if not any(fracs[s][j - 1] is not None and abs(fracs[s][j - 1] - 0.5) < params.band
                                       for j in bj)}
        if margins:
            s_min = min(margins, key=margins.get)
            factors.append(("borderline_ploidy", Q.graded(margins[s_min], Q.PLOIDY_MARGIN),
                            f"{margins[s_min]:.3f}:{region_label(s_min)}"))

        sup = [e["supporting_mut_count"] for e in te if e.get("is_subclonal_mutation_supported")
               and e.get("supporting_mut_count") is not None]
        if sup:
            factors.append(("low_mut_support", Q.graded(min(sup), Q.MIN_SUPPORT), str(min(sup))))

        flags = {f for e in te for f in e["flags"]}
        if "multicluster_parent_child" in flags:
            factors.append(("multicluster_parent_child", Q.FIXED["multicluster_parent_child"], ""))
        if "fill_tie" in flags:
            factors.append(("fill_tie", Q.FIXED["fill_tie"], ""))
        f_fill = fill_info.get(t)
        if isinstance(f_fill, Fill):
            covs = [(v, k) for a in f_fill.added if a.get("relaxed") for k, v in a["coverage"].items()]
            if covs:
                cmin, smin = min(covs)
                factors.append(("minority_coverage", Q.forced(cmin),
                                f"{cmin:.2f}:{region_label(smin)}"))
        if t in h_mismatch:
            factors.append(("ploidy_path_mismatch", Q.FIXED["ploidy_path_mismatch"], "%s vs %s" % h_mismatch[t]))
        n_unplaced = sum(e["placed_clone"] is None and e["placed_as"] != "superseded_by_fill" for e in te)
        if n_unplaced:
            factors.append(("unplaced_events", Q.FIXED["unplaced_events"], str(n_unplaced)))
        f = fill_info.get(t)
        if isinstance(f, Fill):
            n_pgdd = sum(e["placed_as"] == "superseded_by_fill" for e in te)
            # with the prior, a different count is the prior's doing and is already charged
            if len(f.added) != n_pgdd and not prior:
                factors.append(("fill_count_differs", Q.FIXED["fill_count_differs"], f"{len(f.added)} vs {n_pgdd}"))

        ambiguous_clones: set[str] = set()
        amb_detail = ""
        rev: list[str] = []
        gds: dict[str, int] = {}
        if tree is not None:
            gds = _gd_vector(te, cls, tree, int(num_clonal.get(t) or 0))
            r = [(maxab.get((t, c)), c) for c in gds if maxab.get((t, c)) is not None]
            if r:
                rmin, cmin = min(r)
                factors.append(("alpaca_disagrees", Q.graded(rmin, Q.ALPACA_THRESHOLD, 1.0), f"{rmin:.2f}:{cmin}"))
            # reversion: a clone strictly below a doubled clone that looks diploid again
            doubled = set(gds)
            rev = [c for c in tree.parent if any(a in doubled for a in tree.ancestors(c))
                   and med_total.get((t, c)) is not None and med_total[(t, c)] <= 2]
            if rev:
                factors.append(("reversion", Q.FIXED["reversion"], ",".join(rev)))
            # stability: does the fill answer change under small perturbations?
            g = cp_by_tumour.get(t)
            if targets and g is not None and ploidy_rule == "fill":
                props: dict[str, dict[str, float]] = {k: {} for k in targets}
                for k, c, v in g.select("key", "clone", "proportion").iter_rows():
                    if k in props and v > 0:
                        props[k][c] = v
                if all(props.values()):
                    fixed: dict[str, int] = {}
                    for e in te:
                        if e.get("is_subclonal_mutation_supported") and e["placed_clone"] in tree.parent:
                            fixed[e["placed_clone"]] = fixed.get(e["placed_clone"], 0) + 1
                    sc = {c: maxab.get((t, c)) for c in tree.parent}

                    def solve(tg, maj=MAJORITY):
                        res = fill_ploidy(tree, props, tg, fixed, sc, maj)
                        if res.remaining or res.over:
                            return None
                        v = dict(fixed)
                        for a in res.added:
                            v[a["clone"]] = v.get(a["clone"], 0) + 1
                        return v

                    base = solve(targets)
                    alts = [solve(tg) for tg in Q.flips(targets, fracs, flip_margins)] if fracs else []
                    alts += [solve(targets, m) for m in params.ambiguity_majorities]
                    diff = []
                    for a in alts:
                        if a is not None and base is not None and a != base and a not in diff:
                            diff.append(a)
                    for a in diff:
                        ambiguous_clones |= {c for c in set(a) | set(base) | set(gds) if a.get(c, 0) != base.get(c, 0)}
                    if diff:
                        txt = "|".join(",".join(f"{c}={n}" for c, n in sorted(a.items())) or "none" for a in diff[:3])
                        amb_detail = f"alt:{txt}"
                        factors.append(("ambiguous_solution", Q.FIXED["ambiguous_solution"], amb_detail))
        # per-event score: the factors that belong to one GD (its clone, samples, support)
        nth: dict[str, int] = {}
        for e in te:
            c = e["placed_clone"]
            if e["placed_as"] == "superseded_by_fill" or tree is None or c not in tree.parent:
                e["event_quality"], e["event_quality_code"] = None, None
                continue
            ef: list[tuple[str, float, str]] = []
            ev_samples = [sample_key(x) for x in (e.get("samples") or "").split(",") if x] or list(targets)
            # borderline_ploidy at this GD's own level (2026-11-12, LTX0693: two clear clonal GDs and a
            # borderline third): GDs above its clone + its rank on the clone (fixed GDs first, then fill)
            nth[c] = nth.get(c, 0) + 1
            level = sum(gds.get(a, 0) for a in tree.ancestors(c)) + nth[c]
            m = {x: abs(fracs[x][level - 1] - 0.5) for x in ev_samples
                 if level <= 3 and x in fracs and fracs[x][level - 1] is not None}
            if m:
                s_min = min(m, key=m.get)
                ef.append(("borderline_ploidy", Q.graded(m[s_min], Q.PLOIDY_MARGIN), f"{m[s_min]:.3f}:{region_label(s_min)}"))
            if e.get("is_subclonal_mutation_supported") and e.get("supporting_mut_count") is not None:
                ef.append(("low_mut_support", Q.graded(e["supporting_mut_count"], Q.MIN_SUPPORT), str(e["supporting_mut_count"])))
            r = maxab.get((t, c))
            if r is not None:
                ef.append(("alpaca_disagrees", Q.graded(r, Q.ALPACA_THRESHOLD, 1.0), f"{r:.2f}:{c}"))
            for fl in e["flags"]:
                if fl.startswith("minority_coverage("):
                    cov, reg = min((float(x.split(":")[1]), x.split(":")[0])
                                   for x in fl[len("minority_coverage("):-1].split(","))
                    ef.append(("minority_coverage", Q.forced(cov), f"{cov:.2f}:{reg}"))
            if "multicluster_parent_child" in e["flags"]:
                ef.append(("multicluster_parent_child", Q.FIXED["multicluster_parent_child"], ""))
            if "fill_tie" in e["flags"]:
                ef.append(("fill_tie", Q.FIXED["fill_tie"], e.get("fill_alternatives") or ""))
            for fl in e["flags"]:
                if fl.startswith("prior_raised("):
                    d, reg = max((float(x.split(":")[1]), x.split(":")[0])
                                 for x in fl[len("prior_raised("):-1].split(","))
                    ef.append(("prior_raised", Q.prior_factor(d, params.band, params.prior_raised_factor),
                               f"{reg}:{d:.3f}"))
            if c in ambiguous_clones:
                ef.append(("ambiguous_solution", Q.FIXED["ambiguous_solution"], amb_detail))
            # reversion: charged to the nearest doubled ancestor of the reverted clone
            mine = [x for x in rev if next((a for a in tree.ancestors(x) if a in gds), None) == c]
            if mine:
                ef.append(("reversion", Q.FIXED["reversion"], ",".join(mine)))
            e["event_quality"], e["event_quality_code"] = Q.score(ef)
            e["_factors"] = [list(x) for x in ef if x[1] < 1]

        # tumour score = per code, worst over its events and the tumour-level factors (Q.aggregate),
        # so it can be recomputed after rejecting events (step 05 --reject-below)
        tumour_only = [list(x) for x in factors if x[0] in Q.TUMOUR_ONLY and x[1] < 1]
        amb = [x for x in factors if x[0] == "ambiguous_solution"]
        if amb and not any(f[0] == "ambiguous_solution" for e in te for f in e.get("_factors", [])):
            tumour_only.append(list(amb[0]))  # e.g. an alternative adds a GD where there is none
        raised = [x for x in factors if x[0] == "prior_raised"]
        if raised and not any(f[0] == "prior_raised" for e in te for f in e.get("_factors", [])):
            tumour_only.append(list(raised[0]))  # raised sample doubled by a fixed (clonal / mutation) GD
        q, code = Q.aggregate([e.get("_factors", []) for e in te], tumour_only)
        out[t] = (q, code, json.dumps(tumour_only))
    for e in ev:
        e["event_factors"] = json.dumps(e.pop("_factors")) if "_factors" in e else None
        e.setdefault("event_quality", None)
        e.setdefault("event_quality_code", None)
    return out
