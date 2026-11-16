"""Synthetic tumours for the WGD rules deck, run through the real pipeline code (WGD/wgd).

Each case builds the step 04 inputs (PGDD events, tumour row, per-sample num_gds, ALPACA tree
and ratios, clone proportions) and calls `placement.place_events` and `final_table.build_table`,
exactly as steps 04 and 05 do. Every number on the slides comes from here. Ids are synthetic.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

WGD = Path(__file__).resolve().parents[2] / "WGD"
sys.path.insert(0, str(WGD))

from wgd import final_table, placement  # noqa: E402
from wgd import params as P  # noqa: E402

T = "SYN0001-Tumour1"


def sample(s: str) -> str:
    return f"SYN0001_SU_T1-{s}--abc"


@dataclass
class Case:
    key: str
    parent: dict[str, str]                          # clone -> parent ("diploid" for the MRCA)
    cp: dict[str, dict[str, float]]                 # sample -> {clone: proportion}
    num_gds: dict[str, int]                         # sample -> PGDD num_gds (target)
    events: list[tuple]                             # (is_clonal, mut_supported, cluster, samples)
    first_gd: str = "Subclonal"
    ratio: dict[str, float] = field(default_factory=dict)   # max(A, B) vs parent; default 1.0
    frac: dict[str, tuple] | None = None            # sample -> genome_frac_nMaj2/3/4 (needed by the prior)
    mut: dict[str, int] = field(default_factory=dict)       # cluster -> supporting doubled mutations
    rule: str = "fill"
    params: P.Params | None = None
    res: object = None
    table: pl.DataFrame | None = None

    def inputs(self):
        p = self.parent
        rt = pl.DataFrame({"tumour_id": [T] * len(p), "clone": list(p), "parent": list(p.values()),
                           "ratio_A": [self.ratio.get(c, 1.0) for c in p], "ratio_B": [1.0] * len(p)})
        cpdf = pl.DataFrame([{"tumour_id": T, "clone": c, "sample_id": sample(s),
                              "proportion": self.cp[s].get(c, 0.0)} for s in self.cp for c in p])
        ev = pl.DataFrame([{"GD_event": f"GD{i}", "GD_event_id": f"{T}_GD{i}", "tumour_id": T, "is_clonal": c,
                            "is_subclonal_mutation_supported": m, "clusters": cl,
                            "samples": ",".join(sample(s) for s in ss)}
                           for i, (c, m, cl, ss) in enumerate(self.events, 1)], schema_overrides={"clusters": pl.Utf8})
        n_clonal = sum(1 for e in self.events if e[0])
        tu = pl.DataFrame({"tumour_id": [T], "First_GD": [self.first_gd],
                           "num_subclonal_gds": [len(self.events) - n_clonal],
                           "num_second_gd": [int(max(self.num_gds.values()) > 1)], "num_clonal_gds": [n_clonal]})
        mc = pl.DataFrame([{"tumour_id": T, "cluster_id": k, "num_cn2": v, "discover_pass_all": True,
                            "check_pass_all": False} for k, v in self.mut.items()],
                          schema={"tumour_id": pl.Utf8, "cluster_id": pl.Utf8, "num_cn2": pl.Int64,
                                  "discover_pass_all": pl.Boolean, "check_pass_all": pl.Boolean})
        reg = pl.DataFrame({"tumour_id": [T] * len(self.num_gds), "sample_id": [sample(s) for s in self.num_gds],
                            "num_gds": list(self.num_gds.values())})
        summ = None
        if self.frac:
            summ = pl.DataFrame({"tumour_id": [T] * len(self.frac), "sample_name_hash": [sample(s) for s in self.frac],
                                 "genome_frac_nMaj2": [v[0] for v in self.frac.values()],
                                 "genome_frac_nMaj3": [v[1] for v in self.frac.values()],
                                 "genome_frac_nMaj4": [v[2] for v in self.frac.values()],
                                 "num_gds": [self.num_gds[s] for s in self.frac]})
        return ev, tu, mc, rt, cpdf, reg, summ

    def run(self) -> "Case":
        ev, tu, mc, rt, cpdf, reg, summ = self.inputs()
        self.res = placement.place_events(ev, tu, mc, rt, cpdf, summary=summ, score="max_AB", ploidy_rule=self.rule,
                                          regions=reg, params=self.params)
        self.table = final_table.build_table(self.res.events, self.res.tumours, tu, rt)
        return self

    # ---- views
    def gds(self) -> dict[str, int]:
        """GDs per clone, as in wgd_calls.tsv (notebook rule: per placed event)."""
        if self.rule == "notebook":
            out: dict[str, int] = {}
            for lb in self.res.events["placed_as"].to_list():
                c = lb.split("_", 1)[1] if lb and "_" in lb else None
                if c in self.parent:
                    out[c] = out.get(c, 0) + 1
            return out
        return {c: n for c, n in zip(self.table["clone"], self.table["number_of_gds_at_clone"]) if n}

    def order(self) -> list[str]:
        """Clones in the order the GDs were placed (notebook: per event)."""
        return [lb.split("_", 1)[1] for lb in self.res.events["placed_as"].to_list() if lb and "_" in lb
                and lb.split("_", 1)[1] in self.parent]

    @property
    def tumour(self) -> dict:
        return self.res.tumours.row(0, named=True)

    def levels(self, gds: dict[str, int] | None = None) -> dict[str, int]:
        props = {s: dict(v) for s, v in self.cp.items()}
        return placement.sample_levels(placement.Tree(self.parent), props, gds if gds is not None else self.gds())


# ---- the cases ---------------------------------------------------------------------------------

# Nested first + second ploidy-only GD (doc section 5, item 1): GD1 on clone2, GD2 on clone4 (title slide).
NESTED = dict(parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone1", "clone4": "clone2",
                      "clone5": "clone2"},
              cp={"R1": {"clone4": 1.0}, "R2": {"clone5": 1.0}, "R3": {"clone3": 1.0}},
              num_gds={"R1": 2, "R2": 1, "R3": 0},
              events=[(False, False, None, ["R1", "R2"]), (False, False, None, ["R1"])],
              ratio={"clone2": 1.95, "clone4": 2.05})


# Clonal GD + two more GDs in R3 only. The old rule can't reuse a clone, so the second GD goes on
# clone2, which also doubles R2 past its single GD; both ratios > 1.5, so it still says resolved.
STACKED = dict(parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone2", "clone5": "clone3",
                       "clone6": "clone2", "clone4": "clone6", "clone7": "clone1"},
               cp={"R1": {"clone7": 1.0}, "R2": {"clone4": 0.62, "clone6": 0.38}, "R3": {"clone3": 0.6, "clone5": 0.4}},
               num_gds={"R1": 1, "R2": 1, "R3": 3},
               events=[(True, False, None, ["R1", "R2", "R3"]), (False, False, None, ["R3"]),
                       (False, False, None, ["R3"])],
               first_gd="Clonal", ratio={"clone1": 1.9, "clone3": 1.6, "clone2": 1.55})


def stacked_notebook() -> Case:
    return Case("stacked_notebook", rule="notebook", **STACKED).run()


def stacked_fill() -> Case:
    return Case("stacked_fill", **STACKED).run()


def nested_fill() -> Case:
    return Case("nested_fill", **NESTED).run()


# Majority ploidy prior: R1-R3 clearly doubled, R4 just under 0.5. One ploidy-only GD in R1-R3.
PRIOR = dict(parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone1", "clone4": "clone1"},
             cp={"R1": {"clone2": 1.0}, "R2": {"clone3": 1.0}, "R3": {"clone4": 1.0},
                 "R4": {"clone1": 0.4, "clone2": 0.3, "clone3": 0.3}},
             num_gds={"R1": 1, "R2": 1, "R3": 1, "R4": 0},
             events=[(False, False, None, ["R1", "R2", "R3"])],
             ratio={"clone1": 1.9},
             frac={"R1": (0.96, 0.1, 0.0), "R2": (0.95, 0.1, 0.0), "R3": (0.97, 0.1, 0.0), "R4": (0.48, 0.05, 0.0)})


def prior_strict() -> Case:
    return Case("prior_strict", params=P.Params(prior_enabled=False), **PRIOR).run()


def prior_applied() -> Case:
    return Case("prior_applied", **PRIOR).run()


# Relaxed pass: R1 is doubled but its doubled clone is only 40% of it; the clone that makes up the
# other 60% also fills the undoubled R2, so no clone is a clean majority.
RELAXED = dict(parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone1", "clone4": "clone2"},
               cp={"R1": {"clone4": 0.4, "clone2": 0.6}, "R2": {"clone2": 1.0}, "R3": {"clone3": 1.0}},
               num_gds={"R1": 1, "R2": 0, "R3": 0},
               events=[(False, False, None, ["R1"])], ratio={"clone4": 2.1})


def relaxed() -> Case:
    return Case("relaxed", **RELAXED).run()


# Mutation-supported GD on clone2 whose subtree also holds clone5, all of sample R3, which PGDD
# calls undoubled: a clone can't lose a GD, so the tumour fails (reversion detector).
REVERSION = dict(parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone1", "clone4": "clone2",
                         "clone5": "clone2"},
                 cp={"R1": {"clone4": 1.0}, "R2": {"clone3": 1.0}, "R3": {"clone5": 1.0}},
                 num_gds={"R1": 1, "R2": 0, "R3": 0},
                 events=[(False, True, "2", ["R1"])], mut={"2": 25}, ratio={"clone2": 2.0})


def reversion() -> Case:
    return Case("reversion", **REVERSION).run()


# Fixed GDs: a clonal GD and a mutation-supported GD on a parent-child cluster pair (-> the child).
FIXED = dict(parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone1", "clone4": "clone2"},
             cp={"R1": {"clone4": 1.0}, "R2": {"clone3": 1.0}},
             num_gds={"R1": 2, "R2": 1},
             events=[(True, False, None, ["R1", "R2"]), (False, True, "2,4", ["R1"])],
             first_gd="Clonal", mut={"4": 18, "2": 6}, ratio={"clone1": 2.0, "clone4": 1.9})


def fixed() -> Case:
    return Case("fixed", **FIXED).run()


def all_cases() -> dict[str, Case]:
    cs = [nested_fill(), stacked_notebook(), stacked_fill(), prior_strict(), prior_applied(), relaxed(), reversion(), fixed()]
    return {c.key: c for c in cs}


if __name__ == "__main__":
    for k, c in all_cases().items():
        t = c.tumour
        print(f"{k:16s} gds={c.gds()} order={c.order()} status={t['status']} fill={t.get('fill')} "
              f"q={t['quality']:.3f} code={t['quality_code']} prior={t.get('ploidy_prior')}")
