"""Synthetic cases for the deck, run through the real correction (`apc.correct`).

Every number on the slides comes from here: the cases are fed to `correct.run`, and the cost
tables use the same objective as `correct.best_cn` (checked in tests). Tumour ids are synthetic.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

APC = Path(__file__).resolve().parents[2] / "alpaca_post_hoc_correction"
sys.path.insert(0, str(APC))

from apc import correct  # noqa: E402

SEGS = ["1_0_100000000", "2_0_50000000", "3_0_50000000", "5_0_100000000"]
SEG_LABEL = {s: f"chr{s.split('_')[0]}" for s in SEGS}


@dataclass
class Case:
    key: str
    title: str
    parent: dict[str, str]
    gds: dict[str, int]
    present: dict[str, float]                 # clone -> proportion in the (single) sample
    cn: dict[str, list[tuple[int, int]]]      # clone -> (A, B) per segment in SEGS
    res: dict[str, pl.DataFrame] = field(default_factory=dict)

    @property
    def tumour(self) -> str:
        return f"SYN-{self.key}"

    def run(self) -> "Case":
        a = pl.DataFrame([{"tumour_id": self.tumour, "segment": s, "clone": c, "pred_CN_A": x, "pred_CN_B": y}
                          for c, v in self.cn.items() for s, (x, y) in zip(SEGS, v)])
        p = pl.DataFrame({"tumour_id": [self.tumour] * len(self.present), "clone": list(self.present),
                          "sample_id": ["R1"] * len(self.present), "proportion": list(self.present.values())})
        w = pl.DataFrame({"tumour_id": [self.tumour] * len(self.parent), "clone": list(self.parent),
                          "parent": list(self.parent.values()),
                          "number_of_gds_at_clone": [self.gds.get(c, 0) for c in self.parent],
                          "quality": [0.9] * len(self.parent)})
        self.res = correct.run(a, p, w)
        return self

    # ---- views used by the slides
    def corrected_clones(self) -> list[str]:
        """Clones whose CN changed, in the order the correction processed them (deepest first)."""
        log = self.res["corrections"]
        if log.is_empty():
            return []
        return [r["clone"] for r in log.iter_rows(named=True) if r["allele_segments_changed"] > 0]

    def candidates(self) -> list[str]:
        log = self.res["corrections"]
        return [] if log.is_empty() else log["clone"].to_list()

    def cn_after(self, clone: str) -> list[tuple[int, int]]:
        g = self.res["alpaca_corrected"].filter(pl.col("clone") == clone)
        d = {s: (a, b) for s, a, b in g.select("segment", "pred_CN_A", "pred_CN_B").iter_rows()}
        return [d[s] for s in SEGS]

    def ratio(self, clone: str, corrected: bool) -> float:
        """max(ratio_A, ratio_B): the edge statistic the WGD pipeline and the flags use."""
        r = self.res["ratio_comparison"].filter(pl.col("clone") == clone).row(0, named=True)
        sfx = "_corrected" if corrected else ""
        return max(r[f"ratio_A{sfx}"], r[f"ratio_B{sfx}"])

    def flags(self) -> list[tuple[str, str]]:
        return [(r["clone"], r["flag"]) for r in self.res["flags"].iter_rows(named=True)]


def cost_row(child_cn: list[int], mult: list[int], parent_cn: int, xs: range) -> list[int]:
    """sum_c |CN_c - m_c x| + |x - CN_parent| for each candidate x (same objective as best_cn)."""
    return [int(sum(abs(c - m * x) for c, m in zip(child_cn, mult)) + abs(x - parent_cn)) for x in xs]


def cost_table(case: Case, clone: str, allele: int, xs: range = range(0, 5),
               child_cn: dict[str, list[tuple[int, int]]] | None = None) -> list[dict]:
    """One row per segment for `clone`, allele 0 (A) or 1 (B), with ALPACA's CN as the inputs.
    `child_cn` overrides a child's per-segment CN (e.g. a child the correction already changed)."""
    kids = correct.children(case.parent, clone)
    par = case.parent[clone]
    cn = {**case.cn, **(child_cn or {})}
    rows = []
    for i, s in enumerate(SEGS):
        child = [cn[k][i][allele] for k in kids]
        mult = [2 ** case.gds.get(k, 0) for k in kids]
        p = 1 if par == "diploid" else case.cn[par][i][allele]
        costs = cost_row(child, mult, p, xs)
        orig = case.cn[clone][i][allele]
        best = correct.best_cn(child, mult, p, orig)
        rows.append({"segment": SEG_LABEL[s], "children": dict(zip(kids, child)), "mult": dict(zip(kids, mult)),
                     "parent": p, "costs": costs, "alpaca": orig, "new": best})
    return rows


# ---- the cases -------------------------------------------------------------------------------

def parallel() -> Case:
    """The motivating pattern: MRCA absent from every sample, two children doubled independently.
    Truth: clone1 = 1|1 (2|1 on chr5); each child = 2 x clone1, then one private change."""
    return Case("A", "Parallel GDs under an absent MRCA",
                parent={"clone1": "diploid", "clone5": "clone1", "clone6": "clone1"},
                gds={"clone5": 1, "clone6": 1},
                present={"clone5": 0.6, "clone6": 0.4},
                cn={"clone1": [(2, 2), (2, 2), (2, 2), (4, 2)],
                    "clone5": [(2, 2), (2, 2), (3, 2), (4, 2)],
                    "clone6": [(2, 2), (2, 1), (2, 2), (4, 2)]}).run()


def mixed() -> Case:
    """Absent parent with one doubled child and one undoubled child: the undoubled child vouches for
    ALPACA's value, so the rule fires but the minimiser keeps 2|2."""
    return Case("B", "One doubled child, one not",
                parent={"clone1": "diploid", "clone5": "clone1", "clone7": "clone1"},
                gds={"clone5": 1},
                present={"clone5": 0.5, "clone7": 0.5},
                cn={"clone1": [(2, 2)] * 4, "clone5": [(4, 4)] * 4, "clone7": [(2, 2)] * 4}).run()


def chain(clone1_present: bool = False) -> Case:
    """clone1 -> clone2 -> {clone3 GD, clone4 GD}; clone1 and clone2 absent (or clone1 observed)."""
    present = {"clone3": 0.5, "clone4": 0.5}
    if clone1_present:
        present = {"clone1": 0.1, "clone3": 0.5, "clone4": 0.4}
    return Case("C2" if clone1_present else "C", "Chain of missing clones" + (" (clone1 observed)" if clone1_present else ""),
                parent={"clone1": "diploid", "clone2": "clone1", "clone3": "clone2", "clone4": "clone2"},
                gds={"clone3": 1, "clone4": 1},
                present=present,
                cn={c: [(2, 2)] * 4 for c in ("clone1", "clone2", "clone3", "clone4")}).run()


def parent_observed() -> Case:
    c = parallel()
    return Case("D1", "Parent has cells", c.parent, c.gds, {"clone1": 0.1, "clone5": 0.5, "clone6": 0.4}, c.cn).run()


def no_gd() -> Case:
    c = parallel()
    return Case("D2", "No GD on any child", c.parent, {}, c.present, c.cn).run()


def parent_gd() -> Case:
    c = parallel()
    return Case("D3", "The parent itself doubled", c.parent, {"clone1": 1}, c.present, c.cn).run()


def reversion() -> Case:
    """Flag only: clone2 sits below a doubled clone1 but ALPACA has it back at 1|1."""
    return Case("E", "Reversion below a GD",
                parent={"clone1": "diploid", "clone2": "clone1"}, gds={"clone1": 1},
                present={"clone1": 0.5, "clone2": 0.5},
                cn={"clone1": [(2, 2)] * 4, "clone2": [(1, 1)] * 4}).run()


def all_cases() -> dict[str, Case]:
    cs = [parallel(), mixed(), chain(), chain(True), parent_observed(), no_gd(), parent_gd(), reversion()]
    return {c.key: c for c in cs}


if __name__ == "__main__":
    for k, c in all_cases().items():
        print(k, c.title, "| candidates:", c.candidates(), "| changed:", c.corrected_clones(), "| flags:", c.flags())
    a = parallel()
    for al in (0, 1):
        for r in cost_table(a, "clone1", al):
            print("AB"[al], r)
    print({k: (round(a.ratio(k, False), 2), round(a.ratio(k, True), 2)) for k in a.parent})
    print(a.res["ratio_comparison"])
