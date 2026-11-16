"""Build alpaca_correction.pptx: how and when the post-hoc ALPACA correction is applied.

All numbers come from `examples.py`, which runs the synthetic cases through the real
`apc.correct` (../../alpaca_post_hoc_correction). Native, editable shapes, tables and charts.

    pixi run deck        # writes alpaca_correction.pptx next to this file
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_LEGEND_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

import examples as ex

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from deckkit import (AMBER, BODY, CORAL, CW, GRID, H, HEAD, INK, MUTED, MX, PALE_CORAL, TEAL, TINT,  # noqa: E402
                     WHITE, W, box, fill_text, line, rgb, table, text)
from deckkit import new_slide as _new_slide  # noqa: E402

FOOTER = "Synthetic data only · darwins_analysis/decks/alpaca_correction"


def new_slide(prs, title, kicker=None, dark=False, notes=""):
    return _new_slide(prs, title, kicker, dark, notes, footer=FOOTER)


# ------------------------------------------------------------------ tree diagram
def layout_tree(parent: dict[str, str], root_label: bool):
    kids = {}
    for c, p in parent.items():
        kids.setdefault(p, []).append(c)
    for v in kids.values():
        v.sort(key=lambda c: int(c.replace("clone", "")) if c.replace("clone", "").isdigit() else c)
    roots = kids.get("diploid", [])
    order, depth = [], {}

    def walk(c, d):
        depth[c] = d
        order.append(c)
        for k in kids.get(c, []):
            walk(k, d + 1)

    if root_label:
        walk("diploid", 0)
    else:
        for r in roots:
            walk(r, 0)
    leaves = [c for c in order if not kids.get(c)]
    xs = {c: float(i) for i, c in enumerate(leaves)}

    def place(c):
        if c in xs:
            return xs[c]
        xs[c] = sum(place(k) for k in kids[c]) / len(kids[c])
        return xs[c]

    for c in order:
        place(c)
    return order, depth, xs, kids, len(leaves)


def tree(slide, case: ex.Case, x, y, w, h, *, state=None, cn=None, pill=(1.05, 0.4), font=13, cn_font=12,
         root_label=True, badges=None, cn_w=1.5):
    """state: clone -> observed | absent | corrected (default from case.present / corrected clones).
    cn: clone -> label under the pill (str or list of runs). badges: clone -> step number."""
    state = state or {}
    order, depth, xs, kids, n_leaves = layout_tree(case.parent, root_label)
    levels = max(depth.values()) + 1
    pw, ph = pill
    gap_y = (h - ph) / max(levels - 1, 1)
    span = w - pw
    pos = {}
    for c in order:
        cx = x + pw / 2 + (span * xs[c] / max(n_leaves - 1, 1) if n_leaves > 1 else span / 2)
        cy = y + depth[c] * gap_y
        pos[c] = (cx, cy)
    for c in order:
        for k in kids.get(c, []):
            (px, py), (kx, ky) = pos[c], pos[k]
            line(slide, px, py + ph, kx, ky, color=GRID, lw=1.75)
    changed = set(case.corrected_clones())
    for c in order:
        cx, cy = pos[c]
        if c == "diploid":
            b = box(slide, cx - pw / 2, cy, pw, ph, fill=TINT, line=GRID, shape=MSO_SHAPE.ROUNDED_RECTANGLE,
                    radius=0.5, name="node diploid")
            fill_text(b, "diploid", size=font - 1, color=MUTED, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE,
                      margin=0)
            continue
        st = state.get(c) or ("observed" if case.present.get(c) else ("corrected" if c in changed else "absent"))
        if st == "observed":
            b = box(slide, cx - pw / 2, cy, pw, ph, fill=TEAL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5,
                    name=f"node {c}")
            tc = WHITE
        elif st == "corrected":
            b = box(slide, cx - pw / 2, cy, pw, ph, fill=PALE_CORAL, line=CORAL, lw=2.0, dash=True,
                    shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5, name=f"node {c}")
            tc = CORAL
        else:
            b = box(slide, cx - pw / 2, cy, pw, ph, fill=WHITE, line=TEAL, lw=1.75, dash=True,
                    shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5, name=f"node {c}")
            tc = TEAL
        fill_text(b, c, size=font, bold=True, color=tc, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, margin=0)
        if case.gds.get(c):
            d = 0.36
            g = box(slide, cx + pw / 2 - d * 0.55, cy - d * 0.45, d, d, fill=AMBER, shape=MSO_SHAPE.OVAL,
                    name=f"GD {c}")
            fill_text(g, "GD", size=10, bold=True, color=INK, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE,
                      margin=0)
        if badges and c in badges:
            d = 0.36
            g = box(slide, cx - pw / 2 - d * 0.45, cy - d * 0.45, d, d, fill=INK, shape=MSO_SHAPE.OVAL,
                    name=f"step {c}")
            fill_text(g, str(badges[c]), size=11, bold=True, color=WHITE, align=PP_ALIGN.CENTER,
                      anchor=MSO_ANCHOR.MIDDLE, margin=0)
        if cn and c in cn:
            col = CORAL if st == "corrected" else MUTED
            if kids.get(c):          # internal node: label beside the pill, clear of the edges
                dx = 0.28 if case.gds.get(c) else 0.1      # clear of the GD badge
                text(slide, cx + pw / 2 + dx, cy, cn_w, ph, cn[c], size=cn_font, color=col, margin=0,
                     anchor=MSO_ANCHOR.MIDDLE)
            else:
                text(slide, cx - pw / 2 - 0.45, cy + ph + 0.02, pw + 0.9, 0.32, cn[c], size=cn_font, color=col,
                     align=PP_ALIGN.CENTER, margin=0)
    return pos


def legend(slide, x, y):
    items = [("observed (cells in a sample)", dict(fill=TEAL)),
             ("absent from every sample", dict(fill=WHITE, line=TEAL, lw=1.75, dash=True)),
             ("absent, CN corrected", dict(fill=PALE_CORAL, line=CORAL, lw=2.0, dash=True))]
    cx = x
    for lab, st in items:
        box(slide, cx, y + 0.04, 0.42, 0.22, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5, **st)
        text(slide, cx + 0.5, y, 2.6, 0.3, lab, size=12, color=MUTED, margin=0)
        cx += 3.05
    box(slide, cx, y + 0.02, 0.26, 0.26, fill=AMBER, shape=MSO_SHAPE.OVAL)
    text(slide, cx + 0.34, y, 2.6, 0.3, "GD called on this clone", size=12, color=MUTED, margin=0)


def ab(t):
    return f"{t[0]}|{t[1]}"


# ------------------------------------------------------------------ tables / strips
def cost_strip(slide, x, y, costs, chosen, alpaca, *, cell=0.78, label=None):
    """Boxes for x = 0..n with the cost inside; the minimiser in coral; ALPACA's value marked."""
    if label:
        text(slide, x, y - 0.02, 2.2, 0.62, label, size=14, color=INK, margin=0, anchor=MSO_ANCHOR.MIDDLE)
    x0 = x + (2.3 if label else 0)
    for i, c in enumerate(costs):
        bx = x0 + i * (cell + 0.12)
        text(slide, bx, y - 0.34, cell, 0.3, f"x = {i}", size=12, color=MUTED, align=PP_ALIGN.CENTER, margin=0)
        best = i == chosen
        b = box(slide, bx, y, cell, 0.58, fill=CORAL if best else TINT, line=None,
                shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.2)
        fill_text(b, str(c), size=18, bold=best, color=WHITE if best else INK, align=PP_ALIGN.CENTER,
                  anchor=MSO_ANCHOR.MIDDLE, margin=0)
        if i == alpaca:
            text(slide, bx - 0.15, y + 0.6, cell + 0.3, 0.3, "ALPACA", size=11, bold=True, color=TEAL,
                 align=PP_ALIGN.CENTER, margin=0)
    return x0 + len(costs) * (cell + 0.12)


# ------------------------------------------------------------------ slides
def s_title(prs, A):
    s = new_slide(prs, "Post-hoc ALPACA correction", dark=True, notes=(
        "What this deck covers: why ALPACA's copy number needs a correction after parallel genome doublings, "
        "exactly which clones the correction touches (and which it never touches), how the new copy number is "
        "chosen, and worked synthetic examples. Every number on the slides comes from running the real code "
        "(darwins_analysis alpaca_post_hoc_correction, apc/correct.py) on synthetic trees."))
    s.shapes.title.top = Inches(1.8)
    s.shapes.title.width = Inches(7.4)
    s.shapes.title.height = Inches(1.7)
    s.shapes.title.text_frame.paragraphs[0].runs[0].font.size = Pt(44)
    text(s, MX, 3.6, 7.0, 1.2, "How and when ALPACA's copy number is re-inferred for absent ancestors "
         "after parallel genome doublings", size=20, color="C5D3D8", margin=0)
    text(s, MX, 5.2, 7.0, 0.9, ["darwins_analysis · alpaca_post_hoc_correction (PR #3)",
                                "Worked on synthetic trees · 2026-11-16"], size=14, color="9FB0B7", margin=0)
    card = box(s, 8.35, 1.7, 4.35, 4.3, fill="243742", shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06)
    card.name = "tree card"
    tree(s, A, 8.6, 2.1, 3.3, 2.9, pill=(1.1, 0.42), cn_w=1.25,
         cn={c: f"{ab(A.cn[c][0])} → {ab(A.cn_after(c)[0])}" if c in A.corrected_clones() else ab(A.cn[c][0])
             for c in A.parent})
    text(s, 8.65, 5.35, 3.75, 0.5, "chr1 copy number, ALPACA → corrected", size=12, color="C5D3D8",
         align=PP_ALIGN.CENTER, margin=0)


def s_problem(prs, A):
    s = new_slide(prs, "ALPACA cannot see a genome doubling",
                  "Two parallel doublings look, to parsimony, like one gain on the absent parent", notes=(
        "ALPACA infers integer allele-specific copy number per clone by parsimony over the tree; it has no "
        "genome-doubling event. A clone with no cells in any sample is constrained only by its relatives. When "
        "two of its children doubled independently, one shared gain on the parent is cheaper than two parallel "
        "doublings, so ALPACA puts the doubling on the parent. Ratios are max(ratio_A, ratio_B) of each clone vs "
        "its parent, length-weighted over the four synthetic segments of example A."))
    truth = {c: v for c, v in zip(A.parent, ("1|1", "2|2", "2|2"))}
    alp = {c: ab(A.cn[c][0]) for c in A.parent}
    for i, (ttl, labels, st, ratios) in enumerate((
            ("What happened", truth, {"clone1": "absent"}, {c: A.ratio(c, True) for c in A.parent}),
            ("What ALPACA infers", alp, {"clone1": "absent"}, {c: A.ratio(c, False) for c in A.parent}))):
        x0 = MX + i * 4.15
        box(s, x0, 1.95, 3.85, 4.55, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
        text(s, x0 + 0.25, 2.07, 3.4, 0.4, ttl, size=18, bold=True, color=INK, margin=0)
        tree(s, A, x0 + 0.2, 2.7, 3.0, 2.45, state=st, pill=(1.05, 0.4), cn_font=12, cn_w=1.45,
             cn={c: [(labels[c], {"bold": True, "color": INK}), (f"   ratio {ratios[c]:.2f}", {})]
                 for c in A.parent})
        foot = ("clone1 1|1; clone5 and clone6 each doubled (GD) to 2|2" if i == 0
                else "one gain on clone1 (2|2); children unchanged")
        text(s, x0 + 0.25, 5.75, 3.4, 0.6, foot, size=13, color=MUTED, margin=0)
    rx = MX + 8.35
    text(s, rx, 1.95, CW - 8.35, 4.0, [
        [("Parsimony, no GD event. ", {"bold": True}), "ALPACA picks the fewest integer changes over the tree."],
        [("Absent clones float. ", {"bold": True}), "With no cells in any sample, only relatives constrain them."],
        [("Result: ", {"bold": True}), "the parent looks doubled (ratio ≈ 2) and the doubled children look "
                                      "unchanged (≈ 1), the opposite of the GD calls."],
    ], size=16, space_after=12, margin=0)
    b = box(s, rx, 5.4, CW - 8.35, 1.1, fill=PALE_CORAL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
    fill_text(b, "The correction re-infers the absent parent's CN once the WGD pipeline has called the GDs "
                 "on its children.", size=14, color=INK, anchor=MSO_ANCHOR.MIDDLE, margin=0.15)


def s_rule(prs):
    s = new_slide(prs, "When a clone is corrected",
                  "Two tumour-level gates, then three conditions on the clone, all required", notes=(
        "Tumour gates: the tumour must be in the WGD calls and have a tree; with --min-quality, tumours whose "
        "best call quality is below the cut-off are skipped (listed in skipped.tsv, CN left as ALPACA gave it). "
        "Clone conditions (apc/correct.py clones_to_correct): absent = no positive proportion in any sample of "
        "clone_proportions; no GD called on the clone itself; and a GD-carrying clone is reachable from it through "
        "absent, GD-free clones only. Piotr, 2026-11-13: a chain of missing clones above parallel GDs is corrected "
        "as a whole, and an observed clone stops the chain. Qualifying clones are processed deepest first."))
    gates = [("Tumour has WGD calls and a tree", "else: skipped.tsv"),
             ("Call quality ≥ --min-quality", "optional; else skipped.tsv")]
    y = 1.95
    text(s, MX, y, 3.6, 0.35, "TUMOUR", size=13, bold=True, color=MUTED, margin=0)
    for i, (a, b_) in enumerate(gates):
        bx = MX + i * 3.0
        g = box(s, bx, y + 0.42, 2.7, 1.15, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.1)
        fill_text(g, [[(a, {"bold": True})], [(b_, {"size": 12, "color": MUTED})]], size=14,
                  anchor=MSO_ANCHOR.MIDDLE, margin=0.12)
        if i == 0:
            line(s, bx + 2.72, y + 1.0, bx + 2.98, y + 1.0, color=MUTED, lw=1.75, arrow=True)
    y2 = 3.85
    text(s, MX, y2, 6.0, 0.35, "CLONE (all three)", size=13, bold=True, color=MUTED, margin=0)
    conds = [("Absent from every sample", "no cells in any region; observed clones are never changed"),
             ("No GD called on it", "a clone that doubled itself keeps ALPACA's CN"),
             ("Reaches a GD clone via absent, GD-free clones only",
              "a doubled child, or a chain of missing clones above one; an observed clone stops the chain")]
    out_w, gap = 1.9, 0.3
    cw = (CW - out_w - 3 * gap) / 3
    ch = 2.1
    for i, (a, b_) in enumerate(conds):
        bx = MX + i * (cw + gap)
        c = box(s, bx, y2 + 0.42, cw, ch, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
        fill_text(c, [[(a, {"bold": True})], [(b_, {"size": 13, "color": MUTED})]], size=15, margin=0.15,
                  anchor=MSO_ANCHOR.TOP)
        n = box(s, bx + 0.12, y2 + 0.24, 0.38, 0.38, fill=TEAL, shape=MSO_SHAPE.OVAL)
        fill_text(n, str(i + 1), size=13, bold=True, color=WHITE, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE,
                  margin=0)
        c.text_frame.margin_top = Inches(0.32)
    ox = MX + 3 * (cw + gap)
    line(s, ox - gap + 0.02, y2 + 0.42 + ch / 2, ox - 0.02, y2 + 0.42 + ch / 2, color=MUTED, lw=1.75, arrow=True)
    out = box(s, ox, y2 + 0.42, out_w, ch, fill=CORAL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
    fill_text(out, [[("Re-infer CN", {"bold": True})], [("deepest first", {"size": 13})]], size=16,
              color=WHITE, anchor=MSO_ANCHOR.MIDDLE, align=PP_ALIGN.CENTER, margin=0.08)
    text(s, MX + 6.3, 2.1, CW - 6.3, 1.6, [
        "Clones are processed bottom-up, so a corrected child feeds its parent's step.",
        "A qualifying clone can still keep ALPACA's value: the optimisation decides (example B)."],
        size=14, color=INK, bullets=True, space_after=8, margin=0)


def s_gallery(prs, cases):
    s = new_slide(prs, "Which clones qualify: six synthetic trees",
                  "Same rule, different outcomes; outcomes are the code's, not hand-drawn", notes=(
        "Each panel was run through correct.run. A: classic parallel GDs under an absent MRCA, clone1 corrected. "
        "C: chain of two missing clones above parallel GDs: clone2 then clone1. C2: same but clone1 has cells, so "
        "only clone2. D1: parent observed, nothing corrected; the GD calls and ALPACA disagree on clone5/clone6 and "
        "that is flagged (slide 10). D2: no GD on any child, nothing to do. D3: the parent carries the GD itself "
        "(a clonal GD), so it is not a candidate."))
    panels = [("A", "Parallel GDs, absent MRCA"), ("C", "Chain of missing clones"),
              ("C2", "Chain, clone1 observed"), ("D1", "Parent has cells"), ("D2", "No GD on any child"),
              ("D3", "The parent itself doubled")]
    pw, ph, gx, gy = 3.84, 2.3, 0.305, 0.2
    for i, (k, ttl) in enumerate(panels):
        c = cases[k]
        x0, y0 = MX + (i % 3) * (pw + gx), 1.75 + (i // 3) * (ph + gy)
        box(s, x0, y0, pw, ph, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
        text(s, x0 + 0.18, y0 + 0.1, pw - 0.36, 0.32, f"{k}  {ttl}", size=14, bold=True, margin=0)
        levels = max(ex.correct.depth(c.parent, x) for x in c.parent) + 1
        th = 1.3 if levels == 3 else 0.9
        tree(s, c, x0 + 0.6, y0 + 0.58 + (1.3 - th) / 2, pw - 1.2, th, pill=(0.9, 0.3), font=11,
             root_label=False)
        ch = c.corrected_clones()
        verdict = (f"corrected: {' then '.join(ch)}" if ch else "nothing corrected")
        text(s, x0 + 0.18, y0 + ph - 0.38, pw - 0.36, 0.3, verdict, size=13, bold=bool(ch),
             color=CORAL if ch else MUTED, margin=0)
    legend(s, MX, 6.66)


def s_how(prs):
    s = new_slide(prs, "How: a small search per segment and allele",
                  "The new copy number x best explains the children given their GDs, staying close to the parent",
                  notes=(
        "apc/correct.py best_cn. For each segment and each allele separately, try every integer x from 0 to "
        "max(parent, ALPACA, children) + 1 and keep the cheapest; x = 0 is excluded when any child still carries "
        "the allele, since a lost allele cannot come back (Piotr, 2026-11-16). Every child of the clone counts, doubled or not: "
        "a child with g GDs is expected at 2^g times x. The parent term uses the parent's current value (1 for "
        "diploid); in a chain the parent is still ALPACA's value at that point, since ancestors are processed "
        "after. Ties keep ALPACA's value if it is among the minimisers, else the smallest. Afterwards total, "
        "minor/major, LOH and homdel columns are recomputed and the WGD ratios are recomputed exactly as ALPACA's "
        "analysis.calculate_wgd_ratios (matches ALPACA to 1e-9 on all 34 tx_data mock clones)."))
    b = box(s, MX, 1.95, CW, 1.45, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
    fill_text(b, [("cost(x)  =  ", {"italic": False}), ("Σ", {"size": 34}), ("children c", {"size": 14, "sub": True}),
                  ("  | CN", {}), ("c", {"sub": True, "size": 18}), (" − 2", {}),
                  ("g", {"sup": True, "italic": True, "size": 18}), ("c", {"size": 12}), (" · x |   +   | x − CN", {}),
                  ("parent", {"sub": True, "size": 18}), (" |", {})],
              size=30, font=HEAD, italic=False, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    terms = [([("x", {})], "candidate CN, 0 to max + 1; ≥ 1 if a child has the allele"),
             ([("CN", {}), ("c", {"sub": True})], "child c's ALPACA CN, this segment and allele"),
             ([("g", {}), ("c", {"sub": True})], "GDs called on child c (0 if it did not double)"),
             ([("CN", {}), ("parent", {"sub": True})], "parent's current CN; 1 if diploid")]
    for i, (sym, d) in enumerate(terms):
        y = 3.75 + i * 0.62
        sb = box(s, MX, y, 1.45, 0.46, fill=INK, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.3)
        fill_text(sb, sym, size=15, bold=True, color=WHITE, font=HEAD, align=PP_ALIGN.CENTER,
                  anchor=MSO_ANCHOR.MIDDLE, margin=0)
        text(s, MX + 1.65, y, 5.2, 0.46, d, size=15, anchor=MSO_ANCHOR.MIDDLE, margin=0)
    rx = MX + 7.1
    text(s, rx, 3.7, CW - 7.1, 0.4, "Then", size=15, bold=True, color=MUTED, margin=0)
    text(s, rx, 4.1, CW - 7.1, 2.4, [
        "Keep the cheapest x. Ties: ALPACA's value if it is among them, else the smallest.",
        "Recompute total, minor/major, LOH and homdel columns.",
        "Recompute WGD ratios exactly as ALPACA does (port matches to 1e-9 on 34 mock clones).",
    ], size=15, bullets=True, space_after=8, margin=0)


def s_example_a(prs, A):
    s = new_slide(prs, "Example A: cost of each candidate x",
                  "clone1 absent; clone5 and clone6 each carry one GD. Allele A shown; allele B works the same",
                  notes=(
        "Synthetic tumour SYN-A, four segments (chr1 and chr5 100 Mb, chr2 and chr3 50 Mb). Both children doubled, "
        "so each contributes |CN_c - 2x|; the parent is diploid (1). chr1, chr2: 1 explains both children exactly "
        "(cost 0). chr3: clone5 has a private gain after its GD (3), absorbed at cost 1. chr5: a gain before the "
        "GDs (both children 4) is kept as 2 on clone1. ALPACA had 2, 2, 2, 4."))
    rows = ex.cost_table(A, "clone1", 0)
    tree(s, A, MX, 2.15, 2.75, 2.4, pill=(1.05, 0.4), cn_w=0.9,
         cn={c: " ".join(str(v[0]) for v in A.cn[c]) for c in A.parent if c != "diploid"})
    text(s, MX, 5.35, 3.5, 0.6, "Labels: allele A on chr1, 2, 3, 5 (ALPACA)", size=12, color=MUTED, margin=0,
         align=PP_ALIGN.CENTER)
    head = ["Segment", "clone5 (GD)", "clone6 (GD)", "parent"] + [f"x={i}" for i in range(5)] + ["ALPACA → new"]
    body = [[r["segment"], r["children"]["clone5"], r["children"]["clone6"], r["parent"], *r["costs"],
             f"{r['alpaca']} → {r['new']}"] for r in rows]
    hl = {}
    for ri, r in enumerate(rows, start=1):
        hl[(ri, 4 + r["new"])] = (CORAL, WHITE)
        if r["new"] != r["alpaca"]:
            hl[(ri, 9)] = (PALE_CORAL, CORAL)
    cols = [1.0, 1.1, 1.1, 0.8, 0.6, 0.6, 0.6, 0.6, 0.6, 1.3]
    tx = W - MX - sum(cols)
    table(s, tx, 2.05, sum(cols), [head] + body, cols, size=14, row_h=0.56, highlight=hl, bold_cols=(9,))
    text(s, tx, 5.1, sum(cols), 1.3, [
        [("chr3: ", {"bold": True}), "clone5's private gain after its GD (3) costs 1 at x = 1; still the best."],
        [("chr5: ", {"bold": True}), "a gain before both GDs (children at 4) stays on clone1 as 2."],
    ], size=15, space_after=8, margin=0)


def s_result_a(prs, A):
    s = new_slide(prs, "Example A: the ratios now agree with the calls",
                  "Only clone1 changes; the observed clones keep ALPACA's CN, their ratios follow", notes=(
        "Corrected clone1: 1|1, 1|1, 1|1, 2|1 (ALPACA 2|2, 2|2, 2|2, 4|2). Ratios are max(ratio_A, ratio_B) vs "
        "the parent, recomputed with the ALPACA ratio port. 1.5 is the threshold the WGD pipeline and the "
        "gd_without_alpaca_doubling flag use for 'this edge shows a doubling'. Same numbers as the test "
        "test_parallel_gds_absent_parent_halved (clone1 ratio_A 2.67 -> 1.33, clone5 1.08 -> 2.17)."))
    after = A.cn_after("clone1")
    head = ["Segment", "ALPACA", "corrected"]
    body = [[ex.SEG_LABEL[sg], ab(b4), ab(af)] for sg, b4, af in zip(ex.SEGS, A.cn["clone1"], after)]
    hl = {(i + 1, 2): (PALE_CORAL, CORAL) for i in range(len(body))}
    text(s, MX, 2.0, 4.2, 0.4, "clone1 copy number (A|B)", size=16, bold=True, margin=0)
    table(s, MX, 2.5, 4.2, [head] + body, [1.4, 1.4, 1.4], size=16, row_h=0.52, highlight=hl, bold_cols=(2,))
    text(s, MX, 5.3, 4.2, 1.2, "Synthetic truth was 1|1 with a pre-GD gain on chr5 (2|1): recovered exactly.",
         size=14, color=MUTED, margin=0)
    cd = CategoryChartData()
    clones = list(A.parent)
    cd.categories = clones
    cd.add_series("ALPACA", [round(A.ratio(c, False), 2) for c in clones])
    cd.add_series("corrected", [round(A.ratio(c, True), 2) for c in clones])
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(MX + 4.8), Inches(1.9), Inches(CW - 4.8),
                            Inches(4.55), cd)
    ch = gf.chart
    ch.font.name = BODY
    ch.font.size = Pt(13)
    ch.font.color.rgb = rgb(MUTED)
    ch.has_title = True
    ch.chart_title.text_frame.text = "WGD ratio vs parent, max(A, B)"
    tp = ch.chart_title.text_frame.paragraphs[0]
    tp.runs[0].font.size = Pt(16)
    tp.runs[0].font.bold = True
    tp.runs[0].font.color.rgb = rgb(INK)
    ch.has_legend = True
    ch.legend.position = XL_LEGEND_POSITION.TOP
    ch.legend.include_in_layout = False
    ch.legend.font.size = Pt(13)
    plot = ch.plots[0]
    plot.gap_width = 70
    plot.overlap = -10
    plot.has_data_labels = True
    dl = plot.data_labels
    dl.number_format = "0.00"
    dl.number_format_is_linked = False
    dl.position = XL_LABEL_POSITION.OUTSIDE_END
    dl.font.size = Pt(13)
    dl.font.color.rgb = rgb(INK)
    for ser, col in zip(plot.series, (GRID, CORAL)):
        ser.format.fill.solid()
        ser.format.fill.fore_color.rgb = rgb(col)
    va = ch.value_axis
    va.minimum_scale, va.maximum_scale, va.major_unit = 0, 3, 0.5
    va.has_major_gridlines = True
    va.major_gridlines.format.line.color.rgb = rgb("E3E9EC")
    va.format.line.fill.background()
    va.tick_labels.font.size = Pt(12)
    ca = ch.category_axis
    ca.tick_labels.font.size = Pt(13)
    ca.format.line.color.rgb = rgb(GRID)
    text(s, MX + 4.8, 6.5, CW - 4.8, 0.35, "A GD edge is expected near 2; the pipeline's doubling threshold is 1.5.",
         size=12, color=MUTED, margin=0, align=PP_ALIGN.CENTER)


def s_example_b(prs, B):
    s = new_slide(prs, "Example B: an undoubled child holds the value",
                  "clone1 qualifies, but clone7 (no GD, observed at 2|2) vouches for ALPACA's 2|2", notes=(
        "Synthetic SYN-B. clone1 absent; clone5 doubled (4|4), clone7 not doubled and observed (2|2). The rule "
        "fires (clone1 is listed in corrections.tsv with 0 allele-segments changed) but the optimisation keeps "
        "ALPACA's value: x = 2 explains both children exactly and only pays 1 for the step from diploid. This is "
        "why the rule considers every child, not only the doubled ones."))
    r = ex.cost_table(B, "clone1", 0)[0]
    tree(s, B, MX + 0.1, 2.2, 3.6, 2.4, pill=(1.05, 0.4), cn_w=0.8, cn={c: ab(B.cn[c][0]) for c in B.parent})
    text(s, MX, 5.4, 4.3, 0.6, "Copy number on every segment (ALPACA)", size=12, color=MUTED, margin=0,
         align=PP_ALIGN.CENTER)
    rx = MX + 5.0
    text(s, rx, 2.05, CW - 5.0, 0.4, "Allele A, any segment", size=16, bold=True, margin=0)
    text(s, rx, 2.55, 1.2, 0.5, "cost(x) =", size=18, font=HEAD, margin=0, anchor=MSO_ANCHOR.MIDDLE)
    for i, (term, who) in enumerate((("|4 − 2x|", "clone5 (GD)"), ("|2 − x|", "clone7 (no GD)"),
                                     ("|x − 1|", "diploid parent"))):
        bx = rx + 1.3 + i * 2.05
        tb = box(s, bx, 2.55, 1.6, 0.5, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.2)
        fill_text(tb, term, size=18, font=HEAD, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, margin=0)
        text(s, bx - 0.1, 3.08, 1.8, 0.3, who, size=12, color=MUTED, align=PP_ALIGN.CENTER, margin=0)
        if i < 2:
            text(s, bx + 1.6, 2.55, 0.45, 0.5, "+", size=18, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE,
                 margin=0)
    cost_strip(s, rx, 4.0, r["costs"], r["new"], r["alpaca"])
    b = box(s, rx, 5.25, CW - 5.0, 1.1, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
    fill_text(b, [[("Result: ", {"bold": True}), "x = 2, same as ALPACA. Nothing changes; clone1 is still logged "
                                                 "in corrections.tsv with 0 segments changed."]], size=15,
              anchor=MSO_ANCHOR.MIDDLE, margin=0.15)


def s_chain(prs, C, C2):
    s = new_slide(prs, "Example C: chains are corrected deepest first",
                  "clone1 and clone2 both absent; clone3 and clone4 doubled in parallel; everything 2|2 in ALPACA",
                  notes=(
        "Synthetic SYN-C. Step 1, clone2: children clone3 and clone4 doubled (2 each), parent clone1 still at "
        "ALPACA's 2: x = 1 costs 1, x = 2 costs 4, so clone2 -> 1. Step 2, clone1: its child clone2 is now 1 and "
        "its parent is diploid: x = 1 costs 0. If clone1 were processed first, with clone2 still at 2, x = 1 and "
        "x = 2 would tie and the tie rule would keep ALPACA's 2. In SYN-C2, clone1 has cells: only clone2 is "
        "corrected (an observed clone stops the chain)."))
    after = {c: C.cn_after(c) for c in C.parent}
    tree(s, C, MX + 0.2, 1.95, 3.3, 3.9, pill=(1.05, 0.4), badges={"clone2": 1, "clone1": 2}, cn_w=1.3,
         cn={c: (f"2|2 → {ab(after[c][0])}" if c in C.corrected_clones() else "2|2") for c in C.parent})
    rx = MX + 4.6
    s1 = ex.cost_table(C, "clone2", 0)[0]
    s2 = ex.cost_table(C, "clone1", 0, child_cn={"clone2": after["clone2"]})[0]
    s0 = ex.cost_table(C, "clone1", 0)[0]
    cost_strip(s, rx, 2.45, s1["costs"], s1["new"], s1["alpaca"],
               label=[[("1  clone2", {"bold": True})], [("parent at ALPACA's 2", {"size": 12, "color": MUTED})]])
    cost_strip(s, rx, 3.95, s2["costs"], s2["new"], s2["alpaca"],
               label=[[("2  clone1", {"bold": True})], [("child clone2 now 1", {"size": 12, "color": MUTED})]])
    b = box(s, rx, 5.05, CW - 4.6, 1.45, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
    tie = [i for i, v in enumerate(s0["costs"]) if v == min(s0["costs"])]
    fill_text(b, [
        [("Why deepest first: ", {"bold": True}),
         f"with clone2 still at 2, clone1's costs tie at x = {' and '.join(map(str, tie))} and the tie rule would "
         f"keep ALPACA's 2."],
        [("Observed clone stops the chain: ", {"bold": True}),
         f"if clone1 has cells (C2), only {', '.join(C2.corrected_clones())} is corrected."]],
        size=14, space_after=4, anchor=MSO_ANCHOR.MIDDLE, margin=0.15)


def s_never(prs, D1, E):
    s = new_slide(prs, "What it never changes, and what it only flags",
                  "Disagreements outside the rule are reported in flags.tsv for review, not fixed", notes=(
        "Piotr, 2026-11-13: only missing clones are corrected; ALPACA vs call disagreements on observed clones and "
        "reversions are flagged only, since they should never happen. gd_without_alpaca_doubling: a clone with a "
        "called GD whose max(ratio_A, ratio_B) vs its parent is below 1.5 after correction. alpaca_reversion: a "
        "clone below a GD-carrying clone whose length-weighted median total CN is 2 or less. Examples: SYN-D1 is "
        "example A with clone1 observed; SYN-E is clone1 doubled at 2|2 with clone2 back at 1|1."))
    lw = 4.2
    box(s, MX, 1.95, lw, 4.55, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
    text(s, MX + 0.25, 2.1, lw - 0.5, 0.4, "Never changed", size=18, bold=True, margin=0)
    text(s, MX + 0.25, 2.65, lw - 0.5, 3.7, [
        "Observed clones: cells in any sample.",
        "Clones with a GD called on them.",
        "Absent clones with no GD clone below them (through absent clones).",
        "Tumours without calls or a tree, or below --min-quality (skipped.tsv).",
    ], size=15, bullets=True, space_after=10, margin=0)
    fx = MX + lw + 0.35
    fw = (CW - lw - 0.35 - 0.3) / 2
    f1 = [r for r in D1.res["flags"].iter_rows(named=True)]
    for i, (ttl, desc, case, cn, foot) in enumerate((
            ("gd_without_alpaca_doubling", "GD called, but max(ratio A, B) vs parent < 1.5",
             D1, {c: ab(D1.cn[c][0]) for c in D1.parent},
             "D1: clone1 has cells, so it is kept; " + ", ".join(
                 f"{r['clone']} ({D1.ratio(r['clone'], True):.2f})" for r in f1) + " flagged"),
            ("alpaca_reversion", "below a GD clone, median total CN ≤ 2",
             E, {c: ab(E.cn[c][0]) for c in E.parent},
             "E: clone2 back at 1|1 below doubled clone1: flagged"))):
        x0 = fx + i * (fw + 0.3)
        box(s, x0, 1.95, fw, 4.55, fill=WHITE, line=GRID, lw=1.0, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
        tag = box(s, x0 + 0.2, 2.1, fw - 0.4, 0.42, fill=AMBER, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5)
        fill_text(tag, ttl, size=13, bold=True, color=INK, font=BODY, align=PP_ALIGN.CENTER,
                  anchor=MSO_ANCHOR.MIDDLE, margin=0)
        text(s, x0 + 0.2, 2.62, fw - 0.4, 0.6, desc, size=14, color=INK, margin=0)
        tree(s, case, x0 + 0.3, 3.3, fw - 1.1, 1.75, pill=(0.95, 0.36), font=12, cn=cn, cn_w=0.6)
        text(s, x0 + 0.2, 5.75, fw - 0.4, 0.7, foot, size=12, color=MUTED, margin=0)


def s_outputs(prs):
    s = new_slide(prs, "Outputs and where they go",
                  "One run writes seven files; two of them are tx_data sources", notes=(
        "Command from alpaca_post_hoc_correction/README.md. Every input and the output directory are parameters "
        "(Piotr, 2026-11-13). The two tx_data sources were added in tx_data PR #10 (merged 2026-11-16); the source "
        "directory is a build parameter (--source-dir or TX_ALPACA_CORRECTED_DIR)."))
    files = [("alpaca_corrected.parquet", "tx_data alpaca schema + pred_CN_*_original + corrected", True),
             ("wgd_ratios_corrected.tsv", "ALPACA's 12-column wgd_ratio_scores format", True),
             ("corrections.tsv", "one row per candidate clone: GD children, bp changed, ratios", False),
             ("ratio_comparison.tsv", "every clone, ratios before and after", False),
             ("flags.tsv", "disagreements left after the correction", False),
             ("skipped.tsv", "tumours not corrected, with the reason", False),
             ("run.json", "inputs and counts", False)]
    for i, (f, d, txd) in enumerate(files):
        col, row = i % 2, i // 2
        x0, y0 = MX + col * 6.2, 1.95 + row * 0.98
        b = box(s, x0, y0, 5.9, 0.82, fill=PALE_CORAL if txd else TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE,
                radius=0.12)
        fill_text(b, [[(f, {"bold": True, "font": "Courier New", "size": 14,
                            "color": CORAL if txd else INK})], [(d, {"size": 13, "color": MUTED})]],
                  size=14, anchor=MSO_ANCHOR.MIDDLE, margin=0.15)
    x0, y0 = MX + 6.2, 1.95 + 3 * 0.98
    b = box(s, x0, y0, 5.9, 0.82, fill=WHITE, line=CORAL, lw=1.25, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.12)
    fill_text(b, [[("tx_data sources (PR #10, merged)", {"bold": True})],
                  [("alpaca_corrected · alpaca_wgd_ratios_corrected", {"size": 13, "color": MUTED})]],
              size=14, anchor=MSO_ANCHOR.MIDDLE, margin=0.15)
    cmd = box(s, MX, 5.95, CW, 0.75, fill=INK, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.12)
    fill_text(cmd, "pixi run correct --alpaca … --clone-proportions … --tree … --wgd-calls … --out … "
                   "[--min-quality 0.3]", size=14, font="Courier New", color=WHITE, anchor=MSO_ANCHOR.MIDDLE,
              margin=0.2)


def s_status(prs):
    s = new_slide(prs, "Status and next steps", dark=True, notes=(
        "Mock numbers from alpaca_post_hoc_correction/output/mock/run.json and the README. The tx_data mock "
        "has no absent parent above parallel GDs, so the correction is untested end to end on the mock cohort; "
        "the unit tests cover it on synthetic trees. The four motivating real tumours are listed in the analysis "
        "README."))
    stats = [("9", "mock tumours run"), ("0", "clones corrected on the mock"), ("2", "flags (placeholder calls)"),
             ("34/34", "mock clones: ratio port = ALPACA (≤1e-9)")]
    for i, (n, lab) in enumerate(stats):
        x0 = MX + i * 3.1
        text(s, x0, 1.7, 2.9, 1.0, n, size=54, bold=True, font=HEAD, color=CORAL if i == 1 else WHITE, margin=0)
        text(s, x0, 2.75, 2.8, 0.6, lab, size=14, color="C5D3D8", margin=0)
    text(s, MX, 3.75, 5.8, 0.4, "Next", size=18, bold=True, color=WHITE, margin=0)
    text(s, MX, 4.25, 5.8, 2.4, [
        "Realistic tx_data mocks with absent parents above parallel GDs.",
        "Run on real ALPACA output with the released WGD calls.",
        "Check the four motivating tumours (analysis README).",
        "Micropipeline: WGD steps → correction → tx_data build.",
    ], size=15, color=WHITE, bullets=True, space_after=6, margin=0)
    text(s, MX + 6.5, 3.75, CW - 6.5, 0.4, "Open decision", size=18, bold=True, color=WHITE, margin=0)
    text(s, MX + 6.5, 4.25, CW - 6.5, 2.4, [
        "Extend beyond missing clones once live cases are checked? Today observed clones are never changed and "
        "disagreements are only flagged (Piotr, 2026-11-13)."], size=15, color="C5D3D8", margin=0)


def build(out: Path) -> Path:
    cases = ex.all_cases()
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    prs.core_properties.title = "Post-hoc ALPACA correction"
    prs.core_properties.author = "Darwin (darwins_analysis)"
    s_title(prs, cases["A"])
    s_problem(prs, cases["A"])
    s_rule(prs)
    s_gallery(prs, cases)
    s_how(prs)
    s_example_a(prs, cases["A"])
    s_result_a(prs, cases["A"])
    s_example_b(prs, cases["B"])
    s_chain(prs, cases["C"], cases["C2"])
    s_never(prs, cases["D1"], cases["E"])
    s_outputs(prs)
    s_status(prs)
    prs.save(out)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=HERE / "alpaca_correction.pptx")
    print(f"wrote {build(ap.parse_args().out)}")
