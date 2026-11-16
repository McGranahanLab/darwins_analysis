"""Build wgd_pipeline.pptx: the rules and heuristics of the WGD pipeline (darwins_analysis/WGD).

Worked examples come from `examples.py`, which runs synthetic tumours through the real step 04
and 05 code; the step 00 example runs `pgdd_input.compute_num_gds`; quality values are read
from `wgd.params` defaults. Native, editable shapes and tables (python-pptx).

    pixi run deck        # writes wgd_pipeline.pptx next to this file
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import polars as pl
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

import examples as ex

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from deckkit import (AMBER, CORAL, CW, GRID, HEAD, INK, MUTED, MX, PALE_CORAL, TEAL, TINT, WHITE, W,  # noqa: E402
                     box, fill_text, line, table, text)
from deckkit import new_slide as _new_slide  # noqa: E402
from wgd import params as P  # noqa: E402
from wgd import pgdd_input  # noqa: E402

FOOTER = "Synthetic examples · darwins_analysis/decks/wgd_pipeline"
PALE_TEAL = "DDEDEE"
PALE_AMBER = "FCEFD9"


def new_slide(prs, title, kicker=None, dark=False, notes=""):
    return _new_slide(prs, title, kicker, dark, notes, footer=FOOTER)


def card(slide, x, y, w, h, paras, fill=TINT, size=15, anchor=MSO_ANCHOR.TOP, margin=0.15, **kw):
    b = box(slide, x, y, w, h, fill=fill, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.07)
    fill_text(b, paras, size=size, anchor=anchor, margin=margin, **kw)
    return b


def bullets(slide, x, y, w, h, items, size=15, **kw):
    return text(slide, x, y, w, h, items, size=size, bullets=True, space_after=7, margin=0, **kw)


def badge(slide, x, y, label, fill=TEAL, color=WHITE, d=0.4, size=13):
    b = box(slide, x, y, d, d, fill=fill, shape=MSO_SHAPE.OVAL)
    fill_text(b, label, size=size, bold=True, color=color, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, margin=0)
    return b


# ------------------------------------------------------------------ clone tree
def tree(slide, case: ex.Case, x, y, w, h, *, gds: dict[str, int] | None = None, labels=None, pill=(1.0, 0.38),
         font=12, label_font=11, highlight=None, label_color=MUTED):
    """Clones as pills (teal = has cells in a sample, white dashed = absent), GDs as amber badges
    (`GD` or `2×`), `labels` under the leaves / beside internal clones, `highlight` clones in coral."""
    gds = case.gds() if gds is None else gds
    highlight = highlight or set()
    present = {c for v in case.cp.values() for c, p in v.items() if p > 0}
    kids: dict[str, list[str]] = {}
    for c, p in case.parent.items():
        kids.setdefault(p, []).append(c)
    for v in kids.values():
        v.sort(key=lambda c: int(c.replace("clone", "")))
    order, depth = [], {}

    def walk(c, d):
        depth[c] = d
        order.append(c)
        for k in kids.get(c, []):
            walk(k, d + 1)

    for r in kids.get("diploid", []):
        walk(r, 0)
    leaves = [c for c in order if not kids.get(c)]
    xs = {c: float(i) for i, c in enumerate(leaves)}

    def place(c):
        if c not in xs:
            xs[c] = sum(place(k) for k in kids[c]) / len(kids[c])
        return xs[c]

    for c in order:
        place(c)
    pw, ph = pill
    levels = max(depth.values()) + 1
    gy = (h - ph) / max(levels - 1, 1)
    span = w - pw
    pos = {c: (x + pw / 2 + (span * xs[c] / max(len(leaves) - 1, 1) if len(leaves) > 1 else span / 2),
               y + depth[c] * gy) for c in order}
    for c in order:
        for k in kids.get(c, []):
            line(slide, pos[c][0], pos[c][1] + ph, pos[k][0], pos[k][1], color=GRID, lw=1.75)
    for c in order:
        cx, cy = pos[c]
        if c in highlight:
            st = dict(fill=PALE_CORAL, line=CORAL, lw=2.0)
            tc = CORAL
        elif c in present:
            st, tc = dict(fill=TEAL), WHITE
        else:
            st, tc = dict(fill=WHITE, line=TEAL, lw=1.75, dash=True), TEAL
        b = box(slide, cx - pw / 2, cy, pw, ph, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5, name=f"node {c}", **st)
        fill_text(b, c, size=font, bold=True, color=tc, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, margin=0)
        n = gds.get(c, 0)
        if n:
            d = 0.38
            badge(slide, cx + pw / 2 - d * 0.5, cy - d * 0.5, "GD" if n == 1 else f"{n}×", fill=AMBER, color=INK,
                  d=d, size=10)
        if labels and c in labels:
            if kids.get(c):
                text(slide, cx + pw / 2 + (0.3 if n else 0.1), cy, 1.6, ph, labels[c], size=label_font, color=label_color,
                     margin=0, anchor=MSO_ANCHOR.MIDDLE)
            else:
                text(slide, cx - pw / 2 - 0.4, cy + ph + 0.03, pw + 0.8, 0.3, labels[c], size=label_font, color=label_color,
                     align=PP_ALIGN.CENTER, margin=0)
    return pos


def legend(slide, x, y):
    cx = x
    for lab, st in (("clone with cells", dict(fill=TEAL)), ("absent from every sample",
                                                            dict(fill=WHITE, line=TEAL, lw=1.75, dash=True))):
        box(slide, cx, y + 0.04, 0.42, 0.22, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5, **st)
        text(slide, cx + 0.5, y, 2.4, 0.3, lab, size=12, color=MUTED, margin=0)
        cx += 2.9
    badge(slide, cx, y - 0.02, "GD", fill=AMBER, color=INK, d=0.34, size=9)
    text(slide, cx + 0.42, y, 2.6, 0.3, "GD on this clone (2× = two)", size=12, color=MUTED, margin=0)


def leaf_labels(case: ex.Case, levels: dict[str, int] | None = None) -> dict[str, str]:
    """Leaf -> 'R1 · target 2' (sample made of that one clone), with the level reached if given."""
    out = {}
    for s, comp in case.cp.items():
        if len(comp) == 1:
            (c,) = comp
            out[c] = f"{s} · {case.num_gds[s]} GD" + ("s" if case.num_gds[s] != 1 else "")
            if levels is not None and levels[s] != case.num_gds[s]:
                out[c] += f" (got {levels[s]})"
    return out


def gds_text(g: dict[str, int]) -> str:
    return " + ".join(f"{c}" + (f" ×{n}" if n > 1 else "") for c, n in sorted(g.items())) or "none"


# ------------------------------------------------------------------ slides
def s_title(prs, nf):
    s = new_slide(prs, "WGD pipeline: rules and heuristics", dark=True, notes=(
        "darwins_analysis/WGD, steps 00 to 07. The deck walks through every rule that turns PGDD's genome-doubling "
        "calls, ALPACA's clone trees and the QC tables into per-clone GD calls with a quality score. Worked examples are "
        "synthetic tumours run through the pipeline's own code (wgd/placement.py, wgd/final_table.py); real tumours that "
        "shaped the rules are in WGD/docs/previous_analysis_logic.md and the regression cases."))
    t = s.shapes.title
    t.top, t.width, t.height = Inches(1.7), Inches(7.6), Inches(1.8)
    t.text_frame.paragraphs[0].runs[0].font.size = Pt(44)
    text(s, MX, 3.65, 7.2, 1.2, "How PGDD calls, ALPACA trees and QC become per-clone genome-doubling calls with "
         "a quality score", size=20, color="C5D3D8", margin=0)
    text(s, MX, 5.2, 7.2, 0.9, ["darwins_analysis/WGD · steps 00–07", "2026-11-16"], size=14, color="9FB0B7", margin=0)
    box(s, 8.4, 1.6, 4.3, 4.5, fill="243742", shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06)
    tree(s, nf, 8.8, 2.05, 3.5, 2.9, labels=leaf_labels(nf), label_color="C5D3D8")
    text(s, 8.6, 5.45, 3.9, 0.45, "A first GD on clone2, a second on clone4", size=12, color="C5D3D8",
         align=PP_ALIGN.CENTER, margin=0)


def s_flow(prs):
    s = new_slide(prs, "From inputs to per-clone calls", "Eight steps; the placement rules live in step 04", notes=(
        "Steps 00 and 03 are Python ports of Piotr's R scripts and match them on synthetic inputs. PGDD runs between 00 "
        "and 01. Steps 04 and 05 hold almost every rule in this deck. Step 06 compares with the previous release table; "
        "step 07 writes the cohort report."))
    steps = [("00", "Sample ploidy", "num_gds per sample from CN; PGDD input"),
             ("01", "Parse PGDD", "events, regions, clusters, mutation counts"),
             ("02", "Prune to tree", "drop clusters ALPACA has no clone for"),
             ("03", "Summary", "per-sample GDs and doubled-genome fractions"),
             ("04", "Place GDs", "QC, triage, fixed GDs, ploidy fill, prior, quality"),
             ("05", "Final table", "GDs per clone, status, class, per-event table"),
             ("06", "Compare", "against the previous release table"),
             ("07", "Report", "patterns, quality, tumours to review")]
    bw, bh, gx = 2.75, 1.45, 0.38
    for i, (n, ttl, d) in enumerate(steps):
        r, c = divmod(i, 4)
        x0, y0 = MX + c * (bw + gx), 2.3 + r * (bh + 0.6)
        hot = n in ("04", "05")
        card(s, x0, y0, bw, bh, [[(ttl, {"bold": True})], [(d, {"size": 13, "color": MUTED if not hot else INK})]],
             fill=PALE_AMBER if hot else TINT, size=15, margin=0.15)
        s.shapes[-1].text_frame.margin_top = Inches(0.33)
        badge(s, x0 + 0.12, y0 - 0.2, n, fill=CORAL if hot else TEAL, d=0.44, size=12)
        if c < 3:
            line(s, x0 + bw + 0.04, y0 + bh / 2, x0 + bw + gx - 0.04, y0 + bh / 2, color=MUTED, lw=1.75, arrow=True)
    text(s, MX, 6.0, CW, 0.8, [[("Inputs: ", {"bold": True}), "CONIPHER tree + MSAI copy number (00) · PGDD output (01) "
                               "· tx_data ALPACA tree, clone proportions and CN (02, 04, 05) · manual QC table + "
                               "clusterInfo (04)"]], size=14, color=MUTED, margin=0)


def s_step00(prs):
    rows = []
    pos = 0
    segs = {"S1": [(60, (2, 2)), (40, (1, 1))], "S2": [(55, (3, 2)), (45, (2, 1))],
            "S3": [(45, (2, 1)), (55, (1, 1))], "S4": [(52, (2.4, 1)), (48, (1, 1))]}
    for smp, parts in segs.items():
        pos = 0
        for L, (a, b) in parts:
            rows.append({"SAMPLE": smp, "STARTPOS": pos + 1, "ENDPOS": pos + L, "COPY_NUMBER_A": float(a),
                         "COPY_NUMBER_B": float(b)})
            pos += L
    res = pgdd_input.compute_num_gds(pl.DataFrame(rows))
    s = new_slide(prs, "Step 00: how many GDs has a sample had?",
                  "Rounded major copy number over the genome, against three thresholds", notes=(
        "wgd/pgdd_input.py compute_num_gds, port of 00_prepare_pgdd_input.R. MajCN = max(A, B) per segment, rounded half "
        "to even like R's round(), so 2.5 rounds to 2 and the effective cut-off for >= 2 is a raw major CN of 1.5. "
        "Fractions are by base pairs. The cut-offs sit below the theoretical 2, 4, 8 on purpose (Piotr, 2026-11-16): "
        "many losses follow a WGD, so observed major CN falls well short of the doubled values. "
        "The table is computed by that function on four synthetic samples. Single-clone "
        "tumours never go to PGDD: all their GDs are clonal, counted as the mode of num_gds; if samples disagree, "
        "quality is Low (single_clone_num_gds_discordant). Any input problem stops the run (pre-flight)."))
    lad = [("3 GDs", "≥ 50% of the genome at major CN ≥ 4"), ("2 GDs", "≥ 50% at major CN ≥ 3"),
           ("1 GD", "≥ 50% at major CN ≥ 2"), ("0 GDs", "otherwise")]
    text(s, MX, 1.95, 5.2, 0.4, "Rule, first match from the top", size=15, bold=True, margin=0)
    for i, (k, d) in enumerate(lad):
        y = 2.45 + i * 0.68
        b = box(s, MX, y, 1.3, 0.52, fill=INK if i < 3 else GRID, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.3)
        fill_text(b, k, size=15, bold=True, color=WHITE if i < 3 else INK, align=PP_ALIGN.CENTER,
                  anchor=MSO_ANCHOR.MIDDLE, margin=0)
        text(s, MX + 1.45, y, 3.9, 0.52, d, size=15, anchor=MSO_ANCHOR.MIDDLE, margin=0)
    text(s, MX, 5.3, 5.3, 1.3, [
        [("Major CN = max(A, B), ", {"bold": True}), "rounded half to even: 2.5 → 2, 1.5 → 2."],
        [("Single-clone tumours ", {"bold": True}), "skip PGDD: all GDs clonal, the mode of num_gds."]],
        size=14, space_after=6, margin=0)
    head = ["Sample", "segments (bp share: A|B)", "≥ 2", "≥ 3", "≥ 4", "num_gds"]
    desc = {k: ", ".join(f"{L}%: {a:g}|{b:g}" for L, (a, b) in v) for k, v in segs.items()}
    body = [[r["SAMPLE"], desc[r["SAMPLE"]], f"{r['frac_maj_ge2']:.2f}", f"{r['frac_maj_ge3']:.2f}",
             f"{r['frac_maj_ge4']:.2f}", str(r["num_gds"])] for r in res.iter_rows(named=True)]
    cols = [0.95, 2.65, 0.7, 0.7, 0.7, 1.05]
    tx = W - MX - sum(cols)
    table(s, tx, 2.45, sum(cols), [head] + body, cols, size=14, row_h=0.55, bold_cols=(5,),
          highlight={(i + 1, 5): (PALE_AMBER, INK) for i in range(len(body))})
    text(s, tx, 5.3, sum(cols), 1.2, [
        [("S3 ", {"bold": True}), "is 5 points short of a GD; ", ("S4 ", {"bold": True}),
         "gets one because 2.4 rounds to 2. How close a sample is to 0.5 drives the quality score and the ploidy prior "
         "later."]], size=14, margin=0)


def s_pgdd(prs):
    s = new_slide(prs, "Steps 01–02: PGDD events and pruning",
                  "Three kinds of event; only clusters ALPACA has a clone for are kept", notes=(
        "PGDD (amf71/ParallelGDDetect, cb0437b) calls GD events per tumour. Only bootstrap-stable calls are written to "
        "the pgdd_* files (n_boot 500, stability 0.7, binomial test, BH-adjusted alpha 0.05); unstable_calls_* is for "
        "manual review only and never read. PGDD runs on pre-tree clusters, so step 02 drops supporting clusters that "
        "are not ALPACA clones. Cluster ids are strings: cluster 3 is clone3."))
    kinds = [("Clonal", "every sample doubled; the GD predates the MRCA", TEAL),
             ("Mutation-supported", "doubled mutations (2 copies) in a cluster place the GD on that cluster", CORAL),
             ("Ploidy-only", "samples are doubled but no cluster carries doubled mutations", AMBER)]
    cw_ = (CW - 2 * 0.3) / 3
    for i, (k, d, col) in enumerate(kinds):
        x0 = MX + i * (cw_ + 0.3)
        card(s, x0, 1.95, cw_, 1.7, [[(k, {"bold": True, "size": 18})], [(d, {"size": 14, "color": INK})]],
             fill=TINT, margin=0.2)
        box(s, x0 + cw_ - 0.5, 2.1, 0.3, 0.3, fill=col, shape=MSO_SHAPE.OVAL)
    text(s, MX, 3.95, 5.8, 0.4, "Step 01: read PGDD", size=16, bold=True, margin=0)
    bullets(s, MX, 4.4, 5.8, 2.2, [
        "Only bootstrap-stable calls (500 resamples, stability ≥ 0.7, BH α 0.05).",
        "Comma-joined clusters and samples become long tables.",
        "Run settings checked against the agreed PGDD settings; differences warned."])
    x2 = MX + 6.3
    text(s, x2, 3.95, CW - 6.3, 0.4, "Step 02: prune to the tree", size=16, bold=True, margin=0)
    bullets(s, x2, 4.4, CW - 6.3, 2.2, [
        "The tree is the set of clones ALPACA has (clone proportions).",
        "Off-tree supporting clusters are removed.",
        "A mutation-supported event left with no cluster becomes ploidy_only_pruned.",
        "Tumours without a tree: no_alpaca_tree in the output."])


def s_qc(prs):
    s = new_slide(prs, "QC: two levels, both visible in the output",
                  "Failing tumours are not placed, but stay in wgd_calls.tsv with NA calls", notes=(
        "Upstream QC: from the manual QC table, rows with selected_solution; the tumour passes if all of them are "
        "FINALIZED. CCF QC: CONIPHER clusterInfo.txt, clusters that are truncal and on the tree; the lowest meanCCF must "
        "be above 90 (percent). Both are evaluated for every PGDD tumour and written to 04_placed/tumour_qc.tsv. Since "
        "2026-11-16 (Piotr) failing tumours appear as status qc_failed with every WGD call NA; QC-passed tumours without "
        "an ALPACA tree as no_alpaca_tree. darwins_analysis PR #6."))
    items = [("Upstream QC", "Manual QC table: every selected solution FINALIZED", "upstream_qc"),
             ("CCF QC", "Lowest truncal tree-cluster meanCCF > 90%", "ccf_qc")]
    for i, (k, d, col) in enumerate(items):
        y0 = 2.0 + i * 1.45
        card(s, MX, y0, 5.6, 1.2, [[(k, {"bold": True, "size": 17})], [(d, {"size": 14})],
                                   [(f"column: {col}", {"size": 12, "color": MUTED, "font": "Courier New"})]],
             margin=0.18)
    line(s, MX + 5.7, 3.25, MX + 6.35, 3.25, color=MUTED, lw=1.75, arrow=True)
    head = ["status", "when", "WGD calls"]
    body = [["resolved / unresolved / failed", "QC passed, placed", "yes"],
            ["qc_failed", "either QC level failed", "NA"],
            ["no_alpaca_tree", "QC passed, ALPACA has no tree", "NA"]]
    table(s, MX + 6.45, 2.0, CW - 6.45, [head] + body, [2.6, 2.3, CW - 6.45 - 4.9], size=14, row_h=0.6,
          highlight={(2, 0): (PALE_CORAL, CORAL), (3, 0): (PALE_CORAL, CORAL)})
    bullets(s, MX, 5.05, CW, 1.5, [
        "Both levels are evaluated for every PGDD tumour; qc_reason says what failed.",
        "The comparison with the previous table (06) skips NA rows; the report (07) lists them in section 4."])


def s_triage(prs):
    s = new_slide(prs, "Triage: four tumour classes", "From PGDD's per-tumour summary; the class is kept in the output",
                  notes=(
        "wgd/placement.py triage (notebook cells [11], [14]). First_GD and num_subclonal_gds come from "
        "pgdd_GDs_per_tumour; mutation support from any event's is_subclonal_mutation_supported. A tumour whose first "
        "GD is subclonal but which also has a clonal second GD falls in the rest group."))
    head = ["triage", "rule", "class in wgd_calls"]
    body = [["no_gd", "First_GD = No GD", "no_wgd"],
            ["clonal_gd_only", "First_GD = Clonal and no subclonal GD", "clonal_wgd"],
            ["alpaca_called", "the rest, no mutation-supported event", "ploidy_only"],
            ["mut_supported", "the rest, ≥ 1 mutation-supported event", "mut_supported"]]
    table(s, MX, 2.0, 7.6, [head] + body, [2.1, 3.6, 1.9], size=15, row_h=0.62)
    card(s, MX + 8.0, 2.0, CW - 8.0, 3.1, [
        [("sub_class", {"bold": True, "size": 17})],
        [("For ploidy_only and mut_supported tumours with every event placed, over the subclonal GD clones:",
          {"size": 14})],
        [("single subclonal GD · parallel only · sequential only · parallel + sequential",
          {"size": 14, "bold": True, "color": TEAL})],
        [("A fill GD on the MRCA counts as clonal.", {"size": 13, "color": MUTED})]], margin=0.2, space_after=6)
    text(s, MX, 5.4, CW, 1.0, "No GD and clonal-only tumours are resolved at once; the other two go through the "
                              "placement rules on the next slides.", size=15, color=MUTED, margin=0)


def s_fixed(prs, fx):
    t = fx.tumour
    s = new_slide(prs, "Fixed GDs go on first", "Clonal and mutation-supported events have a clone; everything else "
                                                "is filled around them", notes=(
        "Rules F and E of the previous notebook, kept. Clonal events go on the MRCA (the clone whose parent is diploid). "
        "Mutation-supported events go on their cluster's clone. A two-cluster event whose clusters are a direct parent-"
        "child pair goes on the child (Piotr, 2026-11-09): on the parent, the child's doubled mutations would need a "
        "second GD in the child, so the pair means a wrong cluster split or an unseen second GD; it costs 0.7. If a "
        "fixed event can't be put on an ALPACA clone the tumour stays unresolved (fixed_event_not_on_tree). The example "
        f"is run through the code: GDs {gds_text(fx.gds())}, status {t['status']}, quality {t['quality']:.2f} "
        f"({t['quality_code']})."))
    rules = [("Clonal", "→ the MRCA (parent = diploid)"),
             ("Mutation-supported", "→ its cluster's clone"),
             ("Two clusters, parent–child", "→ the child; penalty 0.7"),
             ("Not on the tree", "→ tumour unresolved")]
    for i, (k, d) in enumerate(rules):
        y0 = 2.0 + i * 0.95
        card(s, MX, y0, 6.0, 0.78, [[(k + " ", {"bold": True}), (d, {})]], size=16, anchor=MSO_ANCHOR.MIDDLE)
    box(s, MX + 6.5, 2.0, CW - 6.5, 4.55, fill=WHITE, line=GRID, lw=1.0, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
    text(s, MX + 6.75, 2.12, CW - 7.0, 0.4, "Example: clonal GD + event on clusters 2,4", size=15, bold=True, margin=0)
    tree(s, fx, MX + 7.2, 2.75, 3.6, 2.2, labels=leaf_labels(fx))
    text(s, MX + 6.75, 5.55, CW - 7.0, 0.9, [
        [("Result: ", {"bold": True}), f"{gds_text(fx.gds())}; quality {t['quality']:.2f} "
                                       f"({t['quality_code']})."]], size=14, margin=0)


def s_why_fill(prs, nn, nf):
    tn = nn.tumour
    s = new_slide(prs, "Why ploidy-only GDs are not placed one by one",
                  "The old rule gave each event its highest-ratio clone, never the same clone twice", notes=(
        "Previous notebook rule C: for each ploidy-only event in file order, candidates are the clones present in its "
        "samples plus their ancestors, minus the MRCA and clones already taken; pick the highest ALPACA ratio; the "
        "tumour is resolved if every ratio > 1.5. Here R3 needs two GDs after the clonal one. The old rule takes clone3, "
        "then (clone3 being taken) clone2, which also doubles R2 past its single GD; with ratios 1.6 and 1.55 the old "
        f"notebook called it resolved. Today's reversion check fails it ({tn['quality_code']}). The fill puts both on "
        "clone3. Both are outputs of the code on the same synthetic tumour (--ploidy-rule notebook vs the default fill). "
        "Replaced 2026-11-09 after two real tumours."))
    for i, (case, ttl, verdict, col) in enumerate((
            (nn, "Old rule: per event, highest ALPACA ratio",
             "R2 over-doubled; the old notebook still said resolved (both ratios > 1.5)", CORAL),
            (nf, "Ploidy fill: match every sample", "every sample at its PGDD level", TEAL))):
        x0 = MX + i * 6.2
        box(s, x0, 1.95, 5.9, 4.65, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
        text(s, x0 + 0.25, 2.05, 5.4, 0.4, ttl, size=16, bold=True, margin=0)
        tree(s, case, x0 + 0.25, 2.6, 3.2, 2.45, pill=(0.9, 0.34), font=11)
        lv = case.levels()
        body = [[smp, str(case.num_gds[smp]), str(lv[smp])] for smp in case.num_gds]
        hl = {(r + 1, 2): (PALE_CORAL, CORAL) for r, smp in enumerate(case.num_gds) if lv[smp] != case.num_gds[smp]}
        table(s, x0 + 3.7, 2.6, 1.95, [["sample", "PGDD", "got"]] + body, [0.75, 0.6, 0.6], size=12, row_h=0.42,
              highlight=hl)
        comp = "; ".join(f"{k}: " + " + ".join(f"{c} {p:.0%}" for c, p in v.items()) for k, v in case.cp.items()
                         if len(v) > 1)
        text(s, x0 + 3.7, 4.35, 1.95, 0.8, comp, size=10, color=MUTED, margin=0)
        text(s, x0 + 0.25, 5.25, 5.4, 1.3, [
            [("GDs: ", {"bold": True}), gds_text(case.gds())],
            [(verdict, {"bold": True, "color": col})]], size=14, margin=0, space_after=4)
    text(s, MX, 6.65, CW, 0.3, "ALPACA max(A, B) ratio vs parent: clone1 1.9, clone3 1.6, clone2 1.55, others 1.0.",
         size=11, color=MUTED, margin=0)


def s_fill(prs):
    s = new_slide(prs, "The ploidy fill", "Add GDs until every sample has exactly the number of GDs PGDD gives it",
                  notes=(
        "wgd/placement.py fill_ploidy. The ALPACA clone proportions (cp_table) say which clones make up each sample. A "
        "sample is at level k when clones with k or more GDs on their root path are more than 50% of it, the "
        "sample-level analogue of PGDD's 50%-of-genome rule. Fixed GDs (clonal, mutation-supported) go in first; "
        "PGDD's ploidy-only events are superseded. A clone can take several GDs. resolved when every sample ends exactly "
        "at its num_gds (fill = complete); otherwise unresolved with missing:<sample> or over:<sample>."))
    steps = [("1", "Fixed GDs", "clonal on the MRCA, mutation-supported on their clone"),
             ("2", "Sample level", "largest k with > 50% of the sample's cells in clones carrying ≥ k GDs on their "
                                   "root path"),
             ("3", "Add one GD", "on the best clone (next slide); never push a sample above its target"),
             ("4", "Repeat", "until every sample is at its PGDD num_gds, or nothing fits")]
    for i, (n, k, d) in enumerate(steps):
        y0 = 2.0 + i * 1.08
        card(s, MX + 0.6, y0, 6.4, 0.9, [[(k + "  ", {"bold": True}), (d, {"size": 14})]], size=16,
             anchor=MSO_ANCHOR.MIDDLE)
        badge(s, MX, y0 + 0.23, n, fill=TEAL, d=0.44)
    x2 = MX + 7.5
    text(s, x2, 2.0, CW - 7.5, 0.4, "Outcome", size=16, bold=True, margin=0)
    head = ["status", "fill"]
    body = [["resolved", "complete: every sample on target"], ["unresolved", "missing:<sample>=n"],
            ["unresolved", "over:<sample>=n (fixed GDs too many)"]]
    table(s, x2, 2.5, CW - 7.5, [head] + body, [1.45, CW - 7.5 - 1.45], size=14, row_h=0.55)
    text(s, x2, 4.85, CW - 7.5, 1.6, "Ploidy-only events from PGDD are not placed themselves: the fill replaces "
                                     "them (superseded_by_fill), and a different GD count from theirs costs 0.9.",
         size=14, color=MUTED, margin=0)


def s_ties(prs):
    s = new_slide(prs, "Which clone gets the next GD", "Criteria in order; the first that separates the candidates wins",
                  notes=(
        "fill_ploidy candidate ranking (doc section 5, Ploidy fill, step 3). Cell fractions are summed over samples. (b) "
        "and (c) came from a real tumour (2026-11-10), (d) from another where an extinct clone and its only child cover "
        "the same cells and ALPACA can't separate them: the GD goes as high as the data allow. fill_tie (0.8) is "
        "charged only when the ALPACA score alone separated equally deep clones."))
    crit = [("a", "Raises the most samples still below target", "without pushing any sample above"),
            ("b", "Fewest cells above their sample's target", "no over-doubling"),
            ("c", "Most cells at their sample's target", "also covers minor clones of doubled samples"),
            ("d", "Most ancestral", "when the data can't tell clones apart, go as high as allowed"),
            ("e", "Highest ALPACA max(A, B)", "fill_tie 0.8 if this decided it"),
            ("f", "First in the tree", "deterministic last resort")]
    cw_ = (CW - 0.3) / 2
    for i, (k, a, b) in enumerate(crit):
        r, c = divmod(i, 2)
        x0, y0 = MX + c * (cw_ + 0.3), 2.0 + r * 1.5
        card(s, x0 + 0.55, y0, cw_ - 0.55, 1.25, [[(a, {"bold": True})], [(b, {"size": 14, "color": MUTED})]],
             size=16, margin=0.18, anchor=MSO_ANCHOR.MIDDLE)
        badge(s, x0, y0 + 0.4, k, fill=INK, d=0.44)


def s_passes(prs, rl):
    t = rl.tumour
    s = new_slide(prs, "When clone proportions and ploidy disagree", "Two extra passes, both flagged and penalised",
                  notes=(
        "Relaxed pass (2026-11-11): if a sample is still short and no clone is a majority of it without over-doubling "
        "another sample, the GD goes on the clone covering the largest share (>= 10%) of the short samples, as long as "
        "no cell anywhere ends above target. Flag minority_coverage(<share>:<sample>); factor 0.2 at 0% coverage to 0.5 "
        "at >= 50% (Piotr: place it for transparency, penalise heavily). Minor-clone pass (2026-11-12): after every "
        "sample reached its level, GDs also go on clones still below their sample's level if no cell ends above it; "
        "minor_clone 0.6 for a sliver to 0.9 at >= 50% of a sample. The example is run through the code."))
    card(s, MX, 1.95, 6.0, 2.15, [
        [("Relaxed pass", {"bold": True, "size": 17})],
        [("A sample is short and no clone is a clean majority of it: the GD goes on the clone covering the largest "
          "share (≥ 10%).", {"size": 14})],
        [("minority_coverage: 0.2 at 0% → 0.5 at ≥ 50%", {"size": 13, "color": CORAL, "bold": True})]],
        margin=0.2, space_after=4)
    card(s, MX, 4.3, 6.0, 2.15, [
        [("Minor-clone pass", {"bold": True, "size": 17})],
        [("Every sample is at its level, but some clones in it are not: double them too, if nothing ends above "
          "target.", {"size": 14})],
        [("minor_clone: 0.6 for a sliver → 0.9 at ≥ 50%", {"size": 13, "color": CORAL, "bold": True})]],
        margin=0.2, space_after=4)
    x2 = MX + 6.45
    box(s, x2, 1.95, CW - 6.45, 4.5, fill=WHITE, line=GRID, lw=1.0, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
    text(s, x2 + 0.25, 2.07, CW - 6.95, 0.4, "Example: relaxed pass", size=15, bold=True, margin=0)
    tree(s, rl, x2 + 0.9, 2.65, 4.0, 1.9, labels={"clone2": "R1 60%, R2 100%", "clone3": "R3 · 0 GDs",
                                                   "clone4": "R1 40%"})
    text(s, x2 + 0.25, 4.95, CW - 6.95, 1.4, [
        "R1 is doubled (1 GD); R2 and R3 are not.",
        "clone2 would double R2 too; clone4 is only 40% of R1.",
        [("→ clone4, ", {"bold": True}), (f"quality {t['quality']:.2f} ({t['quality_code']})", {})]],
        size=13, margin=0, space_after=3)


def s_prior(prs, ps, pa):
    tp, ta = ps.tumour, pa.tumour
    s = new_slide(prs, "Borderline samples: the ploidy prior",
                  "A sample just under or over 0.5 may move one level if that makes the solution simpler", notes=(
        "wgd/placement.py borderline_search (prior_mode = search, Piotr 2026-11-13). Samples less than 0.10 from 0.5 at "
        "some GD level may move one level across it; every combination moving fewer than half the samples is filled. "
        "Best: (1) simplest fill (missing/excess GDs, forced GDs, majority-fill GDs), (2) fewest PGDD events split, (3) "
        "smallest total move, (4) fewer minor-clone GDs and moves. Used only if simpler than the strict answer. "
        "prior_raised costs 0.6 at the band edge to 0.85 next to 0.5; prior_lowered (diploid prior, stronger) 0.75 to "
        f"0.95. Example run through the code: strict {gds_text(ps.gds())}, {tp['status']}, quality {tp['quality']:.3f}; "
        f"with the prior {gds_text(pa.gds())}, {ta['status']}, quality {ta['quality']:.3f} ({ta['quality_code']})."))
    head = ["Sample", "frac ≥ 2", "PGDD", "clones"]
    body = [[k, f"{v[0]:.2f}", str(ps.num_gds[k]), ", ".join(f"{c} {p:.0%}" for c, p in ps.cp[k].items())]
            for k, v in ps.frac.items()]
    table(s, MX, 2.0, 5.0, [head] + body, [0.9, 1.0, 0.8, 2.3], size=14, row_h=0.55,
          highlight={(4, 1): (PALE_AMBER, INK)})
    text(s, MX, 4.9, 5.0, 1.6, [
        [("R4 ", {"bold": True}), "is 0.02 under 0.5, so it's borderline; R1–R3 are clearly doubled."],
        [("Band: ", {"bold": True}), "0.10 around 0.5. Fewer than half the samples may move."]],
        size=14, margin=0, space_after=6)
    for i, (case, ttl, t) in enumerate(((ps, "Strict targets", tp), (pa, "Prior: R4 raised to 1", ta))):
        x0 = MX + 5.35 + i * 3.45
        box(s, x0, 2.0, 3.3, 4.5, fill=TINT, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
        text(s, x0 + 0.18, 2.1, 2.95, 0.4, ttl, size=14, bold=True, margin=0)
        tree(s, case, x0 + 0.2, 2.8, 2.9, 1.3, pill=(0.8, 0.32), font=10)
        col = CORAL if t["status"] != "resolved" else TEAL
        text(s, x0 + 0.18, 4.45, 2.7, 1.95, [
            [(f"{gds_text(case.gds())}", {"bold": True})],
            [(f"{t['status']}, quality {t['quality']:.2f}", {"color": col, "bold": True})],
            [(re.sub(r"SYN0001_SU_T1-", "", t["fill"] or "") if t["status"] != "resolved" else
              (t["quality_code"] or ""), {"size": 11, "color": MUTED})]], size=13, margin=0, space_after=3)


def s_fail(prs, rv):
    t = rv.tumour
    s = new_slide(prs, "Impossible solutions fail", "A clone can't lose a doubling; checks against ALPACA flag the rest",
                  notes=(
        "Reversion detector (Piotr, 2026-11-12; placement.sample_reversions): a clone with k GDs on its root path, more "
        "than 50% of whose cells sit in samples PGDD calls with fewer than k GDs, fails the tumour: status failed, "
        "quality 0, code reversion_to_diploid or reversion_to_lower_ploidy. It uses the samples, not clone CN. Other "
        "checks lower the score: ALPACA CN reversion (a clone below a GD with median total CN <= 2) 0.2; GDs on a "
        "root-leaf path not matching the samples' num_gds (rule H) 0.5; PGDD events that could not go on a clone 0.5. "
        f"Example run through the code: status {t['status']}, quality {t['quality']:.1f}, {t['quality_code']}."))
    x1 = MX
    box(s, x1, 1.95, 5.6, 4.6, fill=WHITE, line=GRID, lw=1.0, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
    text(s, x1 + 0.25, 2.07, 5.1, 0.4, "Example: mutation-supported GD on clone2", size=15, bold=True, margin=0)
    tree(s, rv, x1 + 0.5, 2.7, 4.6, 2.2, labels=leaf_labels(rv), highlight={"clone5"})
    text(s, x1 + 0.25, 5.4, 5.1, 1.1, [
        [("clone5 ", {"bold": True}), "inherits the GD, but all of it is R3, which PGDD calls undoubled."],
        [(f"→ {t['status']}, quality {t['quality']:.0f}", {"bold": True, "color": CORAL})]],
        size=13, margin=0, space_after=3)
    x2 = MX + 6.0
    head = ["check", "effect"]
    body = [["clone with k GDs, > 50% of its cells in samples with < k (reversion detector)", "status failed, quality 0"],
            ["ALPACA CN back to ≤ 2 below a GD (reversion)", "× 0.2"],
            ["GDs on a root–leaf path ≠ samples' num_gds (H)", "× 0.5"],
            ["PGDD events not on an ALPACA clone", "× 0.5"]]
    table(s, x2, 1.95, CW - 6.0, [head] + body, [CW - 6.0 - 1.85, 1.85], size=14, row_h=0.85,
          highlight={(1, 1): (PALE_CORAL, CORAL)})


def s_quality(prs):
    p = P.Params()
    f = p.fixed
    s = new_slide(prs, "Quality: a product of evidence factors",
                  "1 = no concern; every factor below 1 is listed in quality_code with what triggered it", notes=(
        "wgd/quality.py; values from wgd.params defaults (params.toml), read when the deck was built. Graded factors run "
        "linearly between the values shown. alpaca_disagrees uses max(ratio_A, ratio_B); for a mutation-supported GD "
        f"its floor rises to {p.alpaca_mut_floor} at >= {p.min_support} doubled mutations (mutations are trusted over "
        "ALPACA). Parallel GDs under an absent parent that ALPACA shows doubled are scored on the parent's ratio "
        "(alpaca_parent). GD order is deliberately not used (Piotr, 2026-11-10). Weights are defaults accepted as a "
        "starting point, to be calibrated on the cohort."))
    graded = [["borderline_ploidy", f"0.5 → 1 at ≥ {p.ploidy_margin:g} from 0.5", "a sample's genome fraction near 0.5"],
              ["low_mut_support", f"0.5 → 1 at ≥ {p.min_support} mutations", "few doubled mutations"],
              ["alpaca_disagrees", f"0.5 at 1 → 1 at ≥ {p.alpaca_threshold:g}", "GD clone's ALPACA ratio low"],
              ["minority_coverage", "0.2 → 0.5 at ≥ 50%", "forced GD (relaxed pass)"],
              ["minor_clone", f"{p.minor_clone_factor[0]:g} → {p.minor_clone_factor[1]:g}", "minor-clone pass GD"],
              ["prior_raised / _lowered", f"{p.prior_raised_factor[0]:g}–{p.prior_raised_factor[1]:g} / "
                                          f"{p.prior_lowered_factor[0]:g}–{p.prior_lowered_factor[1]:g}",
               "sample moved by the prior"]]
    fixed_ = [[k, f"{v:g}"] for k, v in f.items() if k != "reversion_to_diploid"]
    table(s, MX, 1.95, 7.55, [["graded factor", "range", "when"]] + graded, [2.45, 2.55, 2.55], size=13, row_h=0.6)
    table(s, MX + 7.95, 1.95, CW - 7.95, [["fixed factor", "value"]] + fixed_, [CW - 7.95 - 0.9, 0.9], size=13,
          row_h=0.52)
    text(s, MX, 6.3, CW, 0.5, "A failed solution (reversion detector) sets quality to 0. In the fixed-GD example (slide 7), "
                              "multicluster_parent_child alone gives 0.70.", size=13, color=MUTED, margin=0)


def s_events(prs):
    s = new_slide(prs, "Per-event quality and rejection", "Each GD has its own score; one cut-off can drop weak GDs",
                  notes=(
        "Step 05 writes wgd_events.tsv, one row per GD with id tumour:clone:k, its source (clonal, mutation, fill, fill "
        "forced/minor), PGDD event, quality and factors. The tumour score is the product over factor codes of the worst "
        "value among its events and tumour-level factors (path mismatch, unplaced events, fill-count difference are "
        "tumour-only). --reject-below Q drops GDs with event quality < Q and recomputes counts, sub_class and scores, "
        "writing wgd_calls_rejected_below_<Q>.tsv and a report of each rejected GD. Mutation-supported GDs are kept "
        "unless --reject-mutation-supported (Piotr, 2026-11-11: PGDD's bootstrap already triaged them). Rejected GDs "
        "are not penalised; the samples they explained are then under-explained."))
    rows = [("Per event", "wgd_events.tsv: one row per GD (tumour:clone:k), source, factors, event_quality"),
            ("Per clone", "worst_event_quality on each GD-carrying clone in wgd_calls.tsv"),
            ("Per tumour", "quality: worst value per factor over its events + tumour-level factors, multiplied")]
    for i, (k, d) in enumerate(rows):
        y0 = 2.0 + i * 0.95
        card(s, MX, y0, 6.3, 0.8, [[(k + "  ", {"bold": True}), (d, {"size": 14})]], size=16, anchor=MSO_ANCHOR.MIDDLE)
    card(s, MX + 6.7, 2.0, CW - 6.7, 2.7, [
        [("--reject-below Q", {"bold": True, "font": "Courier New", "size": 16})],
        [("Drops GDs with event quality < Q (one cut-off for the cohort) and recomputes the table.", {"size": 14})],
        [("Mutation-supported GDs are kept unless --reject-mutation-supported.", {"size": 14, "bold": True})],
        [("Every drop is listed: rejected_events_below_<Q>.tsv", {"size": 13, "color": MUTED})]],
        fill=PALE_AMBER, margin=0.2, space_after=6)
    text(s, MX, 5.05, CW, 1.4, [
        [("Why per event: ", {"bold": True}), "one weak, forced GD should not hide three solid ones, and a reviewer "
                                             "can filter on the reason."],
        [("Caveat: ", {"bold": True}), "a rejected GD is not penalised, but the samples it explained are left "
                                       "under-explained."]], size=14, margin=0, space_after=6)


def s_outputs(prs):
    s = new_slide(prs, "What comes out", "Step 05 tables and the step 07 report", notes=(
        "wgd_calls.tsv keeps the previous release format (tumour_id, clone, parent, number_of_gds_at_clone, "
        "total_gds_per_tumour, quality, status, class, sub_class) plus quality_code, worst_event_quality, n_unplaced_"
        "events, n_rejected_events, ploidy_prior and the QC columns. modified_tumours.tsv lists every departure from the "
        "plain rules for manual validation (Piotr, 2026-11-12). Report review severities: critical = failed or "
        "unresolved; serious = forced GD, several GDs on a minor clone, ALPACA reversion, path mismatch, unplaced events, "
        "quality < 0.3; warning = ALPACA CN to correct (parallel GDs), prior moves >= 0.05, ambiguous placement."))
    files = [("wgd_calls.tsv", "one row per ALPACA clone: GDs, status, class, sub_class, quality, QC columns"),
             ("wgd_events.tsv", "one row per GD, with source and quality factors"),
             ("modified_tumours.tsv", "every prior move, forced GD, minor-clone GD and failure, strict vs final"),
             ("07_report/report.html", "patterns, quality, worked examples, tumours to review, tumours without calls")]
    for i, (f, d) in enumerate(files):
        y0 = 1.95 + i * 0.95
        card(s, MX, y0, 6.6, 0.8, [[(f, {"bold": True, "font": "Courier New", "size": 14})],
                                   [(d, {"size": 13, "color": MUTED})]], anchor=MSO_ANCHOR.MIDDLE, margin=0.15)
    head = ["review", "reasons"]
    body = [["critical", "failed (reversion), unresolved"],
            ["serious", "forced GD, stacked minor-clone GDs, ALPACA reversion, path mismatch, unplaced events, "
                        "quality < 0.3"],
            ["warning", "ALPACA CN to correct, prior moved a sample ≥ 0.05 from 0.5, ambiguous placement"]]
    table(s, MX + 7.0, 1.95, CW - 7.0, [head] + body, [1.25, CW - 7.0 - 1.25], size=13, row_h=1.05,
          highlight={(1, 0): (PALE_CORAL, CORAL), (2, 0): (PALE_AMBER, INK)})


def s_open(prs, n_cases):
    s = new_slide(prs, "Tuning, safeguards and open questions", dark=True, notes=(
        "Every threshold and weight is in WGD/params.toml (--params), recorded in run.json. Real tumours discussed with "
        "Piotr are frozen in WGD/tests/cases with their agreed answers (agreed.toml); a rule change that breaks one is a "
        "clash to discuss before merging, not a test to update. Open items from WGD/TODO.md and doc sections 5b, 6, 7."))
    stats = [(str(n_cases), "real tumours frozen as regression cases"), ("1", "parameter file (params.toml)"),
             ("0.5", "the cut-off behind both PGDD's ploidy and the sample level")]
    for i, (n, lab) in enumerate(stats):
        x0 = MX + i * 4.1
        text(s, x0, 1.7, 3.8, 1.0, n, size=54, bold=True, font=HEAD, color=CORAL if i == 0 else WHITE, margin=0)
        text(s, x0, 2.75, 3.8, 0.6, lab, size=14, color="C5D3D8", margin=0)
    text(s, MX, 3.7, 5.8, 0.4, "Safeguards", size=18, bold=True, color=WHITE, margin=0)
    bullets(s, MX, 4.2, 5.8, 2.4, [
        "A rule change that breaks an agreed answer is discussed before merging.",
        "Every departure from the plain rules is listed for manual review.",
        "Failing and tree-less tumours stay visible, with NA calls."], color=WHITE)
    text(s, MX + 6.5, 3.7, CW - 6.5, 0.4, "Open", size=18, bold=True, color=WHITE, margin=0)
    bullets(s, MX + 6.5, 4.2, CW - 6.5, 2.4, [
        "Is 1.5 the right ALPACA ratio cut-off for the new ALPACA?",
        "Is alpaca_disagrees informative for subclonal GDs (ratios 1.0–1.3)?",
        "Calibrate factor weights on the cohort.",
        "Second GDs, triage classes, events on > 2 clusters."], color="C5D3D8")


def build(out: Path) -> Path:
    cases = ex.all_cases()
    n_cases = sum(1 for d in (ex.WGD / "tests" / "cases").iterdir() if d.is_dir())
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(7.5)
    prs.core_properties.title = "WGD pipeline: rules and heuristics"
    prs.core_properties.author = "Darwin (darwins_analysis)"
    s_title(prs, cases["nested_fill"])
    s_flow(prs)
    s_step00(prs)
    s_pgdd(prs)
    s_qc(prs)
    s_triage(prs)
    s_fixed(prs, cases["fixed"])
    s_why_fill(prs, cases["stacked_notebook"], cases["stacked_fill"])
    s_fill(prs)
    s_ties(prs)
    s_passes(prs, cases["relaxed"])
    s_prior(prs, cases["prior_strict"], cases["prior_applied"])
    s_fail(prs, cases["reversion"])
    s_quality(prs)
    s_events(prs)
    s_outputs(prs)
    s_open(prs, n_cases)
    prs.save(out)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=HERE / "wgd_pipeline.pptx")
    print(f"wrote {build(ap.parse_args().out)}")
