"""Shared slide primitives for the decks in this folder (python-pptx, native editable shapes).

Palette, fonts and 16:9 geometry; boxes, text, connectors, tables and a slide frame with a
title, kicker, footer and page number. Each deck adds its own diagrams on top.
"""

from __future__ import annotations

from lxml import etree
from pptx.dml.color import RGBColor
from pptx.enum.dml import MSO_LINE_DASH_STYLE
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

# ------------------------------------------------------------------ style
HEAD = "Cambria"
BODY = "Calibri"
INK = "1B2A33"        # text, dark slides
TEAL = "1F6F78"       # observed clones, primary
CORAL = "D1495B"      # corrected values
AMBER = "E8A33D"      # genome doublings
MUTED = "5B6B73"      # captions
GRID = "C9D4D9"       # edges, rules
TINT = "EEF3F4"       # cards
WHITE = "FFFFFF"
PALE_CORAL = "FBE7EA"

W, H = 13.333, 7.5
MX = 0.6
CW = W - 2 * MX


def rgb(h: str) -> RGBColor:
    return RGBColor.from_string(h)


# ------------------------------------------------------------------ primitives
def no_shadow(shape):
    shape.shadow.inherit = False
    return shape


def box(slide, x, y, w, h, *, fill=None, line=None, lw=1.0, dash=False, shape=MSO_SHAPE.RECTANGLE, radius=None,
        name=None):
    s = no_shadow(slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h)))
    if radius is not None:
        s.adjustments[0] = radius
    if fill:
        s.fill.solid()
        s.fill.fore_color.rgb = rgb(fill)
    else:
        s.fill.background()
    if line:
        s.line.color.rgb = rgb(line)
        s.line.width = Pt(lw)
        if dash:
            s.line.dash_style = MSO_LINE_DASH_STYLE.DASH
    else:
        s.line.fill.background()
    if name:
        s.name = name
    return s


def _runs(p, runs, size, color, bold, italic, font):
    if isinstance(runs, str):
        runs = [runs]
    for r in runs:
        opt = {} if isinstance(r, str) else r[1]
        t = r if isinstance(r, str) else r[0]
        run = p.add_run()
        run.text = t
        f = run.font
        f.size = Pt(opt.get("size", size))
        f.bold = opt.get("bold", bold)
        f.italic = opt.get("italic", italic)
        f.name = opt.get("font", font)
        f.color.rgb = rgb(opt.get("color", color))
        if opt.get("sup"):
            run._r.get_or_add_rPr().set("baseline", "30000")
        if opt.get("sub"):
            run._r.get_or_add_rPr().set("baseline", "-25000")


def fill_text(shape, paras, *, size=16, color=INK, bold=False, italic=False, font=BODY, align=PP_ALIGN.LEFT,
              anchor=MSO_ANCHOR.TOP, margin=0.05, space_after=0, bullets=False):
    tf = shape.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for side in ("left", "right", "top", "bottom"):
        setattr(tf, f"margin_{side}", Inches(margin))
    if isinstance(paras, str) or (isinstance(paras, list) and paras and isinstance(paras[0], tuple)):
        paras = [paras]
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        if space_after:
            p.space_after = Pt(space_after)
        if bullets:
            pPr = p._p.get_or_add_pPr()
            pPr.set("marL", str(Inches(0.25)))
            pPr.set("indent", str(-Inches(0.2)))
            bu = etree.SubElement(pPr, qn("a:buChar"))
            bu.set("char", "•")
        _runs(p, para, size, color, bold, italic, font)
    return shape


def text(slide, x, y, w, h, paras, name=None, **kw):
    s = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        s.name = name
    return fill_text(s, paras, **kw)


def line(slide, x1, y1, x2, y2, *, color=GRID, lw=1.5, arrow=False, dash=False):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1), Inches(x2), Inches(y2))
    c.line.color.rgb = rgb(color)
    c.line.width = Pt(lw)
    if dash:
        c.line.dash_style = MSO_LINE_DASH_STYLE.DASH
    if arrow:
        ln = c.line._get_or_add_ln()
        tail = etree.SubElement(ln, qn("a:tailEnd"))
        tail.set("type", "triangle")
        tail.set("w", "med")
        tail.set("h", "med")
    return c


def set_notes(prs, slide, notes: str) -> None:
    """Speaker notes. A deck re-saved by PowerPoint may have a notes master without a body placeholder,
    so new notes pages get none: copy the body placeholder from another slide's notes page."""
    ns = slide.notes_slide
    if ns.notes_text_frame is None:
        import copy
        for other in prs.slides:
            if other is not slide and other.has_notes_slide and other.notes_slide.notes_placeholder is not None:
                ns.shapes._spTree.append(copy.deepcopy(other.notes_slide.notes_placeholder._element))
                break
    ns.notes_text_frame.text = notes


def new_slide(prs, title, kicker=None, dark=False, notes="", footer=""):
    s = prs.slides.add_slide(prs.slide_layouts[5])          # "Title Only"
    if dark:
        s.background.fill.solid()
        s.background.fill.fore_color.rgb = rgb(INK)
    t = s.shapes.title
    t.left, t.top, t.width, t.height = Inches(MX), Inches(0.42), Inches(CW), Inches(0.8)
    t.text_frame.clear()
    fill_text(t, title, size=34, bold=True, font=HEAD, color=WHITE if dark else INK, anchor=MSO_ANCHOR.BOTTOM,
              margin=0)
    if kicker:
        text(s, MX, 1.24, CW, 0.42, kicker, size=17, color="C5D3D8" if dark else MUTED, margin=0)
    n = len(prs.slides)
    text(s, MX, 7.06, 8.0, 0.3, footer, size=10,
         color="9FB0B7" if dark else MUTED, margin=0)
    text(s, W - MX - 1.0, 7.06, 1.0, 0.3, str(n), size=10, color="9FB0B7" if dark else MUTED, margin=0,
         align=PP_ALIGN.RIGHT)
    if notes:
        set_notes(prs, s, notes)
    return s



def table(slide, x, y, w, rows, col_w, *, header_fill=INK, size=14, row_h=0.42, highlight=None, bold_cols=()):
    """rows[0] = header. highlight: {(r, c): (fill, font_color)}."""
    highlight = highlight or {}
    shp = slide.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(w),
                                 Inches(row_h * len(rows)))
    t = shp.table
    for i, cw in enumerate(col_w):
        t.columns[i].width = Inches(cw)
    for r, row in enumerate(rows):
        t.rows[r].height = Inches(row_h)
        for c, val in enumerate(row):
            cell = t.cell(r, c)
            cell.margin_left = cell.margin_right = Inches(0.06)
            cell.margin_top = cell.margin_bottom = Inches(0.02)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            fill, fc = (header_fill, WHITE) if r == 0 else ((TINT if r % 2 == 0 else WHITE), INK)
            if (r, c) in highlight:
                fill, fc = highlight[(r, c)]
            cell.fill.solid()
            cell.fill.fore_color.rgb = rgb(fill)
            cell.text_frame.clear()
            fill_text(cell, val if isinstance(val, list) else str(val), size=size, color=fc,
                      bold=(r == 0 or c in bold_cols), align=PP_ALIGN.CENTER if c else PP_ALIGN.LEFT,
                      anchor=MSO_ANCHOR.MIDDLE, margin=0.04)
    return shp
