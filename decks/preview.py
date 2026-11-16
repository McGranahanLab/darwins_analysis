"""Approximate PNG preview of a .pptx plus a text-overflow report (no LibreOffice needed).

Draws shapes, connectors, text, tables and bar charts with Pillow. Fonts are DejaVu scaled to
roughly Calibri/Cambria widths (WIDTH_SCALE), so fit is approximate: open the deck in PowerPoint
for the final check.

    python ../preview.py <deck>.pptx preview/
"""

from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE
from pptx.util import Emu

DPI = 110
WIDTH_SCALE = 0.9                    # DejaVu is ~10% wider than Calibri/Cambria
FONT_DIR = Path("/usr/share/fonts/truetype/dejavu")
_cache: dict = {}


def px(emu) -> float:
    return Emu(emu).inches * DPI


def font(size_pt: float, bold: bool, serif: bool, mono: bool):
    name = "DejaVuSansMono" if mono else ("DejaVuSerif" if serif else "DejaVuSans")
    name += "-Bold" if bold else ""
    key = (name, round(size_pt * WIDTH_SCALE, 1))
    if key not in _cache:
        _cache[key] = ImageFont.truetype(str(FONT_DIR / f"{name}.ttf"), max(1, int(size_pt * WIDTH_SCALE * DPI / 72)))
    return _cache[key]


def color(fill_or_font, default=None):
    try:
        return "#" + str(fill_or_font.color.rgb) if hasattr(fill_or_font, "color") else "#" + str(fill_or_font.fore_color.rgb)
    except Exception:
        return default


def shape_fill(s):
    try:
        if s.fill.type == 1:                  # solid
            return "#" + str(s.fill.fore_color.rgb)
    except Exception:
        pass
    return None


def shape_line(s):
    try:
        if s.line.fill.type == 1:
            return "#" + str(s.line.color.rgb)
    except Exception:
        pass
    return None


def run_style(r, size_default=18):
    f = r.font
    size = f.size.pt if f.size else size_default
    name = f.name or "Calibri"
    try:
        col = "#" + str(f.color.rgb)
    except Exception:
        col = "#000000"
    return size, bool(f.bold), name in ("Cambria",), name.startswith("Courier"), col


def layout_text(tf, width_px):
    """Return lines [(segments[(text, font, color)], height)] wrapping runs to width."""
    lines = []
    for p in tf.paragraphs:
        indent = 0.25 * DPI if p._p.pPr is not None and p._p.pPr.get("marL") else 0
        words = []
        for r in p.runs:
            size, bold, serif, mono, col = run_style(r)
            f = font(size, bold, serif, mono)
            for i, w in enumerate(r.text.split(" ")):
                words.append((w + (" " if i < len(r.text.split(" ")) - 1 else ""), f, col, size))
        if not words:
            lines.append(([], 14 * DPI / 72 * 1.2))
            continue
        cur, cur_w, h = [], indent, 0
        for w, f, col, size in words:
            ww = f.getlength(w)
            if cur and cur_w + ww > width_px and w.strip():
                lines.append((cur, h))
                cur, cur_w, h = [], indent, 0
            cur.append((w, f, col, ww))
            cur_w += ww
            h = max(h, size * DPI / 72 * 1.2)
        lines.append((cur, h))
        if p.space_after:
            lines.append(([], p.space_after.pt * DPI / 72))
    return lines


def draw_text(d, tf, x, y, w, h, align_center=False, anchor="t", issues=None, label=""):
    ml, mr = px(tf.margin_left), px(tf.margin_right)
    mt, mb = px(tf.margin_top), px(tf.margin_bottom)
    inner_w = w - ml - mr
    lines = layout_text(tf, inner_w)
    total = sum(hh for _, hh in lines)
    widest = max((sum(s[3] for s in segs) for segs, _ in lines), default=0)
    if issues is not None:
        if total > h - mt - mb + 2:
            issues.append(f"{label}: text height {total / DPI:.2f}in > box {(h - mt - mb) / DPI:.2f}in")
        if widest > inner_w + 2:
            issues.append(f"{label}: line width {widest / DPI:.2f}in > box {inner_w / DPI:.2f}in")
    va = str(tf.vertical_anchor) if tf.vertical_anchor is not None else "TOP"
    cy = y + mt
    if "MIDDLE" in va:
        cy = y + (h - total) / 2
    elif "BOTTOM" in va:
        cy = y + h - mb - total
    for p_lines in lines:
        segs, hh = p_lines
        lw = sum(s[3] for s in segs)
        p0 = tf.paragraphs[0]
        cx = x + ml
        if align_center:
            cx = x + (w - lw) / 2
        for t, f, col, ww in segs:
            d.text((cx, cy), t, font=f, fill=col)
            cx += ww
        cy += hh


def is_center(tf):
    a = tf.paragraphs[0].alignment
    return a is not None and "CENTER" in str(a)


def render(pptx_path: Path, out_dir: Path) -> list[str]:
    prs = Presentation(pptx_path)
    W, H = px(prs.slide_width), px(prs.slide_height)
    out_dir.mkdir(parents=True, exist_ok=True)
    report = []
    for n, slide in enumerate(prs.slides, start=1):
        bg = "#FFFFFF"
        try:
            if slide.background.fill.type == 1:
                bg = "#" + str(slide.background.fill.fore_color.rgb)
        except Exception:
            pass
        img = Image.new("RGB", (int(W), int(H)), bg)
        d = ImageDraw.Draw(img)
        issues = []
        for s in slide.shapes:
            x, y, w, h = px(s.left), px(s.top), px(s.width), px(s.height)
            label = f"slide {n} '{s.name}'"
            if x < -1 or y < -1 or x + w > W + 1 or y + h > H + 1:
                issues.append(f"{label}: outside the slide")
            if s.shape_type == MSO_SHAPE_TYPE.LINE or s.__class__.__name__ == "Connector":
                d.line([(px(s.begin_x), px(s.begin_y)), (px(s.end_x), px(s.end_y))], fill=shape_line(s) or "#888",
                       width=2)
                continue
            if s.has_chart if hasattr(s, "has_chart") else False:
                ch = s.chart
                d.rectangle([x, y, x + w, y + h], outline="#DDDDDD")
                plot = ch.plots[0]
                cats = list(plot.categories)
                sers = list(plot.series)
                vmax = ch.value_axis.maximum_scale or max(max(sr.values) for sr in sers)
                bw = (w - 80) / (len(cats) * (len(sers) + 1))
                for ci, _ in enumerate(cats):
                    for si, sr in enumerate(sers):
                        v = sr.values[ci]
                        bx = x + 50 + (ci * (len(sers) + 1) + si) * bw
                        by = y + h - 40 - (h - 110) * v / vmax
                        try:
                            col = "#" + str(sr.format.fill.fore_color.rgb)
                        except Exception:
                            col = "#888888"
                        d.rectangle([bx, by, bx + bw * 0.9, y + h - 40], fill=col)
                        d.text((bx, by - 18), f"{v:.2f}", font=font(12, False, False, False), fill="#1B2A33")
                    d.text((x + 50 + ci * (len(sers) + 1) * bw, y + h - 34), str(cats[ci]),
                           font=font(13, False, False, False), fill="#5B6B73")
                if ch.has_title:
                    d.text((x + w / 2 - 120, y + 8), ch.chart_title.text_frame.text, font=font(16, True, False, False),
                           fill="#1B2A33")
                continue
            if s.has_table if hasattr(s, "has_table") else False:
                t = s.table
                cy = y
                for ri, row in enumerate(t.rows):
                    cx = x
                    rh = px(row.height)
                    for ci, col in enumerate(t.columns):
                        cw = px(col.width)
                        cell = t.cell(ri, ci)
                        try:
                            fc = "#" + str(cell.fill.fore_color.rgb)
                        except Exception:
                            fc = None
                        d.rectangle([cx, cy, cx + cw, cy + rh], fill=fc, outline="#FFFFFF")
                        draw_text(d, cell.text_frame, cx, cy, cw, rh, is_center(cell.text_frame), issues=issues,
                                  label=f"{label} cell {ri},{ci}")
                        cx += cw
                    cy += rh
                continue
            fill, ln = (shape_fill(s), shape_line(s)) if s.shape_type != MSO_SHAPE_TYPE.TEXT_BOX else (None, None)
            kind = ""
            try:
                kind = str(s.auto_shape_type)
            except Exception:
                pass
            if fill or ln:
                if "OVAL" in kind:
                    d.ellipse([x, y, x + w, y + h], fill=fill, outline=ln, width=2)
                elif "ROUNDED" in kind:
                    r = min(w, h) * (s.adjustments[0] if len(s.adjustments) else 0.16)
                    d.rounded_rectangle([x, y, x + w, y + h], radius=r, fill=fill, outline=ln, width=2)
                else:
                    d.rectangle([x, y, x + w, y + h], fill=fill, outline=ln, width=2)
            if s.has_text_frame and s.text_frame.text.strip():
                draw_text(d, s.text_frame, x, y, w, h, is_center(s.text_frame), issues=issues, label=label)
        p = out_dir / f"slide-{n:02d}.png"
        img.save(p)
        report += issues
        print(p)
    return report


if __name__ == "__main__":
    rep = render(Path(sys.argv[1]), Path(sys.argv[2] if len(sys.argv) > 2 else "preview"))
    print("\n".join(rep) if rep else "no overflow found")
