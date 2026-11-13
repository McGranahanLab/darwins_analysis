"""Step 07: cohort summary of the WGD calls (Piotr, 2026-11-13).

Reads steps 04 and 05 and answers three questions, as tables, figures and one HTML report:

1. Quality per WGD pattern, with six worked examples (2 worst, 2 mid, 2 best) among tumours
   with subclonal or parallel GDs, and a sentence for every factor that lowered their score.
2. The cohort: how many tumours have GDs, in which pattern, how many GDs, from what evidence,
   and which quality factors fire how often.
3. Tumours to review by hand: impossible or unresolved solutions, forced GDs, low quality and
   other signs that the copy-number solution (clone proportions, tree, ALPACA CN) is wrong.

WGD pattern (from the GDs per clone; the MRCA is the clone whose parent is `diploid`):
  no GD | clonal only | clonal + subclonal | subclonal, one lineage | subclonal, parallel
"""

from __future__ import annotations

import base64
import html
import io
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

PATTERNS = ("no GD", "clonal only", "clonal + subclonal", "subclonal, one lineage", "subclonal, parallel")
SUBCLONAL_PATTERNS = PATTERNS[2:]
# reference categorical palette, fixed order (validated: dataviz validate_palette.js, light surface)
PATTERN_COLOUR = dict(zip(PATTERNS, ("#898781", "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4")))
SEVERITY_COLOUR = {"critical": "#d03b3b", "serious": "#ec835a", "warning": "#fab219"}
SEVERITY_ORDER = {"critical": 0, "serious": 1, "warning": 2}
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SOURCE_LABEL = {"clonal": "clonal (PGDD)", "mutation": "mutation-supported", "mutation_multicluster": "mutation-supported",
                "fill": "ploidy fill", "fill_forced": "ploidy fill, forced", "fill_minor": "ploidy fill, minor clones"}


# ---- per-tumour table ----------------------------------------------------------------

def _paths(parent: dict[str, str]) -> dict[str, list[str]]:
    out = {}
    for c in parent:
        p, path = c, []
        while p in parent and p not in path:
            path.append(p)
            p = parent[p]
        out[c] = path
    return out


def _factors(tumour_factors: str | None, event_factors: list[str | None]) -> dict[str, tuple[float, str]]:
    """Worst value per quality factor code: {code: (value, detail)}."""
    worst: dict[str, tuple[float, str]] = {}
    for js in [tumour_factors] + list(event_factors):
        for code, f, detail in (json.loads(js) if js else []):
            if code not in worst or f < worst[code][0]:
                worst[code] = (f, detail)
    return worst


def tumour_table(tumours: pl.DataFrame, calls: pl.DataFrame, events: pl.DataFrame,
                 summary: pl.DataFrame | None = None) -> pl.DataFrame:
    """One row per tumour: pattern, GD counts, quality, main penalty, flags."""
    ev_by_t = {t: g for (t,), g in events.group_by("tumour_id")} if events.height else {}
    calls_by_t = {t: g for (t,), g in calls.group_by("tumour_id")}
    n_samples = (dict(summary.group_by("tumour_id").len().iter_rows()) if summary is not None
                 and "tumour_id" in summary.columns else {})
    rows = []
    for r in tumours.iter_rows(named=True):
        t = r["tumour_id"]
        g = calls_by_t.get(t)
        if g is None:
            continue
        parent = dict(zip(g["clone"], g["parent"]))
        gds = {c: n for c, n in zip(g["clone"], g["number_of_gds_at_clone"]) if n}
        mrca = next((c for c, p in parent.items() if p == "diploid"), None)
        paths = _paths(parent)
        n_clonal = gds.get(mrca, 0)
        sub = [c for c in gds if c != mrca]
        n_sub = sum(gds[c] for c in sub)
        lineages = [c for c in sub if not any(a in sub for a in paths[c][1:])]  # top-most subclonal GD clones
        if not gds:
            pattern = "no GD"
        elif not sub:
            pattern = "clonal only"
        elif n_clonal:
            pattern = "clonal + subclonal"
        else:
            pattern = "subclonal, parallel" if len(lineages) > 1 else "subclonal, one lineage"
        te = ev_by_t.get(t)
        flags = ";".join(x for x in (te["flags"].drop_nulls().to_list() if te is not None else []))
        factors = _factors(r.get("tumour_factors"),
                           te["event_factors"].to_list() if te is not None and "event_factors" in te.columns else [])
        main = min(factors.items(), key=lambda kv: kv[1][0]) if factors else None
        rows.append({
            "tumour_id": t, "pattern": pattern, "evidence_class": g["class"][0], "status": r["status"],
            "quality": r["quality"], "quality_code": r["quality_code"],
            "main_penalty": f"{main[0]}({main[1][1]})" if main else None,
            "main_penalty_value": main[1][0] if main else None,
            "n_gds": sum(gds.values()), "clonal_gds": n_clonal, "subclonal_gds": n_sub,
            "parallel_lineages": len(lineages) if not n_clonal else (len(lineages) if len(lineages) > 1 else 0),
            "max_gds_on_path": max((sum(gds.get(a, 0) for a in paths[c]) for c in parent), default=0),
            "gd_clones": ",".join(f"{c}x{n}" if n > 1 else c for c, n in sorted(gds.items())) or None,
            "n_samples": n_samples.get(t),
            "ploidy_prior": r.get("ploidy_prior"),
            "forced_gds": flags.count("minority_coverage("),
            "minor_clone_gds": flags.count("minor_clone_fill("),
            "alpaca_parent": flags.count("alpaca_parent(") > 0,
            "factors_json": json.dumps({k: [round(v[0], 3), v[1]] for k, v in sorted(factors.items())}),
        })
    return pl.DataFrame(rows, infer_schema_length=None).sort("tumour_id")


# ---- 1. examples with commentary ----------------------------------------------------

def explain(code: str, value: float, detail: str) -> str:
    """One sentence per quality factor, for the worked examples and the review list."""
    d = detail or ""
    parts = d.split(":")
    v = f"x{value:.2f}"
    if code == "borderline_ploidy":
        return (f"{v} borderline ploidy: sample {parts[-1]} is only {float(parts[0]):.3f} from the 0.5 cut-off "
                f"at its GD level, so its GD count could go either way.")
    if code == "prior_raised":
        return (f"{v} ploidy prior: sample {parts[0]} ({float(parts[1]):.3f} below 0.5) was raised one GD level "
                f"to fit a simpler placement.")
    if code == "prior_lowered":
        return (f"{v} ploidy prior: sample {parts[0]} ({float(parts[1]):.3f} above 0.5) was lowered one GD level "
                f"to fit a simpler placement.")
    if code == "alpaca_disagrees":
        mut = f", softened by {parts[2][3:]} doubled mutations" if len(parts) > 2 and parts[2].startswith("mut") else ""
        return (f"{v} ALPACA disagrees: {parts[1]} has a copy-number ratio of only {parts[0]} to its parent "
                f"(1.5 expected for a GD){mut}.")
    if code == "low_mut_support":
        return f"{v} weak mutation support: only {d} doubled mutations back the mutation-supported GD."
    if code == "minority_coverage":
        return (f"{v} forced GD: no clone is the majority of sample {parts[-1]} without over-doubling others, "
                f"so a GD went on a clone covering only {float(parts[0]):.0%} of it. Clone proportions and "
                f"sample ploidy disagree.")
    if code == "minor_clone":
        return (f"{v} minor-clone GD: a GD added only to double a clone making up {float(parts[0]):.0%} "
                f"of sample {parts[-1]}.")
    if code == "ambiguous_solution":
        return f"{v} ambiguous: flipping a borderline sample or moving the 50% cut-off gives another placement ({d})."
    if code == "multicluster_parent_child":
        return f"{v} PGDD supported the GD with two clusters (parent and child); it was put on the child."
    if code == "fill_tie":
        return f"{v} tie: equally good clones for a GD were separated only by the ALPACA score ({d})."
    if code == "fill_count_differs":
        return f"{v} the fill added a different number of GDs than PGDD's ploidy-only events ({d})."
    if code == "ploidy_path_mismatch":
        return f"{v} GDs on a root-to-leaf path don't match the samples' maximum GD count ({d})."
    if code == "unplaced_events":
        return f"{v} {d} PGDD event(s) could not be put on an ALPACA clone."
    if code == "reversion":
        return f"{v} ALPACA CN shows clone(s) {d} below a doubled clone looking diploid again."
    if code.startswith("reversion_to_"):
        return (f"{v} impossible: clone(s) below a GD sit mostly in samples with fewer GDs ({d}); "
                f"a clone can't lose a GD, so the placement fails.")
    return f"{v} {code}({d})."


def _commentary(r: dict) -> str:
    t = r["tumour_id"]
    head = (f"{t}: {r['pattern']}, {r['n_gds']} GD(s) on {r['gd_clones'] or 'no clone'} "
            f"({r['clonal_gds']} clonal, {r['subclonal_gds']} subclonal), quality {r['quality']:.3f}.")
    factors = json.loads(r["factors_json"])
    lines = [explain(k, v, d) for k, (v, d) in sorted(factors.items(), key=lambda kv: kv[1][0]) if v < 1]
    if not lines:
        lines = ["Nothing lowered the score: every sample is clear of the 0.5 cut-offs, the placement is "
                 "unique, and ALPACA (or doubled mutations) agree with every GD."]
    if r["ploidy_prior"] and not any(k.startswith("prior_") for k in factors):
        lines.append(f"Ploidy prior considered but not used ({r['ploidy_prior']}).")
    if r["alpaca_parent"]:
        lines.append("Parallel GDs under a parent absent from every sample: ALPACA put their shared gain on the "
                     "parent, so each GD was judged on the parent's ratio (alpaca_parent).")
    return " ".join([head] + lines)


def pick_examples(tt: pl.DataFrame, n: int = 2) -> pl.DataFrame:
    """2 worst, 2 closest to the median, 2 best among tumours with subclonal or parallel GDs."""
    pool = tt.filter(pl.col("pattern").is_in(list(SUBCLONAL_PATTERNS)) & pl.col("quality").is_not_null()).sort("quality")
    if pool.is_empty():
        return pool.with_columns(pl.lit(None, pl.Utf8).alias("band"), pl.lit(None, pl.Utf8).alias("commentary"))
    worst = pool.head(n)["tumour_id"].to_list()
    best = [x for x in pool.tail(n)["tumour_id"].to_list() if x not in worst]
    med = float(pool["quality"].median())
    rest = pool.filter(~pl.col("tumour_id").is_in(worst + best))
    mid = (rest.with_columns((pl.col("quality") - med).abs().alias("_d")).sort("_d", "tumour_id")
           .head(n)["tumour_id"].to_list())
    band = {**{x: "low" for x in worst}, **{x: "mid" for x in mid}, **{x: "high" for x in best}}
    ex = pool.filter(pl.col("tumour_id").is_in(list(band))).with_columns(
        pl.col("tumour_id").replace_strict(band).alias("band"))
    ex = ex.with_columns(pl.Series("commentary", [_commentary(r) for r in ex.iter_rows(named=True)]))
    order = {"low": 0, "mid": 1, "high": 2}
    return ex.with_columns(pl.col("band").replace_strict(order).alias("_o")).sort("_o", "quality").drop("_o")


# ---- 3. review list -------------------------------------------------------------------

@dataclass
class ReviewRule:
    severity: str
    reason: str
    why: str


def review_table(tt: pl.DataFrame, events: pl.DataFrame, review_quality: float = 0.3,
                 prior_review_distance: float = 0.05) -> pl.DataFrame:
    """Tumours to check by hand, one row per reason, worst severity first."""
    flags_by_t = {t: ";".join(g["flags"].drop_nulls().to_list()) for (t,), g in events.group_by("tumour_id")} \
        if events.height else {}
    stacked_minor = set()
    if events.height and "source" in events.columns:
        m = events.filter(pl.col("flags").fill_null("").str.contains("minor_clone_fill"))
        stacked_minor = {t for (t, c), g in m.group_by("tumour_id", "clone") if g.height > 1}
    rows = []
    for r in tt.iter_rows(named=True):
        t, f = r["tumour_id"], json.loads(r["factors_json"])
        fl = flags_by_t.get(t, "")
        hits: list[tuple[ReviewRule, str]] = []
        if r["status"] == "failed":
            hits.append((ReviewRule("critical", "impossible solution (reversion)",
                                     "a clone below a GD sits in samples with fewer GDs; the tree, clone proportions "
                                     "or ploidy calls are wrong"), r["quality_code"] or ""))
        if r["status"] == "unresolved":
            hits.append((ReviewRule("critical", "unresolved", "no placement brings every sample to its PGDD GD "
                                    "level without over-doubling"), r["quality_code"] or ""))
        if r["forced_gds"]:
            hits.append((ReviewRule("serious", "forced GD (clone proportions vs ploidy)",
                                    "a doubled sample has no majority clone that can carry the GD; its clone "
                                    "proportions are probably off"), f.get("minority_coverage", [None, ""])[1]))
        if t in stacked_minor:
            hits.append((ReviewRule("serious", "several GDs on one minor clone",
                                    "a minor clone is more than one GD below its sample, i.e. cells from a less "
                                    "doubled branch; proportions or tree may be off"), r["gd_clones"] or ""))
        if "reversion" in f:
            hits.append((ReviewRule("serious", "ALPACA CN reversion", "ALPACA CN puts a clone below a GD back to "
                                    "diploid"), f["reversion"][1]))
        if "ploidy_path_mismatch" in f:
            hits.append((ReviewRule("serious", "GD path vs ploidy mismatch", "GDs on a path don't add up to the "
                                    "samples' GD count"), f["ploidy_path_mismatch"][1]))
        if "unplaced_events" in f:
            hits.append((ReviewRule("serious", "unplaced PGDD events", "PGDD events on clusters not on the ALPACA "
                                    "tree"), f["unplaced_events"][1]))
        if r["quality"] is not None and r["quality"] < review_quality and r["status"] == "resolved":
            hits.append((ReviewRule("serious", f"low quality (< {review_quality:g})",
                                    f"main penalty {r['main_penalty']}"), f"{r['quality']:.3f}"))
        if r["alpaca_parent"]:
            hits.append((ReviewRule("warning", "ALPACA CN to correct (parallel GDs)",
                                    "ALPACA put the parallel GDs' shared gain on their absent parent; halve its CN "
                                    "(TODO: post-hoc ALPACA correction)"),
                         ",".join(sorted(set(re.findall(r"alpaca_parent\(([^)]*)\)", fl))))))
        pp = r["ploidy_prior"] or ""
        moved = re.findall(r"([^\s,:]+) (\d)>(\d) \(([\d.]+);", pp) if pp.startswith("applied") else []
        far = [m for m in moved if float(m[3]) >= prior_review_distance]
        if far:
            hits.append((ReviewRule("warning", "ploidy call moved by the prior",
                                    f"a sample {prior_review_distance:g} or more from 0.5 was moved one GD level"),
                         ",".join(f"{m[0]} {m[1]}>{m[2]} ({m[3]})" for m in far)))
        if "ambiguous_solution" in f:
            hits.append((ReviewRule("warning", "ambiguous placement", "another placement is almost as good"),
                         f["ambiguous_solution"][1]))
        for rule, evidence in hits:
            rows.append({"tumour_id": t, "severity": rule.severity, "reason": rule.reason, "evidence": evidence,
                         "why": rule.why, "quality": r["quality"], "pattern": r["pattern"], "status": r["status"]})
    schema = {"tumour_id": pl.Utf8, "severity": pl.Utf8, "reason": pl.Utf8, "evidence": pl.Utf8, "why": pl.Utf8,
              "quality": pl.Float64, "pattern": pl.Utf8, "status": pl.Utf8}
    out = pl.DataFrame(rows, schema=schema)
    if out.is_empty():
        return out
    worst = out.group_by("tumour_id").agg(pl.col("severity").replace_strict(SEVERITY_ORDER).min().alias("_t"))
    return (out.join(worst, on="tumour_id")
               .with_columns(pl.col("severity").replace_strict(SEVERITY_ORDER).alias("_s"))
               .sort("_t", "quality", "tumour_id", "_s", nulls_last=True).drop("_t", "_s"))


# ---- figures -------------------------------------------------------------------------

def _style(ax, grid_axis: str | None = "x"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=9, length=0)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_facecolor(SURFACE)


def _fig(w, h):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(w, h), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    return fig, ax


def _save(fig, path: Path) -> str:
    fig.tight_layout()
    fig.savefig(path, facecolor=fig.get_facecolor())
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
    import matplotlib.pyplot as plt
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def fig_quality_by_pattern(tt: pl.DataFrame, path: Path) -> str:
    pats = [p for p in PATTERNS if p != "no GD" and tt.filter(pl.col("pattern") == p).height]
    fig, ax = _fig(7.5, 0.7 + 0.6 * max(len(pats), 1))
    rng = np.random.default_rng(1)
    for i, p in enumerate(pats):
        q = tt.filter((pl.col("pattern") == p) & pl.col("quality").is_not_null())["quality"].to_numpy()
        y = i + rng.uniform(-0.18, 0.18, len(q))
        ax.scatter(q, y, s=28, color=PATTERN_COLOUR[p], edgecolor=SURFACE, linewidth=1.2, zorder=3)
        if len(q):
            med = float(np.median(q))
            ax.plot([med, med], [i - 0.3, i + 0.3], color=INK, linewidth=2, zorder=4)
            ax.text(1.02, i, f"n={len(q)}  median {med:.2f}", va="center", fontsize=8, color=INK2,
                    transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(pats)), pats)
    ax.set_ylim(-0.6, len(pats) - 0.4)
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("tumour quality (product of evidence factors)", color=INK2, fontsize=9)
    ax.set_title("Quality per WGD pattern (dots: tumours; bar: median)", loc="left", fontsize=10, color=INK)
    _style(ax)
    fig.subplots_adjust(right=0.78)
    return _save(fig, path)


def fig_patterns(tt: pl.DataFrame, path: Path) -> str:
    counts = [(p, tt.filter(pl.col("pattern") == p).height) for p in PATTERNS]
    n = max(tt.height, 1)
    fig, ax = _fig(7.5, 2.8)
    for i, (p, c) in enumerate(counts):
        ax.barh(i, c, height=0.6, color=PATTERN_COLOUR[p])
        ax.text(c + n * 0.01, i, f"{c} ({c / n:.0%})", va="center", fontsize=9, color=INK)
    ax.set_yticks(range(len(counts)), [p for p, _ in counts])
    ax.invert_yaxis()
    ax.set_xlim(0, max(c for _, c in counts) * 1.25 + 1)
    ax.set_xlabel("tumours", color=INK2, fontsize=9)
    ax.set_title("Tumours per WGD pattern", loc="left", fontsize=10, color=INK)
    _style(ax)
    return _save(fig, path)


def fig_gd_grid(tt: pl.DataFrame, path: Path) -> str:
    """Clonal x subclonal GD counts: how many tumours have each combination."""
    rows, cols = [0, 1, 2, 3], [0, 1, 2, 3, 4]
    m = np.zeros((len(rows), len(cols)), int)
    for c, s in tt.select("clonal_gds", "subclonal_gds").iter_rows():
        m[min(c, 3), min(s, 4)] += 1
    fig, ax = _fig(6.0, 3.2)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("blue", ["#f0efec", "#86b6ef", "#2a78d6", "#104281"])
    ax.imshow(m, cmap=cmap, vmin=0, vmax=max(m.max(), 1), aspect="auto")
    for i in range(len(rows)):
        for j in range(len(cols)):
            if m[i, j]:
                ax.text(j, i, str(m[i, j]), ha="center", va="center", fontsize=10,
                        color="white" if m[i, j] > m.max() * 0.55 else INK)
    for k in range(len(cols) + 1):   # 2px surface gap between cells
        ax.axvline(k - 0.5, color=SURFACE, linewidth=2)
    for k in range(len(rows) + 1):
        ax.axhline(k - 0.5, color=SURFACE, linewidth=2)
    ax.set_xticks(range(len(cols)), ["0", "1", "2", "3", "4+"])
    ax.set_yticks(range(len(rows)), ["0", "1", "2", "3+"])
    ax.set_xlabel("subclonal GDs", color=INK2, fontsize=9)
    ax.set_ylabel("clonal GDs", color=INK2, fontsize=9)
    ax.set_title("Tumours by clonal and subclonal GD count", loc="left", fontsize=10, color=INK)
    _style(ax, grid_axis=None)
    return _save(fig, path)


def fig_sources(wgd_events: pl.DataFrame, path: Path) -> str:
    src = wgd_events.with_columns(
        pl.when(pl.col("flags").fill_null("").str.contains("minor_clone_fill")).then(pl.lit("fill_minor"))
          .otherwise(pl.col("source")).alias("s")).filter(pl.col("clone").is_not_null())
    order = ["clonal", "mutation", "fill", "fill_minor", "fill_forced"]
    counts = [(SOURCE_LABEL[k], src.filter(pl.col("s").replace({"mutation_multicluster": "mutation"}) == k).height)
              for k in order]
    fig, ax = _fig(7.5, 2.6)
    total = max(sum(c for _, c in counts), 1)
    for i, (lab, c) in enumerate(counts):
        ax.barh(i, c, height=0.6, color="#2a78d6")
        ax.text(c + total * 0.01, i, f"{c} ({c / total:.0%})", va="center", fontsize=9, color=INK)
    ax.set_yticks(range(len(counts)), [lab for lab, _ in counts])
    ax.invert_yaxis()
    ax.set_xlim(0, max(c for _, c in counts) * 1.25 + 1)
    ax.set_xlabel("GDs", color=INK2, fontsize=9)
    ax.set_title("What placed each GD", loc="left", fontsize=10, color=INK)
    _style(ax)
    return _save(fig, path)


def fig_penalties(tt: pl.DataFrame, path: Path) -> tuple[str, pl.DataFrame]:
    """Share of GD-carrying tumours where each quality factor fires, with its median value."""
    gd = tt.filter(pl.col("n_gds") > 0)
    seen: dict[str, list[float]] = {}
    for js in gd["factors_json"]:
        for k, (v, _) in json.loads(js).items():
            if v < 1:
                seen.setdefault(k, []).append(v)
    tab = pl.DataFrame({"factor": list(seen), "tumours": [len(v) for v in seen.values()],
                        "share": [len(v) / max(gd.height, 1) for v in seen.values()],
                        "median_value": [float(np.median(v)) for v in seen.values()]},
                       schema={"factor": pl.Utf8, "tumours": pl.Int64, "share": pl.Float64,
                               "median_value": pl.Float64}).sort("share", descending=True)
    fig, ax = _fig(7.5, 0.6 + 0.32 * max(tab.height, 1))
    for i, r in enumerate(tab.iter_rows(named=True)):
        ax.barh(i, r["share"], height=0.6, color="#2a78d6")
        ax.text(r["share"] + 0.01, i, f"{r['tumours']}  (median x{r['median_value']:.2f})", va="center",
                fontsize=8, color=INK2)
    ax.set_yticks(range(tab.height), tab["factor"].to_list())
    ax.invert_yaxis()
    ax.set_xlim(0, 1.3)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.xaxis.set_major_formatter(__import__("matplotlib.ticker", fromlist=["PercentFormatter"]).PercentFormatter(1.0))
    ax.set_xlabel("share of tumours with GDs where the factor lowers quality", color=INK2, fontsize=9)
    ax.set_title("Which evidence lowers quality", loc="left", fontsize=10, color=INK)
    _style(ax)
    return _save(fig, path), tab


def fig_review(rv: pl.DataFrame, path: Path) -> str:
    counts = (rv.group_by("severity", "reason").agg(pl.col("tumour_id").n_unique().alias("n"))
                .with_columns(pl.col("severity").replace_strict(SEVERITY_ORDER).alias("_s"))
                .sort("_s", "n", descending=[False, True]))
    fig, ax = _fig(7.5, 0.6 + 0.34 * max(counts.height, 1))
    for i, r in enumerate(counts.iter_rows(named=True)):
        ax.barh(i, r["n"], height=0.6, color=SEVERITY_COLOUR[r["severity"]])
        ax.text(r["n"] + 0.05, i, f"{r['n']}  {r['severity']}", va="center", fontsize=8, color=INK2)
    ax.set_yticks(range(counts.height), counts["reason"].to_list())
    ax.invert_yaxis()
    ax.set_xlim(0, (counts["n"].max() or 1) * 1.35 + 0.5)
    ax.set_xlabel("tumours", color=INK2, fontsize=9)
    ax.set_title("Reasons to review by hand", loc="left", fontsize=10, color=INK)
    _style(ax)
    return _save(fig, path)


# ---- report --------------------------------------------------------------------------

def _table_html(df: pl.DataFrame, cols: list[str], fmt: dict | None = None) -> str:
    fmt = fmt or {}
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = []
    for r in df.iter_rows(named=True):
        cells = []
        for c in cols:
            v = r.get(c)
            s = "" if v is None else (fmt[c](v) if c in fmt else str(v))
            cls = f' class="sev-{v}"' if c == "severity" else ""
            cells.append(f"<td{cls}>{html.escape(s)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


CSS = """
body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;background:#f9f9f7;color:#0b0b0b;max-width:1040px;margin:2em auto;padding:0 1em;line-height:1.45}
h1{font-size:1.5em}h2{font-size:1.2em;margin-top:2em;border-bottom:1px solid #e1e0d9;padding-bottom:.2em}
.tiles{display:flex;gap:12px;flex-wrap:wrap}.tile{background:#fcfcfb;border:1px solid #e1e0d9;border-radius:8px;padding:10px 14px;min-width:120px}
.tile .v{font-size:1.6em;font-variant-numeric:tabular-nums}.tile .l{color:#52514e;font-size:.85em}
img{max-width:100%;background:#fcfcfb;border:1px solid #e1e0d9;border-radius:6px;margin:.5em 0}
table{border-collapse:collapse;font-size:.85em;margin:.5em 0;width:100%}th,td{border-bottom:1px solid #e1e0d9;padding:4px 6px;text-align:left;vertical-align:top}
th{color:#52514e;font-weight:600}td{font-variant-numeric:tabular-nums}
.ex{background:#fcfcfb;border:1px solid #e1e0d9;border-left:4px solid #2a78d6;border-radius:6px;padding:8px 12px;margin:8px 0}
.ex.low{border-left-color:#d03b3b}.ex.mid{border-left-color:#fab219}.ex.high{border-left-color:#0ca30c}
.sev-critical::before{content:"\\25CF ";color:#d03b3b}.sev-serious::before{content:"\\25CF ";color:#ec835a}.sev-warning::before{content:"\\25CF ";color:#fab219}
.note{color:#52514e;font-size:.9em}
"""


def write_report(out: Path, tt: pl.DataFrame, ex: pl.DataFrame, rv: pl.DataFrame, pen: pl.DataFrame,
                 imgs: dict[str, str], title: str, review_quality: float) -> Path:
    n = tt.height
    gd = tt.filter(pl.col("n_gds") > 0).height
    cnt = {p: tt.filter(pl.col("pattern") == p).height for p in PATTERNS}
    sub_any = sum(cnt[p] for p in SUBCLONAL_PATTERNS)
    tiles = [("tumours", n), ("with a GD", f"{gd} ({gd / max(n, 1):.0%})"),
             ("clonal GD", f"{cnt['clonal only'] + cnt['clonal + subclonal']}"),
             ("subclonal GD", f"{sub_any}"), ("parallel GDs", f"{cnt['subclonal, parallel']}"),
             ("to review", f"{rv['tumour_id'].n_unique() if rv.height else 0}")]
    q_tab = (tt.filter(pl.col("quality").is_not_null() & (pl.col("pattern") != "no GD"))
               .group_by("pattern").agg(pl.len().alias("tumours"), pl.col("quality").median().alias("median"),
                                        pl.col("quality").quantile(0.25).alias("q25"),
                                        pl.col("quality").quantile(0.75).alias("q75"),
                                        (pl.col("quality") < review_quality).sum().alias(f"below {review_quality:g}"))
               .with_columns(pl.col("pattern").replace_strict({p: i for i, p in enumerate(PATTERNS)}).alias("_o"))
               .sort("_o").drop("_o"))
    f3 = lambda v: f"{v:.3f}"  # noqa: E731
    exs = "".join(f'<div class="ex {r["band"]}"><b>{html.escape(r["band"])}</b> — {html.escape(r["commentary"])}</div>'
                  for r in ex.iter_rows(named=True)) or '<p class="note">No tumour with subclonal or parallel GDs.</p>'
    img = lambda k, alt: f'<img alt="{html.escape(alt)}" src="data:image/png;base64,{imgs[k]}">'  # noqa: E731
    doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(title)}</title><style>{CSS}</style></head><body>
<h1>{html.escape(title)}</h1>
<div class="tiles">{''.join(f'<div class="tile"><div class="v">{v}</div><div class="l">{l}</div></div>' for l, v in tiles)}</div>
<p class="note">WGD pattern from the GDs per clone (MRCA = clonal). Quality is the product of evidence factors in [0, 1]
(wgd/quality.py, params.toml). Tables: tumour_summary.tsv, examples.tsv, quality_factors.tsv, review.tsv.</p>

<h2>1. Quality per WGD pattern</h2>
{img("quality", "Quality per WGD pattern")}
{_table_html(q_tab, q_tab.columns, {c: f3 for c in ("median", "q25", "q75")})}
<h3>Worked examples (subclonal or parallel GDs): 2 low, 2 mid, 2 high</h3>
{exs}

<h2>2. The cohort</h2>
{img("patterns", "Tumours per WGD pattern")}
{img("grid", "Tumours by clonal and subclonal GD count")}
{img("sources", "What placed each GD")}
{img("penalties", "Which evidence lowers quality")}
{_table_html(pen, pen.columns, {"share": lambda v: f"{v:.0%}", "median_value": f3})}

<h2>3. Tumours to review by hand</h2>
<p class="note">critical: the solution is impossible or incomplete; serious: the copy-number solution (clone proportions,
tree, ALPACA CN) is probably wrong or the score is below {review_quality:g}; warning: worth a look (ALPACA CN to correct,
ploidy calls moved by the prior, ambiguous placement).</p>
{img("review", "Reasons to review by hand") if rv.height else '<p class="note">Nothing to review.</p>'}
{_table_html(rv, ["severity", "tumour_id", "reason", "evidence", "why", "quality", "pattern"], {"quality": lambda v: f"{v:.3f}"})}
</body></html>"""
    path = out / "report.html"
    path.write_text(doc)
    return path


def run(out: Path, s04: Path, s05: Path, summary: pl.DataFrame | None = None, review_quality: float = 0.3,
        prior_review_distance: float = 0.05, title: str = "WGD calls: cohort summary") -> dict:
    read = lambda p: pl.read_csv(p, separator="\t", null_values=["NA"], infer_schema_length=None)  # noqa: E731
    tumours, events_placed = read(s04 / "tumours_placed.tsv"), read(s04 / "events_placed.tsv")
    calls, wgd_events = read(s05 / "wgd_calls.tsv"), read(s05 / "wgd_events.tsv")
    for c in ("flags", "event_factors"):
        if c not in wgd_events.columns:
            wgd_events = wgd_events.with_columns(pl.lit(None, pl.Utf8).alias(c))
    out.mkdir(parents=True, exist_ok=True)
    tt = tumour_table(tumours, calls, wgd_events, summary)
    ex = pick_examples(tt)
    rv = review_table(tt, wgd_events, review_quality, prior_review_distance)
    imgs = {"quality": fig_quality_by_pattern(tt, out / "quality_by_pattern.png"),
            "patterns": fig_patterns(tt, out / "patterns.png"),
            "grid": fig_gd_grid(tt, out / "gd_counts.png"),
            "sources": fig_sources(wgd_events, out / "gd_sources.png")}
    imgs["penalties"], pen = fig_penalties(tt, out / "quality_factors.png")
    if rv.height:
        imgs["review"] = fig_review(rv, out / "review_reasons.png")
    tt.drop("factors_json").write_csv(out / "tumour_summary.tsv", separator="\t", null_value="NA")
    ex.select("band", "tumour_id", "pattern", "quality", "gd_clones", "commentary").write_csv(
        out / "examples.tsv", separator="\t", null_value="NA")
    pen.write_csv(out / "quality_factors.tsv", separator="\t", null_value="NA")
    rv.write_csv(out / "review.tsv", separator="\t", null_value="NA")
    report = write_report(out, tt, ex, rv, pen, imgs, title, review_quality)
    return {"tumours": tt, "examples": ex, "review": rv, "factors": pen, "report": report}
