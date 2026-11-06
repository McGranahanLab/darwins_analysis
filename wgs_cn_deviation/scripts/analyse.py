"""Run the single-clone CN deviation analysis on built tx_data tables.

Usage:
  python scripts/analyse.py --tx-data-dir <tx_data>/data [--out output/<name>]
                            [--include-x] [--min-seg-bp N] [--min-frac 0.5]

Writes CSV tables, PNG figures and summary.md into --out.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from wgscn import deviation as d  # noqa: E402


def _git(path) -> str:
    try:
        return subprocess.run(["git", "-C", str(path), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def figures(out, samples, tumours, state, pieces, seg):
    # 1. deviation vs purity
    fig, ax = plt.subplots(figsize=(5, 4))
    sc = ax.scatter(tumours["purity"], tumours["dev_mean"], c=tumours["ploidy"], cmap="viridis",
                    s=30, edgecolor="k", linewidth=0.3)
    ax.axhline(d.RANDOM_DEV, ls="--", c="grey", lw=0.8)
    ax.text(ax.get_xlim()[1], d.RANDOM_DEV, " no integer signal", va="center", ha="right",
            fontsize=7, color="grey", backgroundcolor="white")
    ax.set_xlabel("purity"); ax.set_ylabel("mean |CN - nearest integer|")
    fig.colorbar(sc, label="ploidy"); ax.set_title("Single-clone tumours")
    fig.tight_layout(); fig.savefig(out / "fig1_deviation_vs_purity.png", dpi=150); plt.close(fig)

    # 2. deviation by integer state
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for a, m in (("major", "o"), ("minor", "s")):
        s = state.filter(pl.col("allele") == a)
        ax.plot(s["state"], s["dev_mean"], marker=m, label=a)
    ax.set_xlabel("nearest integer (6 = 6+)"); ax.set_ylabel("mean deviation")
    ax.legend(title="allele"); fig.tight_layout()
    fig.savefig(out / "fig2_deviation_by_state.png", dpi=150); plt.close(fig)

    # 3. recurrent regions along the genome
    fig, ax = plt.subplots(figsize=(10, 3))
    off, ticks, x0 = {}, [], 0
    for c in pieces["chr"].unique(maintain_order=True):
        off[c] = x0
        e = pieces.filter(pl.col("chr") == c)["end"].max()
        ticks.append((x0 + e / 2, c)); x0 += e
        ax.axvline(x0, c="lightgrey", lw=0.5)
    for r in pieces.iter_rows(named=True):
        sig = r["p_fwer"] is not None and r["p_fwer"] < 0.05
        ax.plot([off[r["chr"]] + r["start"], off[r["chr"]] + r["end"]], [r["mean_z"]] * 2,
                c="C3" if sig else ("C0" if r["coverage_frac"] >= 0.5 else "lightgrey"), lw=1.5)
    ax.set_xticks([t for t, _ in ticks], [c for _, c in ticks], fontsize=6)
    ax.axhline(0, c="k", lw=0.5)
    ax.set_ylabel("mean within-sample z")
    ax.set_title("Deviation along the genome across tumours (red: p_fwer < 0.05; grey: <50% of tumours)")
    fig.tight_layout(); fig.savefig(out / "fig3_recurrent_regions.png", dpi=150); plt.close(fig)

    # 4. signed fractional part by purity tertile
    fig, ax = plt.subplots(figsize=(5, 3.5))
    q = np.quantile(tumours["purity"].to_numpy(), [1 / 3, 2 / 3])
    long = d.long_alleles(seg)
    for lab, lo, hi in (("low", -1, q[0]), ("mid", q[0], q[1]), ("high", q[1], 2)):
        f = long.filter((pl.col("purity") > lo) & (pl.col("purity") <= hi))["frac"].to_numpy()
        if f.size:
            ax.hist(f, bins=40, range=(-0.5, 0.5), histtype="step", density=True,
                    label=f"{lab} purity")
    ax.set_xlabel("CN - nearest integer"); ax.set_ylabel("density"); ax.legend()
    fig.tight_layout(); fig.savefig(out / "fig4_signed_fraction.png", dpi=150); plt.close(fig)


def cn_profiles(out, seg, y_max: float = 8.0):
    """Raw CN along the genome: one figure per (tumour, sample), one panel per allele.

    Segments with |CN - nearest integer| > HIGH_DEV are red. Values above
    y_max are drawn at y_max with a triangle. Samples ordered by purity.
    Writes fig5_cn_profiles/<tumour>__<sample>.png and all pages in
    fig5_cn_profiles.pdf.
    """
    pdir = out / "fig5_cn_profiles"
    pdir.mkdir(exist_ok=True)
    chrs = sorted(seg["chr"].unique(), key=lambda c: d.CHR_ORDER.get(c, 99))
    ends = dict(seg.group_by("chr").agg(pl.col("end").max()).rows())
    off, x0 = {}, 0
    for c in chrs:
        off[c] = x0; x0 += ends[c]
    order = (seg.group_by("tumour_id", "sample")
             .agg(pl.col("purity").first(), pl.col("ploidy").first(),
                  pl.col("msai_adjusted").first() if "msai_adjusted" in seg.columns
                  else pl.lit(None).alias("msai_adjusted"))
             .sort("purity", "tumour_id", "sample", nulls_last=True))
    legend = [Line2D([], [], c="C0", lw=2, label=f"|CN - integer| <= {d.HIGH_DEV}"),
              Line2D([], [], c="C3", lw=2, label=f"|CN - integer| > {d.HIGH_DEV}")]
    with PdfPages(out / "fig5_cn_profiles.pdf") as pdf:
        for r in order.iter_rows(named=True):
            g = seg.filter((pl.col("tumour_id") == r["tumour_id"]) & (pl.col("sample") == r["sample"]))
            fig, axes = plt.subplots(2, 1, figsize=(12, 4.5), sharex=True)
            for ax, allele in zip(axes, d.ALLELES):
                for row in g.select("chr", "start", "end", f"cn_{allele}", f"dev_{allele}").iter_rows():
                    c, a, b, cn, dev = row
                    xa, xb = off[c] + a, off[c] + b
                    col = "C3" if dev > d.HIGH_DEV else "C0"
                    ax.plot([xa, xb], [min(cn, y_max)] * 2, c=col, lw=2, solid_capstyle="butt")
                    if cn > y_max:
                        ax.plot((xa + xb) / 2, y_max, marker="^", c=col, ms=4)
                top = min(y_max, max(3.0, np.ceil(g[f"cn_{allele}"].max() or 0) + 0.5))
                for k in range(int(top) + 1):
                    ax.axhline(k, c="lightgrey", lw=0.5, ls=":", zorder=0)
                for c in chrs[1:]:
                    ax.axvline(off[c], c="lightgrey", lw=0.5, zorder=0)
                ax.set_ylim(-0.2, top + 0.2)
                ax.set_ylabel(f"{allele} CN")
            axes[-1].set_xticks([off[c] + ends[c] / 2 for c in chrs], chrs, fontsize=7)
            axes[-1].set_xlim(0, x0)
            adj = {True: "MSAI-adjusted", False: "MSAI not adjusted"}.get(r["msai_adjusted"], "")
            pur = "NA" if r["purity"] is None else f"{r['purity']:.2f}"
            plo = "NA" if r["ploidy"] is None else f"{r['ploidy']:.2f}"
            fig.suptitle(f"{r['tumour_id']}  |  {r['sample']}  |  purity {pur}  |  ploidy {plo}  |  {adj}",
                         fontsize=10)
            axes[0].legend(handles=legend, fontsize=7, loc="upper right", ncol=2, framealpha=0.9)
            fig.tight_layout()
            fig.savefig(pdir / f"{r['tumour_id']}__{r['sample']}.png".replace("/", "_"), dpi=120)
            pdf.savefig(fig)
            plt.close(fig)
    return order.height


def _md(df: pl.DataFrame, n: int = 15, digits: int = 3) -> str:
    df = df.head(n).with_columns(pl.col(pl.Float64).round(digits))
    cols = df.columns
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join("" if v is None else str(v) for v in r) + " |" for r in df.rows()]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tx-data-dir", type=pathlib.Path,
                    default=pathlib.Path(os.environ.get("TX_DATA_DIR", ROOT.parent / "tx_data" / "data")))
    ap.add_argument("--out", type=pathlib.Path)
    ap.add_argument("--include-x", action="store_true")
    ap.add_argument("--min-seg-bp", type=int, default=0)
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="test only pieces covered by at least this fraction of tumours")
    ap.add_argument("--n-perm-regions", type=int, default=1000)
    ap.add_argument("--n-perm", type=int, default=10000)
    a = ap.parse_args()
    out = a.out or ROOT / "output" / dt.date.today().isoformat()
    out.mkdir(parents=True, exist_ok=True)

    seg, meta, clones = d.load(a.tx_data_dir, a.include_x, a.min_seg_bp)
    single = clones.filter(pl.col("n_clones_tree") == 1)["tumour_id"]
    no_seg = sorted(set(single) - set(seg["tumour_id"].unique()))
    if no_seg:
        print(f"warning: {len(no_seg)}/{single.len()} single-clone tumour(s) have no rows in "
              f"wgs_consensus_segments (check the tx_data build log): {no_seg}", file=sys.stderr)
    if seg.is_empty():
        sys.exit(f"no segments for any of the {single.len()} single-clone tumours on the tree")
    missing_purity = seg.filter(pl.col("purity").is_null())["sample"].unique().to_list()
    seg = d.add_deviation(seg)

    samples = d.per_sample(seg)
    tumours = d.per_tumour(samples)
    corr = d.correlations(tumours, a.n_perm)
    state = d.by_state(d.long_alleles(seg))
    length = d.by_length(seg)
    chrom = d.by_chromosome(seg)
    pieces = d.recurrent_pieces(d.cohort_pieces(seg), a.min_frac, a.n_perm_regions)
    regions = d.recurrent_regions(pieces)
    conc = d.region_concordance(seg)

    for name, df in (("per_sample", samples), ("per_tumour", tumours), ("correlations", corr),
                     ("by_state", state), ("by_length", length), ("by_chromosome", chrom),
                     ("recurrent_pieces", pieces), ("recurrent_regions", regions),
                     ("region_concordance", conc)):
        df.with_columns(pl.col(pl.Categorical).cast(pl.String)).write_csv(out / f"{name}.csv")
    figures(out, samples, tumours, state, pieces, seg)
    n_profiles = cn_profiles(out, seg)

    n_all = clones.height
    n_single = clones.filter(pl.col("n_clones_tree") == 1).height
    run = {"date": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "tx_data_dir": str(a.tx_data_dir), "tx_data_commit": _git(a.tx_data_dir.parent),
           "analysis_commit": _git(ROOT), "include_x": a.include_x,
           "min_seg_bp": a.min_seg_bp, "min_frac": a.min_frac, "n_perm": a.n_perm,
           "n_perm_regions": a.n_perm_regions}
    (out / "run.json").write_text(json.dumps(run, indent=1) + "\n")

    md = [
        "# Single-clone WGS tumours: deviation of CN from the nearest integer", "",
        f"Run: {run['date']}; tx_data `{run['tx_data_commit']}`, analysis `{run['analysis_commit']}`. "
        f"Chromosomes: {'1-22, X' if a.include_x else '1-22'}; min segment {a.min_seg_bp} bp.", "",
        f"- Tumours with a tree: {n_all}; single-clone (`n_clones_tree = 1`): {n_single}; "
        f"with segments: {tumours.height} ({samples.height} samples).",
        f"- Deviation = |CN - round(CN)| per allele, length-weighted; {d.RANDOM_DEV} = no integer signal.",
        f"- Cohort median of per-tumour mean deviation: {tumours['dev_mean'].median():.3f} "
        f"(range {tumours['dev_mean'].min():.3f}-{tumours['dev_mean'].max():.3f}).",
    ]
    if missing_purity:
        md.append(f"- Samples without inventory purity (dropped from correlations): {missing_purity[:10]}")
    md.append(f"- MSAI-adjusted segments: {tumours['msai_adjusted'].sum()} tumours; "
              f"NOT_adjusted (single-sample, adjustment can't run): "
              f"{(~tumours['msai_adjusted']).sum()}.")
    if no_seg:
        md.append(f"- Single-clone tumours without segments (not analysed): {len(no_seg)}: {no_seg}")
    md += ["", "## Across tumours (one row per tumour, Spearman, permutation p)", "", _md(corr), "",
           "![](fig1_deviation_vs_purity.png)", "",
           "## Tumours, highest deviation first", "",
           _md(tumours.sort("dev_mean", descending=True).select(
               "tumour_id", "n_samples", "msai_adjusted", "purity", "ploidy", "coverage", "dev_mean",
               "dev_major", "dev_minor", f"frac_genome_dev_gt_{d.HIGH_DEV}",
               "scale_slope", "scale_intercept"), n=25), "",
           "`scale_slope`/`scale_intercept`: fit of observed CN on its nearest integer. Values away "
           "from 1/0 point to a purity/ploidy misfit rather than noise.", "",
           "## Within tumours", "", "### By allele and integer state", "", _md(state, n=20), "",
           "![](fig2_deviation_by_state.png)", "", "### By segment length", "", _md(length), "",
           "### By chromosome (ratio to the sample's own mean)", "", _md(chrom, n=24), "",
           "### Recurrent regions", "",
           f"Genome cut at every breakpoint of every tumour ({pieces.height} pieces); each piece "
           "takes its segment's deviation, z-scored within the sample. Regions = adjacent pieces "
           f"with mean z >= 1, tested where >= {a.min_frac:.0%} of tumours have a segment. "
           "p_fwer: genome-wide, rotating each tumour's profile independently.", "",
           _md(regions) if regions.height else "No piece reaches mean z >= 1.", "",
           "![](fig3_recurrent_regions.png)", "",
           "### Multi-region single-clone tumours: region-to-region concordance", "",
           _md(conc) if conc.height else "No multi-region single-clone tumours.", "",
           "### Signed fractional part", "", "![](fig4_signed_fraction.png)", "",
           "## Raw copy-number profiles", "",
           f"{n_profiles} samples, ordered by purity: `fig5_cn_profiles.pdf` (one page each) and "
           f"`fig5_cn_profiles/<tumour>__<sample>.png`. Red = |CN - nearest integer| > {d.HIGH_DEV}; "
           "dotted lines = integers; triangles = CN above the axis limit (8).", ""]
    (out / "summary.md").write_text("\n".join(md))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
