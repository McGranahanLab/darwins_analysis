"""Deviation of fractional copy number from the nearest integer.

All functions take/return polars frames. Inputs are tx_data canonical tables:
  wgs_consensus_segments  (sample x segment, MSAI-adjusted fractional CN)
  wgs_tree_clones         (single-clone = n_clones_tree == 1)
  wgs_battenberg_inventory (purity, ploidy, coverage per sample)

Per allele, `frac = cn - round(cn)` (signed, in [-0.5, 0.5]) and
`dev = |frac|` (0 = integer, 0.5 = half-way). With no integer signal
`dev` is uniform on [0, 0.5], mean 0.25. Genome summaries are weighted by
segment length.
"""

from __future__ import annotations

import numpy as np
import polars as pl

ALLELES = ("major", "minor")
RANDOM_DEV = 0.25          # mean |frac| if CN carried no integer signal
HIGH_DEV = 0.2             # threshold for "far from integer"
CHR_ORDER = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24}


# ---------------------------------------------------------------- loading
def load(tx_dir, include_x: bool = False, min_seg_bp: int = 0):
    """Return (segments of single-clone tumours, per-sample metadata, tumour list)."""
    seg = pl.read_parquet(tx_dir / "wgs_consensus_segments.parquet")
    tree = pl.read_parquet(tx_dir / "wgs_tree_clones.parquet")
    inv = pl.read_parquet(tx_dir / "wgs_battenberg_inventory.parquet")

    clones = tree.group_by("tumour_id").agg(pl.col("n_clones_tree").first())
    single = clones.filter(pl.col("n_clones_tree") == 1)["tumour_id"]
    meta = inv.select("tumour_id", "patient_id", "sample_name_hash", "region",
                      "purity", "ploidy", "coverage")

    seg = seg.filter(pl.col("tumour_id").is_in(single.implode()))
    if not include_x:
        seg = seg.filter(pl.col("chr").is_in([str(i) for i in range(1, 23)]))
    seg = seg.with_columns((pl.col("end") - pl.col("start") + 1).alias("bp"))
    if min_seg_bp:
        seg = seg.filter(pl.col("bp") >= min_seg_bp)
    seg = seg.join(meta.drop("tumour_id", "patient_id"), on="sample_name_hash", how="left")
    return seg, meta, clones


def add_deviation(seg: pl.DataFrame) -> pl.DataFrame:
    out = []
    for a in ALLELES:
        cn = pl.col(f"cn_{a}")
        out += [cn.round().cast(pl.Int64).alias(f"int_{a}"),
                (cn - cn.round()).alias(f"frac_{a}"),
                (cn - cn.round()).abs().alias(f"dev_{a}")]
    seg = seg.with_columns(out)
    return seg.with_columns(((pl.col("dev_major") + pl.col("dev_minor")) / 2).alias("dev"))


def long_alleles(seg: pl.DataFrame) -> pl.DataFrame:
    """One row per (sample, segment, allele)."""
    keep = ["tumour_id", "sample", "chr", "start", "end", "bp", "purity", "ploidy"]
    return pl.concat([
        seg.select(*keep, pl.lit(a).alias("allele"), pl.col(f"cn_{a}").alias("cn"),
                   pl.col(f"int_{a}").alias("int"), pl.col(f"frac_{a}").alias("frac"),
                   pl.col(f"dev_{a}").alias("dev"))
        for a in ALLELES
    ])


def _wmean(x: str, w: str = "bp") -> pl.Expr:
    return (pl.col(x) * pl.col(w)).sum() / pl.col(w).sum()


# ---------------------------------------------------------------- per sample
def per_sample(seg: pl.DataFrame) -> pl.DataFrame:
    """Genome-wide deviation summary per sample, plus a scale-error fit.

    `scale_slope`/`scale_intercept`: weighted least-squares fit of observed CN
    on its nearest integer (both alleles). A wrong purity/ploidy fit stretches
    or shifts every state, giving slope != 1 or intercept != 0, while noise
    alone doesn't. This is only a diagnostic: rounding pulls the fit towards 1.
    """
    rows = []
    for (tid, smp), g in seg.group_by(["tumour_id", "sample"], maintain_order=True):
        w = np.concatenate([g["bp"].to_numpy()] * 2).astype(float)
        x = np.concatenate([g["int_major"].to_numpy(), g["int_minor"].to_numpy()]).astype(float)
        y = np.concatenate([g["cn_major"].to_numpy(), g["cn_minor"].to_numpy()])
        slope, intercept = _wls(x, y, w)
        rows.append({"tumour_id": tid, "sample": smp, "scale_slope": slope,
                     "scale_intercept": intercept})
    fit = pl.DataFrame(rows, schema={"tumour_id": pl.String, "sample": pl.String,
                                     "scale_slope": pl.Float64, "scale_intercept": pl.Float64})
    summ = seg.group_by("tumour_id", "sample").agg(
        pl.col("patient_id").first(),
        pl.col("purity").first(), pl.col("ploidy").first(), pl.col("coverage").first(),
        pl.len().alias("n_segments"),
        pl.col("bp").sum().alias("genome_bp"),
        _wmean("dev").alias("dev_mean"),
        _wmean("dev_major").alias("dev_major"),
        _wmean("dev_minor").alias("dev_minor"),
        _wmean("frac_major").alias("frac_major_signed"),
        _wmean("frac_minor").alias("frac_minor_signed"),
        ((pl.col("dev") > HIGH_DEV).cast(pl.Float64) * pl.col("bp")).sum().truediv(pl.col("bp").sum())
        .alias(f"frac_genome_dev_gt_{HIGH_DEV}"),
        pl.col("dev").median().alias("dev_median_unweighted"),
    )
    return summ.join(fit, on=["tumour_id", "sample"]).sort("dev_mean", descending=True)


def _wls(x, y, w):
    if np.unique(x).size < 2:
        return float("nan"), float("nan")
    W = w / w.sum()
    mx, my = (W * x).sum(), (W * y).sum()
    slope = (W * (x - mx) * (y - my)).sum() / (W * (x - mx) ** 2).sum()
    return float(slope), float(my - slope * mx)


# ---------------------------------------------------------------- across tumours
def spearman(x, y, n_perm: int = 10000, seed: int = 1):
    """Spearman rho with a two-sided permutation p-value (no scipy needed)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = x.size
    if n < 4:
        return float("nan"), float("nan"), n
    rx, ry = _rank(x), _rank(y)
    rho = np.corrcoef(rx, ry)[0, 1]
    rng = np.random.default_rng(seed)
    perm = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(n_perm)])
    p = (np.sum(np.abs(perm) >= abs(rho) - 1e-12) + 1) / (n_perm + 1)
    return float(rho), float(p), n


def _rank(a):
    order = a.argsort(kind="mergesort")
    r = np.empty(a.size)
    r[order] = np.arange(a.size)
    for v in np.unique(a):  # average ties
        m = a == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def per_tumour(samples: pl.DataFrame) -> pl.DataFrame:
    """One row per tumour (mean over its samples), to avoid counting multi-region tumours twice."""
    return samples.group_by("tumour_id").agg(
        pl.len().alias("n_samples"),
        *[pl.col(c).mean() for c in ("purity", "ploidy", "coverage", "dev_mean", "dev_major",
                                     "dev_minor", f"frac_genome_dev_gt_{HIGH_DEV}",
                                     "scale_slope", "scale_intercept")],
    ).sort("tumour_id")


def correlations(tumours: pl.DataFrame, n_perm: int = 10000) -> pl.DataFrame:
    y = tumours["dev_mean"].to_numpy()
    rows = []
    for cov in ("purity", "ploidy", "coverage"):
        rho, p, n = spearman(tumours[cov].to_numpy(), y, n_perm)
        rows.append({"outcome": "dev_mean", "covariate": cov, "spearman_rho": rho,
                     "p_perm": p, "n_tumours": n})
    return pl.DataFrame(rows)


# ---------------------------------------------------------------- within tumours
def by_state(long: pl.DataFrame) -> pl.DataFrame:
    """Deviation by allele and nearest-integer state (states >= 6 pooled)."""
    return (long.with_columns(pl.col("int").clip(upper_bound=6).alias("state"))
            .group_by("allele", "state").agg(
                pl.col("tumour_id").n_unique().alias("n_tumours"),
                pl.len().alias("n_segments"),
                (pl.col("bp").sum() / 1e6).alias("Mb"),
                _wmean("dev").alias("dev_mean"),
                _wmean("frac").alias("frac_signed_mean"))
            .sort("allele", "state"))


def by_length(seg: pl.DataFrame, n_bins: int = 5) -> pl.DataFrame:
    """Deviation by segment-length quantile (short segments are noisier)."""
    qs = np.quantile(seg["bp"].to_numpy(), np.linspace(0, 1, n_bins + 1)[1:-1])
    return (seg.with_columns(pl.col("bp").cut(list(qs)).alias("length_bin"))
            .group_by("length_bin").agg(pl.len().alias("n_segments"),
                                        pl.col("bp").median().alias("median_bp"),
                                        pl.col("dev").mean().alias("dev_mean_unweighted"))
            .sort("median_bp"))


def by_chromosome(seg: pl.DataFrame) -> pl.DataFrame:
    """Per-chromosome deviation relative to each sample's genome-wide mean.

    ratio > 1: this chromosome is further from integer than the rest of the
    sample. Recurrent if it holds in many tumours.
    """
    smp = seg.group_by("tumour_id", "sample").agg(_wmean("dev").alias("smp_dev"))
    chrom = (seg.group_by("tumour_id", "sample", "chr").agg(_wmean("dev").alias("chr_dev"))
             .join(smp, on=["tumour_id", "sample"])
             .with_columns((pl.col("chr_dev") / pl.col("smp_dev")).alias("ratio")))
    tum = chrom.group_by("tumour_id", "chr").agg(pl.col("ratio").mean())
    return (tum.group_by("chr").agg(
                pl.len().alias("n_tumours"),
                pl.col("ratio").median().alias("median_ratio"),
                (pl.col("ratio") > 1.5).mean().alias("frac_tumours_ratio_gt_1.5"))
            .with_columns(pl.col("chr").replace_strict(CHR_ORDER, default=99).alias("_o"))
            .sort("_o").drop("_o"))


def _segment_z(seg: pl.DataFrame) -> pl.DataFrame:
    """Per-segment deviation z-scored within its sample (length-weighted mean/sd).

    Removes each sample's genome-wide level (driven by purity), leaving where
    in its own genome a sample is furthest from integer.
    """
    mu = _wmean("dev").over("sample")
    var = ((pl.col("dev") - mu) ** 2 * pl.col("bp")).sum().over("sample") / pl.col("bp").sum().over("sample")
    return seg.with_columns(((pl.col("dev") - mu) / var.sqrt()).alias("z"))


def cohort_pieces(seg: pl.DataFrame) -> pl.DataFrame:
    """Cut the genome at every segment breakpoint of every tumour.

    Each piece lies inside at most one segment per sample, so it takes that
    segment's value unchanged: no bin size, no averaging across segments.
    Returns one row per (piece, tumour), the tumour's value being the mean z
    over its samples (multi-region tumours count once).
    """
    z = _segment_z(seg)
    bps = pl.concat([z.select("chr", pl.col("start").alias("p")),
                     z.select("chr", (pl.col("end") + 1).alias("p"))]).unique().sort("chr", "p")
    pieces = (bps.with_columns(pl.col("p").shift(-1).over("chr").alias("next"))
              .drop_nulls("next")
              .select("chr", pl.col("p").alias("start"), (pl.col("next") - 1).alias("end")))
    out = []
    for smp, g in z.group_by("sample"):
        g = g.select("tumour_id", "chr", pl.col("start").alias("seg_start"),
                     pl.col("end").alias("seg_end"), "z", "dev").sort("seg_start")
        h = (pieces.sort("start").join_asof(g, left_on="start", right_on="seg_start", by="chr",
                                            strategy="backward",
                                            check_sortedness=False)
             .filter(pl.col("seg_end").is_not_null() & (pl.col("end") <= pl.col("seg_end"))))
        out.append(h.select("chr", "start", "end", "tumour_id", "z", "dev"))
    return (pl.concat(out).group_by("chr", "start", "end", "tumour_id")
            .agg(pl.col("z").mean(), pl.col("dev").mean()))


def recurrent_pieces(pt: pl.DataFrame, min_frac: float = 0.5, n_perm: int = 1000,
                     seed: int = 1) -> pl.DataFrame:
    """Across tumours, per piece: mean z, fraction of tumours with z > 1, and a
    genome-wide (FWER) permutation p for mean z.

    Null: each tumour's z profile is rotated by an independent random offset
    along the concatenated genome, keeping its segment structure and value
    distribution but breaking positional alignment between tumours. p_fwer =
    share of permutations whose maximum mean z (over pieces covered by
    >= min_frac of tumours) reaches the observed value.
    """
    n_t = pt["tumour_id"].n_unique()
    summ = (pt.group_by("chr", "start", "end").agg(
                pl.len().alias("n_tumours"), pl.col("z").mean().alias("mean_z"),
                (pl.col("z") > 1).mean().alias("frac_tumours_z_gt_1"),
                pl.col("dev").median().alias("median_dev"))
            .with_columns((pl.col("n_tumours") / n_t).alias("coverage_frac"),
                          (pl.col("end") - pl.col("start") + 1).alias("bp"),
                          pl.col("chr").replace_strict(CHR_ORDER, default=99).alias("_o"))
            .sort("_o", "start"))
    # Concatenated-genome coordinates.
    ends = summ.group_by("_o").agg(pl.col("end").max()).sort("_o")
    offs = dict(zip(ends["_o"], np.concatenate([[0], np.cumsum(ends["end"].to_numpy())[:-1]])))
    G = float(ends["end"].sum())
    mid = (summ["_o"].replace_strict(offs).to_numpy()
           + (summ["start"].to_numpy() + summ["end"].to_numpy()) / 2)
    keep = summ["coverage_frac"].to_numpy() >= min_frac

    profiles = []
    for _, g in pt.with_columns(pl.col("chr").replace_strict(CHR_ORDER, default=99).alias("_o")) \
                  .group_by("tumour_id"):
        g = g.sort("_o", "start")
        a = g["_o"].replace_strict(offs).to_numpy() + g["start"].to_numpy()
        b = g["_o"].replace_strict(offs).to_numpy() + g["end"].to_numpy()
        profiles.append((a.astype(float), b.astype(float), g["z"].to_numpy()))

    def mean_z_at(pos_shift):
        tot = np.zeros(mid.size); cnt = np.zeros(mid.size)
        for (a, b, zz), sh in zip(profiles, pos_shift):
            q = (mid - sh) % G
            i = np.searchsorted(a, q, side="right") - 1
            ok = (i >= 0) & (q <= b[np.clip(i, 0, None)])
            tot[ok] += zz[i[ok]]; cnt[ok] += 1
        with np.errstate(invalid="ignore"):
            return tot / cnt, cnt

    rng = np.random.default_rng(seed)
    obs = summ["mean_z"].to_numpy()
    null_max = np.empty(n_perm)
    for k in range(n_perm):
        m, c = mean_z_at(rng.uniform(0, G, len(profiles)))
        ok = keep & (c / n_t >= min_frac)
        null_max[k] = np.nanmax(m[ok]) if ok.any() else -np.inf
    p = np.array([(np.sum(null_max >= o) + 1) / (n_perm + 1) if kp else np.nan
                  for o, kp in zip(obs, keep)])
    return summ.with_columns(pl.Series("p_fwer", p)).drop("_o")


def recurrent_regions(pieces: pl.DataFrame, z_min: float = 1.0) -> pl.DataFrame:
    """Merge adjacent pieces with mean z >= z_min into regions, for reporting."""
    s = pieces.filter(pl.col("p_fwer").is_not_null()).with_columns(
        (pl.col("mean_z") >= z_min).alias("hi"))
    s = s.with_columns(((pl.col("hi") != pl.col("hi").shift(1).over("chr"))
                        | (pl.col("start") != pl.col("end").shift(1).over("chr") + 1))
                       .fill_null(True).cum_sum().alias("rid"))
    return (s.filter("hi").group_by("chr", "rid").agg(
                pl.col("start").min(), pl.col("end").max(), pl.len().alias("n_pieces"),
                pl.col("n_tumours").min().alias("min_tumours"),
                ((pl.col("mean_z") * pl.col("bp")).sum() / pl.col("bp").sum()).alias("mean_z"),
                pl.col("mean_z").max().alias("max_z"), pl.col("p_fwer").min().alias("min_p_fwer"))
            .with_columns(((pl.col("end") - pl.col("start") + 1) / 1e6).alias("Mb"))
            .sort("min_p_fwer", pl.col("max_z"), descending=[False, True]).drop("rid"))


def region_concordance(seg: pl.DataFrame) -> pl.DataFrame:
    """Multi-region single-clone tumours: do the same segments deviate in every region?

    Consensus segmentation gives each region identical segments, so per-segment
    deviation can be correlated directly between regions (Spearman).
    """
    rows = []
    for tid, g in seg.group_by("tumour_id"):
        samples = sorted(g["sample"].unique())
        if len(samples) < 2:
            continue
        wide = g.pivot(on="sample", index=["chr", "start"], values="dev")
        for i, a in enumerate(samples):
            for b in samples[i + 1:]:
                rho, p, n = spearman(wide[a].to_numpy(), wide[b].to_numpy(), n_perm=2000)
                rows.append({"tumour_id": tid[0], "sample_a": a, "sample_b": b,
                             "n_segments": n, "spearman_rho": rho, "p_perm": p})
    return pl.DataFrame(rows, schema={"tumour_id": pl.String, "sample_a": pl.String,
                                      "sample_b": pl.String, "n_segments": pl.Int64,
                                      "spearman_rho": pl.Float64, "p_perm": pl.Float64})
