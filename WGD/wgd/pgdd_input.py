"""Build the PGDD input table from CONIPHER outputs and MSAI-adjusted CN segments.

Python port of `scripts/00_prepare_pgdd_input.R` (Piotr, 2026-11-09). Per tumour:

  <conipher_dir>/<tumour>/<tumour>.SCoutput.CLEAN.tsv   mutations (MUT_COPY, CN at locus, PLOIDY)
  <conipher_dir>/<tumour>/treeTable.tsv                  mutations that made it onto the tree (treeCLUSTER)
  <conipher_dir>/<tumour>/allTrees.txt                   tree edges (parent \\t child)
  <msai_dir>/<tumour>/<tumour>_copynumber_MSAI_adjusted_segments_tree_input.tsv   CN segments per sample

Mutations are restricted to those in treeTable.tsv and take their tree cluster.
`num_gds` per sample is called from the segments: the highest k in 3, 2, 1 with
>= 50% of the genome (bp) at round(MajCN) >= k + 1, else 0, where
MajCN = max(COPY_NUMBER_A, COPY_NUMBER_B). `frac_maj_ge2` is the fraction of
the genome that looks doubled.

Differences from the R script:
* allTrees.txt: only the default (first) tree is read. The R script strips the
  `# tree N` headers and pools the edges of all alternative trees.
* Tumours with a single-cluster tree are skipped, as in R (reported in `skipped`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

MUT_COLS = ("CASE_ID", "CHR", "POS", "REF", "ALT", "CLUSTER", "SAMPLE",
            "MUT_COPY", "COPY_NUMBER_A", "COPY_NUMBER_B", "mutation_id")
CN_COLS = ("SAMPLE", "STARTPOS", "ENDPOS", "COPY_NUMBER_A", "COPY_NUMBER_B")

# Column order of pgdd_input_all_tumours.tsv, as written by the R script.
OUTPUT_COLS = ("tumour_id", "chromosome", "position", "ref", "alt", "cluster_id",
               "is_clonal_cluster", "sample_id", "mut_cpn", "MajCN", "num_gds",
               "frac_maj_ge2", "frac_maj_ge3", "frac_maj_ge4",
               "conf_gd1", "conf_gd2", "conf_gd3", "ploidy")


class SkipTumour(Exception):
    """Input problem that makes one tumour unusable; the run continues."""


@dataclass
class TumourFiles:
    tumour_id: str
    mutations: Path
    tree_table: Path
    all_trees: Path
    cn: Path

    @classmethod
    def locate(cls, conipher_dir: Path, msai_dir: Path, tumour_id: str) -> "TumourFiles":
        t = Path(conipher_dir) / tumour_id
        return cls(
            tumour_id=tumour_id,
            mutations=t / f"{tumour_id}.SCoutput.CLEAN.tsv",
            tree_table=t / "treeTable.tsv",
            all_trees=t / "allTrees.txt",
            cn=Path(msai_dir) / tumour_id / f"{tumour_id}_copynumber_MSAI_adjusted_segments_tree_input.tsv",
        )

    def check(self) -> None:
        for label, p in (("mutation file", self.mutations), ("CN segment file", self.cn),
                         ("tree file", self.all_trees), ("tree table file", self.tree_table)):
            if not p.exists():
                raise SkipTumour(f"missing {label}: {p}")
            if p.stat().st_size == 0:
                raise SkipTumour(f"empty {label}: {p}")


# ---- tree ------------------------------------------------------------------

def default_tree_edges(all_trees: Path) -> list[tuple[int, int]]:
    """Edges of the first tree in a CONIPHER allTrees.txt."""
    lines = Path(all_trees).read_text().splitlines()
    starts = [i for i, ln in enumerate(lines) if ln.startswith("# tree ")]
    if starts:
        end = starts[1] if len(starts) > 1 else len(lines)
        lines = lines[starts[0] + 1:end]
    edges = [ln.split("\t") for ln in lines if ln.strip() and not ln.startswith("#")]
    if not edges:
        raise SkipTumour(f"no edges found in {all_trees}")
    if len(edges) == 1:
        raise SkipTumour(f"single cluster tree found in {all_trees}")
    return [(int(p), int(c)) for p, c, *_ in edges]


def clonal_cluster(edges: list[tuple[int, int]]) -> int:
    """Root of the tree: a parent that is never a child. Smallest if several."""
    parents = {p for p, _ in edges}
    children = {c for _, c in edges}
    roots = sorted(parents - children)
    if not roots:
        raise SkipTumour("could not identify the root")
    return roots[0]


def tree_clusters(edges: list[tuple[int, int]]) -> set[int]:
    return {n for e in edges for n in e}


# ---- CN: num_gds and doubled-genome fractions -------------------------------

def _level_confidence(a: np.ndarray, b: np.ndarray, w: np.ndarray, target: float) -> float:
    """Length-weighted closeness to (target | target): 0 at (target/2 | target/2), 1 at (target | target).

    CN is capped at `target` so that very high segments are not over-penalised.
    """
    ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(w) & (w > 0)
    if not ok.any():
        return float("nan")
    a = np.clip(a[ok], 0, target)
    b = np.clip(b[ok], 0, target)
    w = w[ok]
    deficit = np.sum(((target - a) + (target - b)) * w) / np.sum(w)
    return float(np.clip(1 - deficit / target, 0, 1))


def compute_num_gds(cn: pl.DataFrame) -> pl.DataFrame:
    """Per SAMPLE: frac_maj_ge2/3/4, num_gds (0-3) and conf_gd1/2/3."""
    maj = cn.select(pl.max_horizontal("COPY_NUMBER_A", "COPY_NUMBER_B").cast(pl.Float64))
    cn = cn.with_columns(
        (pl.col("ENDPOS") - pl.col("STARTPOS") + 1).alias("seg_len"),
        # round half to even, like R's round(); missing CN counts in the denominator only
        pl.Series("MajCN_round", np.round(maj.to_series().to_numpy())).fill_nan(None),  # polars sorts NaN above all numbers
    )
    frac = {k: (pl.col("seg_len").filter(pl.col("MajCN_round") >= k).sum() / pl.col("seg_len").sum())
            .alias(f"frac_maj_ge{k}") for k in (2, 3, 4)}
    out = (cn.group_by("SAMPLE", maintain_order=True)
             .agg(*frac.values(), pl.col("COPY_NUMBER_A"), pl.col("COPY_NUMBER_B"), pl.col("seg_len")))
    out = out.with_columns(
        pl.when(pl.col("frac_maj_ge4") >= 0.5).then(3)
          .when(pl.col("frac_maj_ge3") >= 0.5).then(2)
          .when(pl.col("frac_maj_ge2") >= 0.5).then(1)
          .otherwise(0).cast(pl.Int32).alias("num_gds"))
    conf = {f"conf_gd{i}": [] for i in (1, 2, 3)}
    for a, b, w in out.select("COPY_NUMBER_A", "COPY_NUMBER_B", "seg_len").iter_rows():
        a, b, w = (np.asarray(x, dtype=float) for x in (a, b, w))
        for i, target in ((1, 2.0), (2, 4.0), (3, 8.0)):
            conf[f"conf_gd{i}"].append(_level_confidence(a, b, w, target))
    return (out.drop("COPY_NUMBER_A", "COPY_NUMBER_B", "seg_len")
               .with_columns(*[pl.Series(k, v, dtype=pl.Float64).fill_nan(None) for k, v in conf.items()])
               .select("SAMPLE", "frac_maj_ge2", "frac_maj_ge3", "frac_maj_ge4", "num_gds",
                       "conf_gd1", "conf_gd2", "conf_gd3"))


# ---- per tumour --------------------------------------------------------------

def _read(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", null_values=["NA"], infer_schema_length=None,
                       schema_overrides={"SAMPLE": pl.Utf8, "mutation_id": pl.Utf8, "CASE_ID": pl.Utf8,
                                         "CHR": pl.Utf8})


def _require(df: pl.DataFrame, cols, label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SkipTumour(f"missing column(s) in {label}: {', '.join(missing)}")


def prepare_tumour(files: TumourFiles, remove_nontree_clusters: bool = True) -> pl.DataFrame:
    """PGDD input rows for one tumour. Raises SkipTumour on unusable input."""
    files.check()
    muts = _read(files.mutations)
    _require(muts, MUT_COLS, "mutations")
    tree_table = _read(files.tree_table)
    _require(tree_table, ("mutation_id", "SAMPLE", "treeCLUSTER"), "tree table")

    # keep mutations that made it onto the tree; their cluster is the tree cluster
    key = pl.concat_str("mutation_id", "SAMPLE", separator=":").alias("mutation_id_sample")
    muts = (muts.with_columns(key)
                .join(tree_table.select(key, "treeCLUSTER"), on="mutation_id_sample", how="inner")
                .with_columns(pl.col("treeCLUSTER").alias("CLUSTER")))

    edges = default_tree_edges(files.all_trees)
    root = clonal_cluster(edges)
    on_tree = tree_clusters(edges)

    cn = _read(files.cn)
    _require(cn, CN_COLS, "CN segments")
    gds = compute_num_gds(cn)

    mut_samples = {s for s in muts["SAMPLE"].drop_nulls() if s.strip()}
    cn_samples = {s for s in gds["SAMPLE"].drop_nulls() if s.strip()}
    if mut_samples != cn_samples:
        raise SkipTumour(
            "SAMPLE ids differ between mutations and CN: "
            f"only in mutations {sorted(mut_samples - cn_samples)[:10]}, "
            f"only in CN {sorted(cn_samples - mut_samples)[:10]}")

    muts = muts.join(gds, on="SAMPLE", how="left")
    if remove_nontree_clusters:
        muts = muts.filter(pl.col("CLUSTER").is_in(list(on_tree)))

    # row order as in R: merge() sorts by SAMPLE, then by mutation_id_sample
    muts = muts.sort("SAMPLE", "mutation_id_sample")
    ploidy = pl.col("PLOIDY").cast(pl.Float64) if "PLOIDY" in muts.columns else pl.lit(None, pl.Float64)
    out = muts.select(
        pl.col("CASE_ID").alias("tumour_id"),
        (pl.lit("chr") + pl.col("CHR").cast(pl.Utf8)).alias("chromosome"),
        pl.col("POS").alias("position"),
        pl.col("REF").alias("ref"),
        pl.col("ALT").alias("alt"),
        pl.col("CLUSTER").alias("cluster_id"),
        (pl.col("CLUSTER") == root).alias("is_clonal_cluster"),
        pl.col("SAMPLE").alias("sample_id"),
        pl.col("MUT_COPY").alias("mut_cpn"),
        pl.max_horizontal("COPY_NUMBER_A", "COPY_NUMBER_B").alias("MajCN"),
        "num_gds", "frac_maj_ge2", "frac_maj_ge3", "frac_maj_ge4",
        "conf_gd1", "conf_gd2", "conf_gd3",
        ploidy.alias("ploidy"),
    )
    if out.is_empty():
        raise SkipTumour("no PGDD input rows after filtering")
    return out


@dataclass
class PrepareResult:
    table: pl.DataFrame
    skipped: dict[str, str] = field(default_factory=dict)


def write_tsv(table: pl.DataFrame, path: Path) -> None:
    """Write for R: booleans as TRUE/FALSE, missing as NA."""
    table.with_columns(
        pl.when(pl.col(c)).then(pl.lit("TRUE")).when(~pl.col(c)).then(pl.lit("FALSE")).alias(c)
        for c, t in table.schema.items() if t == pl.Boolean
    ).write_csv(path, separator="\t", null_value="NA")


def canonical_tumour_id(col: str = "tumour_id") -> pl.Expr:
    """CONIPHER CASE_ID `LTX0001_1` -> `LTX0001-Tumour1` (R: gsub('_', '-Tumour', ...))."""
    return pl.col(col).str.replace_all("_", "-Tumour", literal=True)


def prepare_all(conipher_dir: Path, msai_dir: Path, tumours: list[str] | None = None,
                remove_nontree_clusters: bool = True, log=print) -> PrepareResult:
    """PGDD input for every tumour directory under `conipher_dir` (or the `tumours` subset)."""
    dirs = sorted(p.name for p in Path(conipher_dir).iterdir() if p.is_dir())
    if not dirs:
        raise FileNotFoundError(f"no tumour directories in {conipher_dir}")
    if tumours is not None:
        missing = sorted(set(tumours) - set(dirs))
        if missing:
            log(f"WARNING: selected tumours not found: {', '.join(missing)}")
        dirs = [d for d in dirs if d in set(tumours)]

    tables, skipped = [], {}
    for t in dirs:
        try:
            df = prepare_tumour(TumourFiles.locate(conipher_dir, msai_dir, t), remove_nontree_clusters)
        except SkipTumour as e:
            skipped[t] = str(e)
            log(f"[{t}] skipped: {e}")
            continue
        log(f"[{t}] samples={df['sample_id'].n_unique()} mutations={df.height} clusters={df['cluster_id'].n_unique()}")
        tables.append(df)

    if not tables:
        raise RuntimeError("no PGDD input rows produced; see skipped tumours")
    table = pl.concat(tables, how="vertical_relaxed").with_columns(canonical_tumour_id())
    return PrepareResult(table.select(OUTPUT_COLS), skipped)
