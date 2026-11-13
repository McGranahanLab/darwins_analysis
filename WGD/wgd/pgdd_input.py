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
* Single-clone tumours (CONIPHER tree `1 \t 1`) are not skipped. They are not
  sent to PGDD either: their GDs are all clonal on the one clone and are counted
  from `num_gds` (Piotr, 2026-11-09), see `single_clone_events`. The R script
  skipped every tree with exactly one edge, which also dropped real two-clone
  trees (`1 \t 2`); those now go to PGDD.
* A pre-flight pass checks all tumours first and stops the run, listing every
  problem, if any tumour has a missing/empty file, missing columns, an unreadable
  tree, or SAMPLE ids that differ between SCoutput, treeTable and the CN segments
  (Piotr, 2026-11-09). The R script skipped such tumours.
"""

from __future__ import annotations

from dataclasses import dataclass
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


SINGLE_CLONE_SCHEMA = {
    "tumour_id": pl.Utf8, "sample_id": pl.Utf8, "cluster_id": pl.Int64, "clone": pl.Utf8, "num_gds": pl.Int32,
    "frac_maj_ge2": pl.Float64, "frac_maj_ge3": pl.Float64, "frac_maj_ge4": pl.Float64,
    "conf_gd1": pl.Float64, "conf_gd2": pl.Float64, "conf_gd3": pl.Float64, "ploidy": pl.Float64,
}


class TumourInputError(Exception):
    """Input problem in one tumour. Collected by the pre-flight, then fatal (InputError)."""


class InputError(RuntimeError):
    """One or more tumours have unusable inputs: missing/empty files, missing columns,
    an unreadable tree, or SAMPLE ids that differ between inputs. The run stops."""


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
                raise TumourInputError(f"missing {label}: {p}")
            if p.stat().st_size == 0:
                raise TumourInputError(f"empty {label}: {p}")


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
        raise TumourInputError(f"no edges found in {all_trees}")
    return [(int(p), int(c)) for p, c, *_ in edges]


def is_single_clone(edges: list[tuple[int, int]]) -> bool:
    """CONIPHER writes a single-cluster tree as one self-loop edge, `1 \\t 1`."""
    return len(tree_clusters(edges)) == 1


def clonal_cluster(edges: list[tuple[int, int]]) -> int:
    """Root of the tree: a parent that is never a child. Smallest if several."""
    if is_single_clone(edges):
        return edges[0][0]
    parents = {p for p, _ in edges}
    children = {c for _, c in edges}
    roots = sorted(parents - children)
    if not roots:
        raise TumourInputError("could not identify the root")
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


# ---- pre-flight ----------------------------------------------------------------

def _samples(path: Path) -> set[str]:
    col = pl.read_csv(path, separator="\t", columns=["SAMPLE"], schema_overrides={"SAMPLE": pl.Utf8})["SAMPLE"]
    return {s for s in col.drop_nulls() if s.strip()}


def _check_header(path: Path, cols, label: str) -> None:
    header = path.open().readline().rstrip("\n").split("\t")
    missing = [c for c in cols if c not in header]
    if missing:
        raise TumourInputError(f"missing column(s) in {label}: {', '.join(missing)}")


def _preflight_tumour(files: TumourFiles) -> None:
    files.check()
    _check_header(files.mutations, MUT_COLS, "mutations")
    _check_header(files.tree_table, ("mutation_id", "SAMPLE", "treeCLUSTER"), "tree table")
    _check_header(files.cn, CN_COLS, "CN segments")
    clonal_cluster(default_tree_edges(files.all_trees))
    sets = {"mutations": _samples(files.mutations), "treeTable": _samples(files.tree_table),
            "CN": _samples(files.cn)}
    union = set.union(*sets.values())
    missing = {k: sorted(union - v) for k, v in sets.items() if union - v}
    if missing:
        raise TumourInputError("SAMPLE ids differ; missing from " +
                               "; ".join(f"{k}: {v}" for k, v in missing.items()))


def preflight(conipher_dir: Path, msai_dir: Path, tumours: list[str]) -> None:
    """Check every tumour before any processing; raise InputError listing all problems.

    Checked: files present and non-empty, required columns, a readable default tree,
    and identical SAMPLE sets in SCoutput, treeTable and the CN segments.
    """
    problems = []
    for t in tumours:
        try:
            _preflight_tumour(TumourFiles.locate(conipher_dir, msai_dir, t))
        except TumourInputError as e:
            problems.append(f"  {t}: {e}")
    if problems:
        raise InputError(f"input problems in {len(problems)} tumour(s):\n" + "\n".join(problems))


# ---- per tumour --------------------------------------------------------------

def _read(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", null_values=["NA"], infer_schema_length=None,
                       schema_overrides={"SAMPLE": pl.Utf8, "mutation_id": pl.Utf8, "CASE_ID": pl.Utf8,
                                         "CHR": pl.Utf8})


def _require(df: pl.DataFrame, cols, label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise TumourInputError(f"missing column(s) in {label}: {', '.join(missing)}")


def prepare_tumour(files: TumourFiles, remove_nontree_clusters: bool = True
                   ) -> tuple[pl.DataFrame | None, pl.DataFrame | None]:
    """(PGDD input rows, None) for a tumour with a tree, or (None, per-sample rows) for a single-clone tumour.

    Single-clone tumours do not go to PGDD: all their GDs are clonal, on the one
    clone, and are counted from the CN-based `num_gds` (Piotr, 2026-11-09).
    Raises TumourInputError on unusable input.
    """
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
    if mut_samples != cn_samples:  # pre-flight should have caught this
        raise TumourInputError(
            f"{files.tumour_id}: SAMPLE ids differ between tree mutations and CN: "
            f"only in mutations {sorted(mut_samples - cn_samples)[:10]}, "
            f"only in CN {sorted(cn_samples - mut_samples)[:10]}")

    ploidy = pl.col("PLOIDY").cast(pl.Float64) if "PLOIDY" in muts.columns else pl.lit(None, pl.Float64)

    if is_single_clone(edges):
        per_sample = muts.group_by("SAMPLE").agg(ploidy.first().alias("ploidy"),
                                                 pl.col("CASE_ID").first().alias("tumour_id"))
        out = (gds.join(per_sample, on="SAMPLE", how="left")
                  .select("tumour_id", pl.col("SAMPLE").alias("sample_id"),
                          pl.lit(root).alias("cluster_id"), pl.lit(f"clone{root}").alias("clone"),
                          "num_gds", "frac_maj_ge2", "frac_maj_ge3", "frac_maj_ge4",
                          "conf_gd1", "conf_gd2", "conf_gd3", "ploidy")
                  .sort("sample_id"))
        return None, out

    muts = muts.join(gds, on="SAMPLE", how="left")
    if remove_nontree_clusters:
        muts = muts.filter(pl.col("CLUSTER").is_in(list(on_tree)))

    # row order as in R: merge() sorts by SAMPLE, then by mutation_id_sample
    muts = muts.sort("SAMPLE", "mutation_id_sample")
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
        raise TumourInputError("no PGDD input rows after filtering")
    return out, None


@dataclass
class PrepareResult:
    table: pl.DataFrame                 # PGDD input (tumours with a tree)
    single_clone: pl.DataFrame          # per-sample rows of single-clone tumours (not run through PGDD)


# Quality codes (`quality_code`, ";"-joined) explain a quality below High.
QUALITY_CODES = {
    "single_clone_num_gds_discordant": "single-clone tumour whose samples disagree on num_gds; the mode is used",
    "num_gds_mode_tie": "no unique mode of num_gds across samples; the lowest tied value is used",
}


def single_clone_events(single_clone: pl.DataFrame) -> pl.DataFrame:
    """One row per single-clone tumour: `n_clonal_gds` clonal GDs on its only clone (the MRCA).

    The clone is the same in every sample, so samples should agree on `num_gds`.
    If they don't (Piotr, 2026-11-09): `n_clonal_gds` = the mode (lowest value on a
    tie), quality = Low, and `quality_code` says why.
    """
    rows = []
    for (tumour, clone), g in single_clone.group_by("tumour_id", "clone", maintain_order=True):
        counts = g["num_gds"].value_counts(sort=True)
        top = counts["count"].max()
        modes = sorted(counts.filter(pl.col("count") == top)["num_gds"].to_list())
        codes = []
        if g["num_gds"].n_unique() > 1:
            codes.append("single_clone_num_gds_discordant")
        if len(modes) > 1:
            codes.append("num_gds_mode_tie")
        rows.append({
            "tumour_id": tumour, "clone": clone, "n_clonal_gds": modes[0],
            "n_samples": g.height,
            "num_gds_per_sample": ",".join(f"{s}:{n}" for s, n in g.sort("sample_id").select("sample_id", "num_gds").iter_rows()),
            "quality": "Low" if codes else "High",
            "quality_code": ";".join(codes) or None,
        })
    schema = {"tumour_id": pl.Utf8, "clone": pl.Utf8, "n_clonal_gds": pl.Int32, "n_samples": pl.Int64,
              "num_gds_per_sample": pl.Utf8, "quality": pl.Utf8, "quality_code": pl.Utf8}
    return pl.DataFrame(rows, schema=schema).sort("tumour_id")


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
    """PGDD input for every tumour directory under `conipher_dir` (or the `tumours` subset).

    Runs `preflight` first: any input problem stops the run before processing.
    """
    dirs = sorted(p.name for p in Path(conipher_dir).iterdir() if p.is_dir())
    if not dirs:
        raise FileNotFoundError(f"no tumour directories in {conipher_dir}")
    if tumours is not None:
        missing = sorted(set(tumours) - set(dirs))
        if missing:
            raise InputError(f"selected tumours not found in {conipher_dir}: {', '.join(missing)}")
        dirs = [d for d in dirs if d in set(tumours)]

    preflight(conipher_dir, msai_dir, dirs)

    tables, single = [], []
    for t in dirs:
        try:
            df, sc = prepare_tumour(TumourFiles.locate(conipher_dir, msai_dir, t), remove_nontree_clusters)
        except TumourInputError as e:
            raise InputError(f"{t}: {e}") from e
        if sc is not None:
            log(f"[{t}] single clone: samples={sc.height} num_gds={sorted(set(sc['num_gds']))} (not sent to PGDD)")
            single.append(sc)
        else:
            log(f"[{t}] samples={df['sample_id'].n_unique()} mutations={df.height} clusters={df['cluster_id'].n_unique()}")
            tables.append(df)

    table = (pl.concat(tables, how="vertical_relaxed").with_columns(canonical_tumour_id()).select(OUTPUT_COLS)
             if tables else pl.DataFrame(schema={c: pl.Utf8 for c in OUTPUT_COLS}))
    single_clone = (pl.concat(single, how="vertical_relaxed").with_columns(canonical_tumour_id())
                    if single else pl.DataFrame(schema=SINGLE_CLONE_SCHEMA))
    return PrepareResult(table, single_clone)
