"""Read ParallelGDDetect (PGDD) outputs and reshape them into long tables.

PGDD (amf71/ParallelGDDetect, bootstrap version) writes one directory per run:

  pgdd_GDs_per_tumour.tsv     one row per tumour
  pgdd_GDs_per_region.tsv     one row per region (sample)
  pgdd_GD_events.tsv          one row per GD event
  pgdd_mut_counts.tsv         per tumour x sample x cluster mutation counts
  pgdd_bootstrap_discovery_summary.tsv / _settings.tsv

Only the stability-filtered (`pgdd_*`) tables are read. `unstable_calls_*`
files are for manual review of rejected calls and are never ingested.

Multi-valued cells are comma-joined strings: `GD_events.clusters` ("3,6" =
event supported by clusters 3 and 6), `GD_events.samples`,
`GDs_per_region.gd_clusters` and `GDs_per_region.gd_events`. A single value
("5") looks numeric, so these columns are always read as strings and
exploded here, never parsed as numbers.

PGDD clusters are pre-tree mutation clusters: a cluster id may not exist on
the final CONIPHER tree. Pruning against the tree is a separate step.
"""

from __future__ import annotations

from pathlib import Path, PurePosixPath

import polars as pl

# Real location of the WES run (old run; to be regenerated with new hg19 ALPACA, then hg38).
PGDD_NEMO_DIR = PurePosixPath(
    "/nemo/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842/WGD/release_tx842/pgdd"
)

FILES = {
    "tumours": "pgdd_GDs_per_tumour.tsv",
    "regions": "pgdd_GDs_per_region.tsv",
    "events": "pgdd_GD_events.tsv",
    "mut_counts": "pgdd_mut_counts.tsv",
    "bootstrap": "pgdd_bootstrap_discovery_summary.tsv",
    "settings": "pgdd_bootstrap_discovery_settings.tsv",
}

# Columns that must stay strings (comma lists or ids that look numeric).
STRING_COLS = {
    "tumour_id", "sample_id", "cluster_id", "clusters", "samples",
    "gd_clusters", "gd_events", "GD_event", "GD_event_id", "GD_statuses",
    "First_GD", "Second_GD",
}


def pgdd_dir(nemo_root: Path) -> Path:
    """Re-root the Nemo PGDD directory under `nemo_root` (`/nemo` on the HPC, `nemo_mock` locally)."""
    return Path(nemo_root) / PGDD_NEMO_DIR.relative_to("/nemo")


def _read_tsv(path: Path) -> pl.DataFrame:
    header = path.open().readline().rstrip("\n").split("\t")
    overrides = {c: pl.Utf8 for c in header if c in STRING_COLS}
    df = pl.read_csv(path, separator="\t", schema_overrides=overrides,
                     null_values=["NA"], infer_schema_length=None)
    # fwrite writes empty strings as "" -> normalise empty to null
    return df.with_columns(
        [pl.when(pl.col(c).str.len_chars() == 0).then(None).otherwise(pl.col(c)).alias(c)
         for c in df.columns if df.schema[c] == pl.Utf8]
    )


def read_pgdd(directory: Path) -> dict[str, pl.DataFrame]:
    """Read the stability-filtered PGDD tables from one run directory."""
    directory = Path(directory)
    missing = [f for f in FILES.values() if not (directory / f).exists()]
    if missing:
        raise FileNotFoundError(f"PGDD files missing in {directory}: {missing}")
    return {key: _read_tsv(directory / fname) for key, fname in FILES.items()}


def _explode_list(df: pl.DataFrame, col: str, new: str) -> pl.DataFrame:
    return (df.with_columns(pl.col(col).str.split(",").alias(new))
              .explode(new)
              .with_columns(pl.col(new).str.strip_chars()))


def event_clusters(events: pl.DataFrame) -> pl.DataFrame:
    """One row per (tumour, GD event, supporting cluster).

    Events without mutation support (clonal or ploidy-only) keep one row
    with `cluster_id` null, so every event is present exactly once when
    grouped by GD_event_id.
    """
    out = _explode_list(events, "clusters", "cluster_id")
    return out.select(
        "tumour_id", "GD_event", "GD_event_id", "is_clonal",
        "is_subclonal_mutation_supported", "cluster_id",
        pl.col("clusters").alias("clusters_raw"),
    ).sort("tumour_id", "GD_event", "cluster_id", nulls_last=True)


def event_regions(events: pl.DataFrame) -> pl.DataFrame:
    """One row per (tumour, GD event, region carrying the event)."""
    out = _explode_list(events, "samples", "sample_id")
    return out.select(
        "tumour_id", "GD_event", "GD_event_id", "is_clonal", "sample_id",
        "num_regions_gd_event", "total_regions",
    ).sort("tumour_id", "GD_event", "sample_id")


def region_clusters(regions: pl.DataFrame) -> pl.DataFrame:
    """One row per (region, doubled subclonal cluster) from `GDs_per_region.gd_clusters`.

    Regions with no doubled cluster are dropped (they carry no cluster evidence).
    """
    out = _explode_list(regions.filter(pl.col("gd_clusters").is_not_null()),
                        "gd_clusters", "cluster_id")
    return out.select("tumour_id", "sample_id", "num_gds", "cluster_id").sort(
        "tumour_id", "sample_id", "cluster_id")


def summarise(pgdd: dict[str, pl.DataFrame]) -> dict[str, int]:
    ev = pgdd["events"]
    ec = event_clusters(ev)
    return {
        "tumours": pgdd["tumours"].height,
        "regions": pgdd["regions"].height,
        "events": ev.height,
        "clonal_events": ev.filter(pl.col("is_clonal")).height,
        "subclonal_events": ev.filter(~pl.col("is_clonal")).height,
        "subclonal_mut_supported": ev.filter(~pl.col("is_clonal") & pl.col("is_subclonal_mutation_supported")).height,
        "multi_cluster_events": ev.filter(pl.col("clusters").str.contains(",", literal=True)).height,
        "event_cluster_rows": ec.filter(pl.col("cluster_id").is_not_null()).height,
    }
