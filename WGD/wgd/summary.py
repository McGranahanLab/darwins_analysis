"""Per-sample PGDD cohort summary (`pgdd_cohort_gd_summary.tsv`).

Python port of `scripts/03_generate_summary.R` (Piotr, 2026-11-09). Joins the
per-sample CN summary carried in the PGDD input (step 00) with PGDD's
per-region calls (`pgdd_GDs_per_region.tsv`):

  sample_name_hash      sample_id
  GD predicted          PGDD num_gds for the region
  Second GD predicted   max(0, GD predicted - 1)
  Ploidy                PLOIDY from CONIPHER
  genome_frac_nMaj2/3/4 fraction of the genome (bp) with round(MajCN) >= 2/3/4;
                        nMaj2 is the doubled-genome fraction
  num_gds               ploidy-based GD count from the CN segments (step 00)
  tumour_id

Differences from the R script: `tumour_id` falls back to the PGDD input for
samples missing from the regions table, and rows are sorted by tumour and sample.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

SAMPLE_CN = {  # summary column <- PGDD input column
    "Ploidy": "ploidy",
    "genome_frac_nMaj2": "frac_maj_ge2",
    "genome_frac_nMaj3": "frac_maj_ge3",
    "genome_frac_nMaj4": "frac_maj_ge4",
    "num_gds": "num_gds",
}
COLUMNS = ("sample_name_hash", "GD predicted", "Second GD predicted", "Ploidy",
           "genome_frac_nMaj2", "genome_frac_nMaj3", "genome_frac_nMaj4", "num_gds", "tumour_id")


def read_tsv(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", null_values=["NA", ""], infer_schema_length=None,
                       schema_overrides={"sample_id": pl.Utf8, "tumour_id": pl.Utf8,
                                         "gd_clusters": pl.Utf8, "gd_events": pl.Utf8})


def cohort_summary(pgdd_input: pl.DataFrame, pgdd_regions: pl.DataFrame) -> pl.DataFrame:
    """One row per sample. Raises if a sample has more than one CN summary in the input."""
    sample_cn = pgdd_input.select(
        pl.col("sample_id").alias("sample_name_hash"),
        *[(pl.col(src) if src in pgdd_input.columns else pl.lit(None, pl.Float64)).alias(dst)
          for dst, src in SAMPLE_CN.items()],
    ).unique(maintain_order=True)
    dup = sample_cn.filter(pl.col("sample_name_hash").is_duplicated())
    if not dup.is_empty():
        raise ValueError(f"per-sample CN values are not constant for: {dup['sample_name_hash'].unique().to_list()[:10]}")

    regions = pgdd_regions.select(
        pl.col("sample_id").alias("sample_name_hash"),
        "tumour_id",
        pl.col("num_gds").cast(pl.Int64).alias("GD predicted"),
    ).with_columns((pl.col("GD predicted") - 1).clip(lower_bound=0).alias("Second GD predicted"))

    input_tumour = pgdd_input.select(pl.col("sample_id").alias("sample_name_hash"),
                                     pl.col("tumour_id").alias("_tumour_input")).unique()
    return (sample_cn.join(regions, on="sample_name_hash", how="full", coalesce=True)
                     .join(input_tumour, on="sample_name_hash", how="left")
                     .with_columns(pl.coalesce("tumour_id", "_tumour_input").alias("tumour_id"))
                     .select(COLUMNS)
                     .sort("tumour_id", "sample_name_hash", nulls_last=True))
