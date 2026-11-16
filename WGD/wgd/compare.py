"""Step 06: compare the new per-clone table with the previous release table.

The previous table went through manual review in wgd_viewer, so some
differences are expected (status, and placements decided by hand). The point
is to see whether the automatic logic was ported faithfully: tumours whose
calls the notebook made automatically (no_wgd, clonal_wgd, ploidy_only with
all ratios > 1.5) should agree exactly.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

KEY = ("tumour_id", "clone")
FIELDS = ("parent", "number_of_gds_at_clone", "total_gds_per_tumour", "quality", "status", "class", "sub_class")


def read_table(path: Path) -> pl.DataFrame:
    """Read the old table (TSV or CSV; a leading unnamed index column is dropped)."""
    path = Path(path)
    sep = "," if path.suffix == ".csv" else "\t"
    df = pl.read_csv(path, separator=sep, infer_schema_length=None, null_values=["NA", "NaN", "nan", ""])
    first = df.columns[0]
    if first in ("", "index", "Unnamed: 0") or first.startswith("column_"):
        df = df.drop(first)
    return df


def _same(a: pl.Expr, b: pl.Expr) -> pl.Expr:
    """Equal, treating null == null as equal; numbers compared as strings of ints."""
    return (a.cast(pl.Utf8) == b.cast(pl.Utf8)).fill_null(a.is_null() & b.is_null())


def compare(old: pl.DataFrame, new: pl.DataFrame) -> dict[str, pl.DataFrame | dict]:
    fields = [f for f in FIELDS if f in old.columns and f in new.columns]
    old_t, new_t = set(old["tumour_id"]), set(new["tumour_id"])
    both = old_t & new_t

    o = old.filter(pl.col("tumour_id").is_in(list(both))).select(*KEY, *fields)
    n = new.filter(pl.col("tumour_id").is_in(list(both))).select(*KEY, *fields)
    j = o.join(n, on=list(KEY), how="full", suffix="_new", coalesce=True)
    j = j.with_columns(
        pl.col(fields[0]).is_null().and_(pl.col(f"{fields[0]}_new").is_not_null()).alias("_only_new"),
        pl.col(f"{fields[0]}_new").is_null().and_(pl.col(fields[0]).is_not_null()).alias("_only_old"),
        *[_same(pl.col(f), pl.col(f"{f}_new")).alias(f"{f}_match") for f in fields],
    )
    clone_diff = j.filter(~pl.all_horizontal([pl.col(f"{f}_match") for f in fields if f != "status"])
                          | pl.col("_only_new") | pl.col("_only_old"))

    per_tumour = (j.group_by("tumour_id").agg(
        pl.col("class").first().alias("class_old"), pl.col("class_new").first(),
        pl.col("status").first().alias("status_old"), pl.col("status_new").first(),
        pl.col("total_gds_per_tumour").first().alias("total_gds_old"),
        pl.col("total_gds_per_tumour_new").first().alias("total_gds_new"),
        (pl.col("_only_new") | pl.col("_only_old")).any().alias("clone_sets_differ"),
        *[pl.col(f"{f}_match").all() for f in fields],
    ).with_columns(
        pl.all_horizontal("number_of_gds_at_clone_match", ~pl.col("clone_sets_differ")).alias("placement_match"),
    ).sort("tumour_id"))

    def rate(col: str) -> dict:
        return {"match": int(per_tumour[col].sum()), "of": per_tumour.height}

    summary = {
        "tumours_old": len(old_t), "tumours_new": len(new_t), "tumours_both": len(both),
        "tumours_only_old": sorted(old_t - new_t), "tumours_only_new": sorted(new_t - old_t),
        "placement_match": rate("placement_match"),
        **{f"{f}_match": rate(f"{f}_match") for f in fields},
        # status_new = resolved: the notebook decided automatically, so these must match.
        # unresolved: the old flow sent them to manual review in wgd_viewer; differences are expected.
        "placement_match_by_new_status": {
            str(st): {"match": int(g["placement_match"].sum()), "of": g.height}
            for (st,), g in per_tumour.group_by("status_new")},
        "placement_match_by_old_class": {
            str(c): {"match": int(g["placement_match"].sum()), "of": g.height}
            for (c,), g in per_tumour.group_by("class_old")},
    }
    return {"summary": summary, "per_tumour": per_tumour, "clone_diff": clone_diff.drop("_only_new", "_only_old")}
