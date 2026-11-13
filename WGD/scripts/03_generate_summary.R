suppressPackageStartupMessages({
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 03_pgdd_cohort_summary.R <pgdd_input_file> <pgdd_regions_file> <output_dir>")
}
pgdd_input_file   = args[1]
pgdd_regions_file = args[2]
output_dir        = args[3]

message("Reading PGDD input: ", pgdd_input_file)
pgdd_input = fread(pgdd_input_file)

message("Reading PGDD per-region output: ", pgdd_regions_file)
pgdd_regions = fread(pgdd_regions_file)

col_or_na = function(dt, col, type = NA_real_) {
  if (col %in% names(dt)) dt[[col]] else rep(type, nrow(dt))
}

# Per-sample CN summary — these values are constant across mutations for a given sample
sample_cn = unique(pgdd_input[, .(
  sample_name_hash  = sample_id,
  Ploidy            = col_or_na(pgdd_input, "ploidy"),
  genome_frac_nMaj2 = col_or_na(pgdd_input, "frac_maj_ge2"),
  genome_frac_nMaj3 = col_or_na(pgdd_input, "frac_maj_ge3"),
  genome_frac_nMaj4 = col_or_na(pgdd_input, "frac_maj_ge4"),
  num_gds           = col_or_na(pgdd_input, "num_gds")
)])

# add tumour_id to sample_cn:
sample_to_tumour_map = unique(pgdd_regions[, .(sample_id, tumour_id)])
sample_cn = merge(sample_cn, sample_to_tumour_map, by.x = "sample_name_hash", by.y = "sample_id", all.x = TRUE)

# Per-sample GD predictions from pgdd algorithm
sample_gd = pgdd_regions[, .(
  sample_name_hash      = sample_id,
  `GD predicted`        = num_gds,
  `Second GD predicted` = pmax(0L, num_gds - 1L)
)]

summary_dt = merge(sample_cn, sample_gd, by = "sample_name_hash", all = TRUE)

setcolorder(summary_dt, c(
  "sample_name_hash",
  "GD predicted",
  "Second GD predicted",
  "Ploidy",
  "genome_frac_nMaj2",
  "genome_frac_nMaj3",
  "genome_frac_nMaj4",
  "num_gds"
))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file = file.path(output_dir, "pgdd_cohort_gd_summary.tsv")
fwrite(summary_dt, output_file, sep = "\t")

message("Written: ", output_file)
message("Samples: ", nrow(summary_dt))
