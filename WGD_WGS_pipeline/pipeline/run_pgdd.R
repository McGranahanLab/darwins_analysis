#!/usr/bin/env Rscript
# Run ParallelGDDetect (detect_par_gd) on the step 00 input with the settings from config.toml.
#
#   Rscript run_pgdd.R <pgdd_src> <pgdd_input_all_tumours.tsv> <pgdd_settings.json> <out_dir>
#
# Every setting comes from [pgdd.settings] (written to JSON by wgs.py), so nothing relies on
# PGDD's defaults. Output file names are the ones the lab's run wrapper writes and step 01 reads;
# unstable_calls_* are for manual review only.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: run_pgdd.R <pgdd_src> <pgdd_input.tsv> <settings.json> <out_dir>")
pgdd_src <- path.expand(args[1])
source(file.path(pgdd_src, "R", "parrallel_gd_caller.R"))

input <- fread(args[2], sep = "\t")
settings <- fromJSON(args[3], simplifyVector = TRUE)
out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("PGDD: ", uniqueN(input$tumour_id), " tumours, ", nrow(input), " rows; settings: ",
        paste(names(settings), unlist(settings), sep = "=", collapse = ", "))
# PGDD returns nothing when no sample has a GD (`all(num_gds == 0)`, also for an empty input: e.g. only
# single-clone tumours). Write what it writes for no-GD tumours in a normal run instead: "No GD" rows
# per tumour and sample, other tables empty with PGDD's columns (2026-11-18).
HEADERS <- list(
  GDs_per_tumour = c("tumour_id", "First_GD", "Second_GD", "num_first_gd", "num_second_gd", "First_GD_homogen",
                     "Second_GD_homogen", "GD_status_homogen", "GD_statuses", "frac_0_gd_regions",
                     "frac_1_gd_regions", "frac_2_gd_regions", "num_clonal_gds", "num_subclonal_gds"),
  GDs_per_region = c("sample_id", "tumour_id", "gd_clusters", "num_gds", "gd_events", "First_GD", "Second_GD"),
  GDs_events = c("GD_event", "GD_event_id", "tumour_id", "is_clonal", "num_regions_gd_event", "total_regions",
                 "is_subclonal_mutation_supported", "clusters", "samples"),
  mut_counts = c("tumour_id", "sample_id", "cluster_id", "num_muts", "num_cn2", "num_cn2_all", "num_muts_present",
                 "perc_cn2", "num_cn2_all_check", "perc_cn2_all_check", "num_total_present", "gd_p_value",
                 "is_subcl_gd_discovery", "is_subcl_gd", "is_subcl_gd_any_region", "discover_mut_cpn_2_threshold",
                 "check_mut_cpn_2_threshold", "discover_num_muts_value", "discover_num_muts_threshold",
                 "discover_num_muts_pass", "discover_frac_2_cpn_muts_value", "discover_frac_2_cpn_muts_threshold",
                 "discover_frac_2_cpn_muts_pass", "discover_num_2_cpn_muts_value", "discover_num_2_cpn_muts_threshold",
                 "discover_num_2_cpn_muts_pass", "check_frac_2_cpn_muts_value", "check_frac_2_cpn_muts_threshold",
                 "check_frac_2_cpn_muts_pass", "check_num_2_cpn_muts_value", "check_num_2_cpn_muts_threshold",
                 "check_num_2_cpn_muts_pass", "discover_pass_all", "check_pass_all"),
  bootstrap_discovery = c("tumour_id", "sample_id", "cluster_id", "num_muts", "perc_cn2", "num_cn2", "k", "B",
                          "stability", "ci_low", "ci_high", "keep_call", "p_value", "q_value", "is_stable"))
empty_table <- function(cols) setnames(data.table(matrix(character(), 0, length(cols))), cols)

no_gd <- nrow(input) == 0 || input[, all(num_gds == 0)]
if (no_gd) {
  message("PGDD: no sample has a GD; writing No GD rows (PGDD itself returns nothing in this case)")
  samples <- unique(input[, .(sample_id, tumour_id)])
  tumours <- unique(samples$tumour_id)
  per_tumour <- if (length(tumours)) data.table(tumour_id = tumours, First_GD = "No GD", Second_GD = "No GD",
    num_first_gd = 0L, num_second_gd = 0L, First_GD_homogen = TRUE, Second_GD_homogen = TRUE,
    GD_status_homogen = TRUE, GD_statuses = "0", frac_0_gd_regions = 1, frac_1_gd_regions = 0,
    frac_2_gd_regions = 0, num_clonal_gds = 0L, num_subclonal_gds = 0L) else empty_table(HEADERS$GDs_per_tumour)
  per_region <- if (nrow(samples)) samples[, .(sample_id, tumour_id, gd_clusters = "", num_gds = 0L, gd_events = "",
    First_GD = "No GD", Second_GD = "No GD")] else empty_table(HEADERS$GDs_per_region)
  out <- list(GDs_per_tumour = per_tumour, GDs_per_region = per_region,
              GDs_events = empty_table(HEADERS$GDs_events), mut_counts = empty_table(HEADERS$mut_counts),
              unstable_calls_GDs_per_tumour = empty_table(HEADERS$GDs_per_tumour),
              unstable_calls_GDs_per_region = empty_table(HEADERS$GDs_per_region),
              unstable_calls_GDs_events = empty_table(HEADERS$GDs_events),
              unstable_calls_mut_counts = empty_table(HEADERS$mut_counts),
              bootstrap_discovery = empty_table(HEADERS$bootstrap_discovery),
              bootstrap_discovery_settings = data.table(n_boot = settings$n_boot, seed = settings$seed,
                stability_threshold = settings$stability_threshold, stability_alpha = settings$stability_alpha,
                stability_use_binom_test = settings$stability_use_binom_test,
                stability_adjust_method = settings$stability_adjust_method,
                bootstrap_unit = "tumour_id x sample_id x cluster_id", phase = "not_run_no_gd_samples",
                check_phase_requires_stable_anchor_in_cluster = TRUE))
} else {
  out <- do.call(detect_par_gd, c(list(input = input), settings))
}

files <- c(GDs_per_tumour = "pgdd_GDs_per_tumour.tsv",
           GDs_per_region = "pgdd_GDs_per_region.tsv",
           GDs_events = "pgdd_GD_events.tsv",
           mut_counts = "pgdd_mut_counts.tsv",
           unstable_calls_GDs_per_tumour = "unstable_calls_pgdd_GDs_per_tumour.tsv",
           unstable_calls_GDs_per_region = "unstable_calls_pgdd_GDs_per_region.tsv",
           unstable_calls_GDs_events = "unstable_calls_pgdd_GD_events.tsv",
           unstable_calls_mut_counts = "unstable_calls_pgdd_mut_counts.tsv",
           bootstrap_discovery = "pgdd_bootstrap_discovery_summary.tsv",
           bootstrap_discovery_settings = "pgdd_bootstrap_discovery_settings.tsv")
for (el in names(files)) {
  if (!is.null(out[[el]])) fwrite(out[[el]], file.path(out_dir, files[[el]]), sep = "\t")
}
writeLines(toJSON(list(settings = settings, pgdd_src = pgdd_src, no_gd_samples = no_gd,
                       pgdd_commit = tryCatch(system(paste("git -C", shQuote(pgdd_src), "rev-parse --short HEAD"),
                                                     intern = TRUE), error = function(e) NA)),
                  auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "run.json"))
message("wrote ", out_dir)
