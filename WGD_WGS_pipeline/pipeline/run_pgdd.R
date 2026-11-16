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
out <- do.call(detect_par_gd, c(list(input = input), settings))

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
writeLines(toJSON(list(settings = settings, pgdd_src = pgdd_src,
                       pgdd_commit = tryCatch(system(paste("git -C", shQuote(pgdd_src), "rev-parse --short HEAD"),
                                                     intern = TRUE), error = function(e) NA)),
                  auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "run.json"))
message("wrote ", out_dir)
