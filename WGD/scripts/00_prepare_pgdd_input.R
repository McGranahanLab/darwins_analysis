suppressPackageStartupMessages({
  library(data.table)
})

# ---- verbose diagnostics helpers -------------------------
check_file_readable = function(file_path, label, tumour_id) {
  if (!file.exists(file_path)) {
    warning("  [", tumour_id, "] Missing ", label, ": ", file_path)
    return(FALSE)
  }

  file_info = file.info(file_path)
  if (is.na(file_info$size) || file_info$size <= 0) {
    warning("  [", tumour_id, "] Empty ", label, ": ", file_path)
    return(FALSE)
  }

  message("  [", tumour_id, "] OK ", label, " (", file_info$size, " bytes): ", file_path)
  return(TRUE)
}

print_dt_summary = function(dt, label, tumour_id, n = 5L) {
  if (is.null(dt)) {
    warning("  [", tumour_id, "] ", label, " is NULL")
    return(invisible(NULL))
  }

  if (!is.data.table(dt)) {
    warning("  [", tumour_id, "] ", label, " is not a data.table (class: ",
            paste(class(dt), collapse = ", "), ")")
    return(invisible(NULL))
  }

  message("  [", tumour_id, "] ", label,
          " -> rows=", nrow(dt),
          ", cols=", ncol(dt))
  if (ncol(dt) > 0) {
    message("  [", tumour_id, "] ", label,
            " columns: ", paste(names(dt), collapse = ", "))
  }

  if (nrow(dt) > 0) {
    print(head(dt, n))
  } else {
    message("  [", tumour_id, "] ", label, " is empty")
  }
}

require_columns = function(dt, required_cols, label, tumour_id) {
  missing_cols = setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    warning("  [", tumour_id, "] Missing required column(s) in ", label, ": ",
            paste(missing_cols, collapse = ", "))
    return(FALSE)
  }
  message("  [", tumour_id, "] ", label, " contains required columns")
  return(TRUE)
}

diagnose_sample_id_alignment = function(muts, num_gds_dt, tumour_id) {
  skip = FALSE
  mut_samples = unique(as.character(muts$SAMPLE))
  cn_samples  = unique(as.character(num_gds_dt$SAMPLE))

  mut_empty = sum(is.na(mut_samples) | trimws(mut_samples) == "")
  cn_empty  = sum(is.na(cn_samples) | trimws(cn_samples) == "")
  message("  [", tumour_id, "] SAMPLE summary: muts unique=", length(mut_samples),
          " (empty/NA=", mut_empty, "), num_gds unique=", length(cn_samples),
          " (empty/NA=", cn_empty, ")")

  mut_samples_clean = unique(mut_samples[!is.na(mut_samples) & trimws(mut_samples) != ""])
  cn_samples_clean  = unique(cn_samples[!is.na(cn_samples) & trimws(cn_samples) != ""])

  overlap_exact = intersect(mut_samples_clean, cn_samples_clean)
  overlap_trimmed = intersect(trimws(mut_samples_clean), trimws(cn_samples_clean))
  overlap_upper_trim = intersect(toupper(trimws(mut_samples_clean)), toupper(trimws(cn_samples_clean)))

  message("  [", tumour_id, "] SAMPLE overlap: exact=", length(overlap_exact),
          ", trimmed=", length(overlap_trimmed),
          ", case-insensitive-trimmed=", length(overlap_upper_trim))

  if (length(overlap_exact) == 0 && length(overlap_trimmed) > 0) {
    warning("  [", tumour_id, "] No exact SAMPLE overlap, but overlap appears after trimws(). This suggests leading/trailing whitespace mismatch.")
  }
  if (length(overlap_exact) == 0 && length(overlap_upper_trim) > 0) {
    warning("  [", tumour_id, "] No exact SAMPLE overlap, but overlap appears after case-insensitive comparison. This suggests case mismatch.")
  }

  only_in_muts = setdiff(mut_samples_clean, cn_samples_clean)
  only_in_cn   = setdiff(cn_samples_clean, mut_samples_clean)
  unique_to_either = length(only_in_muts) + length(only_in_cn)
  if (length(only_in_muts) > 0) {
    message("  [", tumour_id, "] SAMPLE present only in mutations (first 10): ",
            paste(head(only_in_muts, 10), collapse = ", "))
  }
  if (length(only_in_cn) > 0) {
    message("  [", tumour_id, "] SAMPLE present only in num_gds/CN (first 10): ",
            paste(head(only_in_cn, 10), collapse = ", "))
  }
  # throw error and exit if at least one sample is not matching:
  if (unique_to_either > 0) {
    print('-----------------')
    print(unique_to_either)
    message("  [", tumour_id, "] SAMPLE IDs do not match between mutations and CN data. This will lead to NA num_gds values after merge. Please check the above diagnostics for clues on the mismatch.")
    skip = TRUE
  }
  return(skip)
}

parse_boolean = function(arg) {
  if (is.logical(arg)) return(arg)
  if (is.character(arg)) {
    arg_lower = tolower(arg)
    if (arg_lower %in% c("true", "t", "1")) return(TRUE)
    if (arg_lower %in% c("false", "f", "0")) return(FALSE)
  }
  stop("Invalid boolean argument: ", arg)
}
# ---- pars arguments ------------------------------------
if ((length(args) == 0)) {
  # assume we are in interactive debugging mode and set default arguments for testing
  args = c(conipher_dir, msai_dir, output_dir, "TRUE")
} else {
  args = commandArgs(trailingOnly = TRUE)
}


conipher_dir = args[1]
msai_dir    = args[2]
output_dir  = args[3]
remove_nontree_cluesters = ifelse(length(args) >= 4, parse_boolean(args[4]), TRUE)
# printed parsed arguments for confirmation
message("Parsed arguments:")
message("  conipher_dir: ", conipher_dir)
message("  msai_dir:     ", msai_dir)
message("  output_dir:   ", output_dir)
message("  remove_nontree_cluesters: ", remove_nontree_cluesters)
# --------------------------------------------------------

# ---- ad-hoc params -------------------------------------
selected_tumours = NULL # set to NULL to run all
# --------------------------------------------------------

# ---- helper: parse allTrees.txt and return clonal cluster --
# The clonal cluster is the root node: present as a parent
# but never appearing as a child in any edge of the tree.
get_clonal_cluster = function(tree_file, tumour_id) {
  tree_lines = readLines(tree_file)
  # drop comment/header lines (starting with #) and blank lines
  tree_lines = tree_lines[!grepl("^#", tree_lines) & nzchar(trimws(tree_lines))]
  if (length(tree_lines) == 0) stop("No edges found in ", tree_file)
  if (length(tree_lines) == 1) stop("Single cluster tree found in ", tree_file)
  edges = do.call(rbind, strsplit(tree_lines, "\t"))
  parents  = as.integer(edges[, 1])
  children = as.integer(edges[, 2])
  root = setdiff(unique(parents), unique(children))
  if (length(root) == 0) stop("Could not identify root in ", tree_file)
  if (length(root) > 1) {
    warning("Multiple root candidates for ", tumour_id, ": ",
            paste(root, collapse = ", "), ". Using the smallest.")
    root = min(root)
  }
  return(root)
}


# ---- helper: parse treeTable.tsv and return a list of tree clusters --
get_tree_clusters = function(tree_file, tumour_id) {
  tree_lines = readLines(tree_file)
  # drop comment/header lines (starting with #) and blank lines
  tree_lines = tree_lines[!grepl("^#", tree_lines) & nzchar(trimws(tree_lines))]
  if (length(tree_lines) == 0) stop("No edges found in ", tree_file)
  if (length(tree_lines) == 1) stop("Single cluster tree found in ", tree_file)
  edges = do.call(rbind, strsplit(tree_lines, "\t"))
  parents  = as.integer(edges[, 1])
  children = as.integer(edges[, 2])
  list_of_clusters = unique(c(parents, children))
  # ensure the list contains integers:
  list_of_clusters = as.integer(list_of_clusters)
  return(list_of_clusters)
}


# ---- helper: compute num_gds per sample from CN segments ---
# Thresholds follow README recommendations:
#   num_gds = 3  if >=50% of genome (by bp) has round(MajCN) >= 4
#   num_gds = 2  if >=50% of genome (by bp) has round(MajCN) >= 3
#   num_gds = 1  if >=50% of genome (by bp) has round(MajCN) >= 2
#   num_gds = 0  otherwise
# MajCN = max(COPY_NUMBER_A, COPY_NUMBER_B) per segment
# Additionally compute the confidence level so that we can distinguish between samples were GD called on values close to threshold vs close to theoretical expected
compute_level_confidence = function(cnA, cnB, seg_len, target_cn) {
  valid = is.finite(cnA) & is.finite(cnB) & is.finite(seg_len) & seg_len > 0
  if (!any(valid)) return(NA_real_)

  # cap at each GD level to avoid over-penalizing segments with very high CN
  a = pmin(pmax(cnA[valid], 0), target_cn)
  b = pmin(pmax(cnB[valid], 0), target_cn)
  w = seg_len[valid]

  # 0 at previous baseline (target/2 | target/2), 1 at perfect (target | target)
  weighted_deficit = sum(((target_cn - a) + (target_cn - b)) * w) / sum(w)
  score = 1 - (weighted_deficit / target_cn)

  return(pmax(0, pmin(1, score)))
}

compute_num_gds = function(cn) {
  cn = copy(cn)
  cn[, seg_len := ENDPOS - STARTPOS + 1L]
  cn[, MajCN   := pmax(COPY_NUMBER_A, COPY_NUMBER_B)]

  num_gds_dt = cn[, {
    total_len    = sum(seg_len, na.rm = TRUE)
    frac_maj_ge2 = sum(seg_len[round(MajCN) >= 2L], na.rm = TRUE) / total_len
    frac_maj_ge3 = sum(seg_len[round(MajCN) >= 3L], na.rm = TRUE) / total_len
    frac_maj_ge4 = sum(seg_len[round(MajCN) >= 4L], na.rm = TRUE) / total_len

    num_gds = fcase(
      frac_maj_ge4 >= 0.5, 3L,
      frac_maj_ge3 >= 0.5, 2L,
      frac_maj_ge2 >= 0.5, 1L,
      default = 0L
    )

    conf_gd1 = compute_level_confidence(COPY_NUMBER_A, COPY_NUMBER_B, seg_len, target_cn = 2)
    conf_gd2 = compute_level_confidence(COPY_NUMBER_A, COPY_NUMBER_B, seg_len, target_cn = 4)
    conf_gd3 = compute_level_confidence(COPY_NUMBER_A, COPY_NUMBER_B, seg_len, target_cn = 8)

    .(
      frac_maj_ge2 = frac_maj_ge2,
      frac_maj_ge3 = frac_maj_ge3,
      frac_maj_ge4 = frac_maj_ge4,
      num_gds = num_gds,
      conf_gd1 = conf_gd1,
      conf_gd2 = conf_gd2,
      conf_gd3 = conf_gd3
    )
  }, by = SAMPLE]

  return(num_gds_dt)
}

# ---- discover tumour directories ---------------------------
tumour_dirs = list.dirs(conipher_dir, full.names = TRUE, recursive = FALSE)
if (length(tumour_dirs) == 0) stop("No tumour directories found in: ", conipher_dir)
message("Found ", length(tumour_dirs), " tumour(s): ",
        paste(basename(tumour_dirs), collapse = ", "))

# If selected_tumours is set, filter to those and warn if any are missing
if (!is.null(selected_tumours)) {
  missing_tumours = setdiff(selected_tumours, basename(tumour_dirs))
  if (length(missing_tumours) > 0) {
    warning("The following selected tumours were not found and will be skipped: ",
            paste(missing_tumours, collapse = ", "))
  }
  tumour_dirs = tumour_dirs[basename(tumour_dirs) %in% selected_tumours]
  message("Processing ", length(tumour_dirs), " tumour(s) after filtering: ",
          paste(basename(tumour_dirs), collapse = ", "))
}

# ---- process each tumour -----------------------------------

all_pgdd_input = lapply(tumour_dirs, function(tumour_dir) {
  tumour_id = basename(tumour_dir)
  message("\nProcessing: ", tumour_id)
  msai_tumour_dir = file.path(msai_dir, tumour_id)
  # --- check required files exist --
  mut_file  = file.path(tumour_dir, paste0(tumour_id, ".SCoutput.CLEAN.tsv"))
  tree_table_file = file.path(tumour_dir, "treeTable.tsv")
  cn_file   = file.path(msai_tumour_dir, paste0(tumour_id, "_copynumber_MSAI_adjusted_segments_tree_input.tsv"))
  tree_file = file.path(tumour_dir, "allTrees.txt")

  files_ok = c(
    check_file_readable(mut_file, "mutation file", tumour_id),
    check_file_readable(cn_file, "CN segment file", tumour_id),
    check_file_readable(tree_file, "tree file", tumour_id),
    check_file_readable(tree_table_file, "tree table file", tumour_id)
  )
  if (!all(files_ok)) {
    warning("Skipping ", tumour_id, " due to missing/empty required input file(s)")
    return(NULL)
  }

  # --- mutations ---
  muts = tryCatch(
    fread(mut_file),
    error = function(e) {
      warning("  [", tumour_id, "] Failed reading mutation file: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(muts)) return(NULL)
  print_dt_summary(muts, "muts (raw)", tumour_id)

  required_mut_cols = c(
    "CASE_ID", "CHR", "POS", "REF", "ALT", "CLUSTER", "SAMPLE",
    "MUT_COPY", "COPY_NUMBER_A", "COPY_NUMBER_B"
  )
  if (!require_columns(muts, required_mut_cols, "muts (raw)", tumour_id)) {
    warning("Skipping ", tumour_id, " due to missing mutation columns")
    return(NULL)
  }
  if (nrow(muts) == 0) {
    warning("  [", tumour_id, "] Mutation table is empty")
  }

  # --- tree table ---
  tree_table = tryCatch(
    fread(tree_table_file),
    error = function(e) {
      warning("  [", tumour_id, "] Failed reading tree table file: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(tree_table)) return(NULL)
  print_dt_summary(tree_table, "tree_table (raw)", tumour_id)

  # filter muts by tree table mutations (i.e. mutations which made it to the tree an passed all the conipher filters)
  tree_table$mutation_id_sample = paste(tree_table$mutation_id, tree_table$SAMPLE, sep = ":")
  muts$mutation_id_sample = paste(muts$mutation_id, muts$SAMPLE, sep = ":")

  tree_table_merge_sub_df = tree_table[, .(mutation_id_sample, treeCLUSTER)]
  muts = merge(
    muts,
    tree_table_merge_sub_df,
    by = "mutation_id_sample",
    all.x = FALSE,
    all.y = FALSE
  )
  muts$CLUSTER = muts$treeCLUSTER
  
  # --- clonal cluster from tree ---
  clonal_cluster = tryCatch(
    get_clonal_cluster(tree_file, tumour_id),
    error = function(e) {
      warning("  [", tumour_id, "] Failed parsing tree file: ", conditionMessage(e))
      return(NA_integer_)
    }
  )
  if (is.na(clonal_cluster)) {
    warning("Skipping ", tumour_id, " due to invalid clonal cluster")
    return(NULL)
  }
  message("  Clonal cluster: ", clonal_cluster)


  # --- tree cluster from tree ---
  tree_clusters = tryCatch(
    get_tree_clusters(tree_file, tumour_id),
    error = function(e) {
      warning("  [", tumour_id, "] Failed parsing tree file: ", conditionMessage(e))
      return(NA_integer_)
    }
  )
  if (any(is.na(tree_clusters))) {
    warning("Skipping ", tumour_id, " due to invalid tree clusters")
    return(NULL)
  }
  message("  Tree clusters: ", paste(tree_clusters, collapse = ", "))


  # --- num_gds from CN segments ---
  cn = tryCatch(
    fread(cn_file),
    error = function(e) {
      warning("  [", tumour_id, "] Failed reading CN segment file: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(cn)) return(NULL)
  print_dt_summary(cn, "cn (raw)", tumour_id)

  required_cn_cols = c("SAMPLE", "STARTPOS", "ENDPOS", "COPY_NUMBER_A", "COPY_NUMBER_B")
  if (!require_columns(cn, required_cn_cols, "cn (raw)", tumour_id)) {
    warning("Skipping ", tumour_id, " due to missing CN columns")
    return(NULL)
  }
  if (nrow(cn) == 0) {
    warning("  [", tumour_id, "] CN segment table is empty")
  }
  
  num_gds_dt = tryCatch(
    compute_num_gds(cn),
    error = function(e) {
      warning("  [", tumour_id, "] Failed computing num_gds: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(num_gds_dt)) return(NULL)
  print_dt_summary(num_gds_dt, "num_gds_dt", tumour_id)

  if (!require_columns(num_gds_dt, c("SAMPLE", "num_gds"), "num_gds_dt", tumour_id)) {
    warning("Skipping ", tumour_id, " due to malformed num_gds table")
    return(NULL)
  }

  message("  num_gds per sample:")
  print(num_gds_dt)

  skip = diagnose_sample_id_alignment(muts, num_gds_dt, tumour_id)
  if (skip) {
    warning("Skipping ", tumour_id, " due to SAMPLE ID mismatch between mutations and CN data.")
    return(NULL)
  }
  # --- merge num_gds into mutations ---
  muts = merge(muts, num_gds_dt, by = "SAMPLE", all.x = TRUE)
  print_dt_summary(muts, "muts (after merge)", tumour_id)
  missing_gds = muts[is.na(num_gds), unique(SAMPLE)]
  if (length(missing_gds) > 0)
    warning("  No CN data for sample(s): ", paste(missing_gds, collapse = ", "))
  if (nrow(muts) > 0 && all(is.na(muts$num_gds))) {
    warning("  [", tumour_id, "] All num_gds values are NA after merge. This usually indicates SAMPLE ID mismatch between mutation and CN inputs.")
    message("  [", tumour_id, "] First 10 SAMPLE IDs in merged table: ",
            paste(head(unique(as.character(muts$SAMPLE)), 10), collapse = ", "))
  }

  # --- build PGDD input table ---
  # MajCN at mutation locus = max(COPY_NUMBER_A, COPY_NUMBER_B) from mutations file

  if (remove_nontree_cluesters) {
    muts = muts[CLUSTER %in% tree_clusters]
    message("  [", tumour_id, "] Removed non-tree clusters. Remaining clusters: ",
            paste(sort(unique(muts$CLUSTER)), collapse = ", "))
  }

  pgdd_input = muts[, .(
    tumour_id         = CASE_ID,
    chromosome        = paste0("chr", CHR),
    position          = POS,
    ref               = REF,
    alt               = ALT,
    cluster_id        = CLUSTER,
    is_clonal_cluster = (CLUSTER == clonal_cluster),
    sample_id         = SAMPLE,
    mut_cpn           = MUT_COPY,
    MajCN             = pmax(COPY_NUMBER_A, COPY_NUMBER_B),
    num_gds           = num_gds,
    frac_maj_ge2      = frac_maj_ge2,
    frac_maj_ge3      = frac_maj_ge3,
    frac_maj_ge4      = frac_maj_ge4,
    conf_gd1          = conf_gd1,
    conf_gd2          = conf_gd2,
    conf_gd3          = conf_gd3,
    ploidy            = if ("PLOIDY" %in% names(muts)) PLOIDY else NA_real_
  )]


  print_dt_summary(pgdd_input, "pgdd_input (per tumour)", tumour_id)
  if (nrow(pgdd_input) == 0) {
    warning("  [", tumour_id, "] pgdd_input has 0 rows after transformation")
    return(NULL)
  }

  message("  Samples: ", pgdd_input[, length(unique(sample_id))],
          "  |  Mutations: ", nrow(pgdd_input),
          "  |  Clusters: ", pgdd_input[, length(unique(cluster_id))])
  return(pgdd_input)
})
browser()
# ---- combine and write output ------------------------------
all_pgdd_input = rbindlist(Filter(Negate(is.null), all_pgdd_input))
message("Combined non-NULL tumour tables: ", nrow(all_pgdd_input), " rows")
if (nrow(all_pgdd_input) > 0) {
  print_dt_summary(all_pgdd_input, "all_pgdd_input (combined)", "GLOBAL")
}

if (nrow(all_pgdd_input) == 0) stop("No PGDD input rows produced — check warnings above.")

# ensure uniform convention for tumour_id:
all_pgdd_input[, tumour_id := gsub("_", "-Tumour", tumour_id)]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}
output_file = file.path(output_dir, "pgdd_input_all_tumours.tsv")
fwrite(all_pgdd_input, output_file, sep = "\t")

message("\n=== Summary ===")
message("Output: ",   output_file)
message("Tumours: ",  all_pgdd_input[, length(unique(tumour_id))])
message("Samples: ",  all_pgdd_input[, length(unique(sample_id))])
message("Mutations: ", nrow(all_pgdd_input))
message("Done.")
