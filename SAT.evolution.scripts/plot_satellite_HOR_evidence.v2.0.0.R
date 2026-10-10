#!/usr/bin/env Rscript

# ==============================================================================
# SAT113/SAT483 higher-order repeat (HOR) candidate scan and evidence plots
#
# This script reads the output of the SAT113/SAT483 variant-origin workflow.  In
# automatic mode it scans every source-array block, reports statistically
# supported HOR-like periodicity, and produces one evidence page per candidate.
# Manual mode retains the block/config interface of the SAT113 v1.0.2 script.
#
# Direct evidence shown on every page:
#   A. Ordered recurrence of the main-threshold variant clusters
#   B. Pairwise monomer identity heatmap (when the fixed-reference alignment is
#      available), otherwise a clearly labelled multi-threshold cluster matrix
#   C. Lag concordance at 80%, 85%, 90%, and 95% clustering thresholds
#   D. Same-HOR-phase versus other-lag monomer similarity
#
# Main input may be a complete run directory, its results/ directory, or a ZIP.
# The script accepts --input, --primary-dir, or one positional input path.
#
# Automatic SAT483 scan and candidate plots:
#   Rscript plot_satellite_HOR_evidence.v2.0.0.R \
#     --input SAT483.results.zip --mode auto --dpi 600
#
# Select individual blocks (repeat --block as needed):
#   Rscript plot_satellite_HOR_evidence.v2.0.0.R RESULTS_OR_ZIP \
#     --mode manual \
#     --block "FCA,ChrF1__28252,2" \
#     --block "PPA,ChrE1__21831,66" \
#     --block "PTI,ChrA1__694,81"
#
# A --block value has these comma-separated fields:
#   SPECIES,SOURCE_ARRAY_KEY,EXPECTED_PERIOD[,COPY_START[,COPY_END[,LABEL]]]
# Use "auto" for EXPECTED_PERIOD when it is unknown.
#
# For many/custom selections, use --config.  Required column: species.  Supply
# source_array_key, or supply region_key and leave source_array_key blank/ALL to
# expand the region into separate source blocks.  Optional columns:
#   expected_period, copy_start, copy_end, label, enabled
#
# Required plotting packages: ggplot2, patchwork, scales
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

SCRIPT_VERSION <- "2026-09-22.v2.0.0"
SCRIPT_NAME <- "plot_satellite_HOR_evidence.v2.0.0.R"

usage <- function(status = 0L) {
  cat(
    SCRIPT_NAME, " ", SCRIPT_VERSION, "\n\n",
    "Usage:\n",
    "  Rscript ", SCRIPT_NAME, " INPUT [options]\n",
    "  Rscript ", SCRIPT_NAME, " --input INPUT [options]\n\n",
    "Input:\n",
    "  INPUT                        Run root, results/ directory, or ZIP\n",
    "  --input PATH                 Named input form\n",
    "  --primary-dir PATH           Alias of --input\n",
    "  --alignment FASTA            Optional fixed-reference haplotype alignment\n",
    "                               (auto-detected in a complete run directory)\n\n",
    "Family and mode:\n",
    "  --satellite auto|SAT113|SAT483  Family label [auto]\n",
    "  --mode auto|manual            Scan all blocks or inspect selected blocks [auto]\n\n",
    "Automatic candidate scan:\n",
    "  --min-copies INT             Minimum copies required to test a block [16]\n",
    "  --min-period INT             Smallest tested HOR period [2]\n",
    "  --max-period INT             Largest tested HOR period [120]\n",
    "  --min-units FLOAT            Minimum putative units represented [2]\n",
    "  --min-lag-pairs INT          Minimum valid pairs at a lag [8]\n",
    "  --min-pair-fraction FLOAT    Minimum pair count as block fraction [0.25]\n",
    "  --support-excess FLOAT       Per-threshold support cutoff [0.10]\n",
    "  --screen-excess FLOAT        Permutation-screen effect cutoff [0.10]\n",
    "  --screen-support INT         Thresholds required for permutations [2]\n",
    "  --high-fdr FLOAT             HIGH FDR cutoff [0.05]\n",
    "  --high-excess FLOAT          HIGH robust-effect cutoff [0.20]\n",
    "  --high-support INT           Thresholds required for HIGH [3]\n",
    "  --moderate-fdr FLOAT         MODERATE FDR cutoff [0.10]\n",
    "  --moderate-excess FLOAT      MODERATE robust-effect cutoff [0.10]\n",
    "  --moderate-support INT       Thresholds required for MODERATE [2]\n",
    "  --candidate-calls TEXT       Calls plotted in auto mode [HIGH,MODERATE]\n",
    "                               With --permutations 0, the implicit default is\n",
    "                               EFFECT_ONLY_HIGH,EFFECT_ONLY_MODERATE\n\n",
    "Manual block selection:\n",
    "  --preset strong|extended|all|none  SAT113 legacy candidate presets [none]\n",
    "  --block SPEC                 Repeatable inline block selector\n",
    "  --config TSV                 Block-selection table\n",
    "  --list-blocks                Write available blocks and exit\n",
    "  --min-block-copies INT       Minimum copies when a whole region is expanded [10]\n\n",
    "Analysis and output:\n",
    "  --identity FLOAT             Main cluster identity [0.85]\n",
    "  --max-lag INT                Manual-plot maximum lag [120]\n",
    "  --permutations INT           Joint-row max-lag permutations [999]\n",
    "                               Use 99 for exploration; 0 gives effect-only calls\n",
    "  --max-matrix-copies INT      Maximum copies rendered in a matrix [600]\n",
    "  --heatmap-min FLOAT          Lower exact-identity colour limit [0.50]\n",
    "  --outdir DIR                 Output directory [<family>_HOR_analysis]\n",
    "  --prefix TEXT                Output filename prefix [<family>]\n",
    "  --dpi INT                    Per-block PNG resolution [600]\n",
    "  --font-family TEXT           Plot font family [sans]\n",
    "  --seed INT                   Random seed [20260731]\n",
    "  --no-png                     Do not write per-block PNG files\n",
    "  --no-pdf                     Do not write the multipage PDF\n",
    "  --no-plots                   Write automatic-scan tables only\n",
    "  --version                    Print version and exit\n",
    "  --help                       Show this help\n\n",
    "Preset strong blocks:\n",
    "  FCA ChrF1__28235 (8), 28244 (8), 28252 (2), 28258 (8),\n",
    "      28259 (21), 28271 (8)\n",
    "  PPA ChrE1__21788 (62), 21808 (65), 21814 (65),\n",
    "      21830 (66), 21831 (66)\n",
    "  PTI ChrA1__694 (81)\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

parse_args <- function(args) {
  opt <- list(
    input = NULL,
    alignment = NULL,
    satellite = "AUTO",
    mode = NULL,
    mode_explicit = FALSE,
    preset = NULL,
    preset_explicit = FALSE,
    blocks = character(),
    config = NULL,
    list_blocks = FALSE,
    min_block_copies = 10L,
    min_copies = 16L,
    min_period = 2L,
    max_period = 120L,
    min_units = 2,
    min_lag_pairs = 8L,
    min_pair_fraction = 0.25,
    support_excess = 0.10,
    screen_excess = 0.10,
    screen_support = 2L,
    fundamental_abs_tolerance = 0.02,
    fundamental_rel_tolerance = 0.05,
    high_fdr = 0.05,
    high_excess = 0.20,
    high_support = 3L,
    moderate_fdr = 0.10,
    moderate_excess = 0.10,
    moderate_support = 2L,
    candidate_calls = "HIGH,MODERATE",
    candidate_calls_explicit = FALSE,
    main_identity = 0.85,
    max_lag = 120L,
    permutations = 999L,
    max_matrix_copies = 600L,
    heatmap_min = 0.50,
    outdir = NULL,
    prefix = NULL,
    dpi = 600L,
    font_family = "sans",
    seed = 20260731L,
    write_png = TRUE,
    write_pdf = TRUE
  )

  value_options <- c(
    "--input", "--primary-dir", "--alignment", "--satellite", "--mode",
    "--preset", "--block", "--config", "--min-block-copies",
    "--min-copies", "--min-period", "--max-period", "--min-units",
    "--min-lag-pairs", "--min-pair-fraction", "--support-excess",
    "--screen-excess", "--screen-support", "--fundamental-abs-tolerance",
    "--fundamental-rel-tolerance", "--high-fdr", "--high-excess",
    "--high-support", "--moderate-fdr", "--moderate-excess",
    "--moderate-support", "--candidate-calls", "--identity", "--max-lag",
    "--permutations", "--max-matrix-copies", "--heatmap-min",
    "--outdir", "--prefix", "--dpi", "--font-family", "--seed"
  )

  i <- 1L
  positional <- character()
  while (i <= length(args)) {
    token <- args[[i]]

    if (token %in% c("--help", "-h")) usage(0L)
    if (token == "--version") {
      cat(SCRIPT_NAME, " ", SCRIPT_VERSION, "\n", sep = "")
      quit(save = "no", status = 0L)
    }
    if (token == "--list-blocks") {
      opt$list_blocks <- TRUE
      i <- i + 1L
      next
    }
    if (token == "--no-png") {
      opt$write_png <- FALSE
      i <- i + 1L
      next
    }
    if (token == "--no-pdf") {
      opt$write_pdf <- FALSE
      i <- i + 1L
      next
    }
    if (token == "--no-plots") {
      opt$write_pdf <- FALSE
      opt$write_png <- FALSE
      i <- i + 1L
      next
    }

    matched_equals <- FALSE
    for (name in value_options) {
      prefix <- paste0(name, "=")
      if (startsWith(token, prefix)) {
        args <- append(args, substring(token, nchar(prefix) + 1L), after = i)
        token <- name
        matched_equals <- TRUE
        break
      }
    }

    if (token %in% value_options) {
      if (i == length(args)) stop("Missing value after ", token)
      value <- args[[i + 1L]]
      if (token %in% c("--input", "--primary-dir")) opt$input <- value
      if (token == "--alignment") opt$alignment <- value
      if (token == "--satellite") opt$satellite <- toupper(value)
      if (token == "--mode") {
        opt$mode <- tolower(value)
        opt$mode_explicit <- TRUE
      }
      if (token == "--preset") {
        opt$preset <- tolower(value)
        opt$preset_explicit <- TRUE
      }
      if (token == "--block") opt$blocks <- c(opt$blocks, value)
      if (token == "--config") opt$config <- value
      if (token == "--min-block-copies") opt$min_block_copies <- as.integer(value)
      if (token == "--min-copies") opt$min_copies <- as.integer(value)
      if (token == "--min-period") opt$min_period <- as.integer(value)
      if (token == "--max-period") opt$max_period <- as.integer(value)
      if (token == "--min-units") opt$min_units <- as.numeric(value)
      if (token == "--min-lag-pairs") opt$min_lag_pairs <- as.integer(value)
      if (token == "--min-pair-fraction") opt$min_pair_fraction <- as.numeric(value)
      if (token == "--support-excess") opt$support_excess <- as.numeric(value)
      if (token == "--screen-excess") opt$screen_excess <- as.numeric(value)
      if (token == "--screen-support") opt$screen_support <- as.integer(value)
      if (token == "--fundamental-abs-tolerance") opt$fundamental_abs_tolerance <- as.numeric(value)
      if (token == "--fundamental-rel-tolerance") opt$fundamental_rel_tolerance <- as.numeric(value)
      if (token == "--high-fdr") opt$high_fdr <- as.numeric(value)
      if (token == "--high-excess") opt$high_excess <- as.numeric(value)
      if (token == "--high-support") opt$high_support <- as.integer(value)
      if (token == "--moderate-fdr") opt$moderate_fdr <- as.numeric(value)
      if (token == "--moderate-excess") opt$moderate_excess <- as.numeric(value)
      if (token == "--moderate-support") opt$moderate_support <- as.integer(value)
      if (token == "--candidate-calls") {
        opt$candidate_calls <- value
        opt$candidate_calls_explicit <- TRUE
      }
      if (token == "--identity") opt$main_identity <- as.numeric(value)
      if (token == "--max-lag") opt$max_lag <- as.integer(value)
      if (token == "--permutations") opt$permutations <- as.integer(value)
      if (token == "--max-matrix-copies") opt$max_matrix_copies <- as.integer(value)
      if (token == "--heatmap-min") opt$heatmap_min <- as.numeric(value)
      if (token == "--outdir") opt$outdir <- value
      if (token == "--prefix") opt$prefix <- value
      if (token == "--dpi") opt$dpi <- as.integer(value)
      if (token == "--font-family") opt$font_family <- value
      if (token == "--seed") opt$seed <- as.integer(value)
      i <- i + 2L
      next
    }

    if (startsWith(token, "-")) stop("Unknown option: ", token)
    positional <- c(positional, token)
    i <- i + 1L
  }

  if (length(positional) > 1L) {
    stop("Unexpected positional arguments: ", paste(positional[-1L], collapse = " "))
  }
  if (is.null(opt$input) && length(positional) == 1L) opt$input <- positional[[1L]]
  if (!is.null(opt$input) && length(positional) == 1L && positional[[1L]] != opt$input) {
    stop("Input was supplied both positionally and by option with different values.")
  }
  if (is.null(opt$input)) stop("No input supplied. Use --input PATH or one positional PATH.")

  finite_positive_integer <- function(x, name, allow_zero = FALSE) {
    lower_ok <- if (allow_zero) x >= 0 else x >= 1
    if (length(x) != 1L || !is.finite(x) || !lower_ok) stop(name, " is invalid: ", x)
  }
  finite_positive_integer(opt$min_block_copies, "--min-block-copies")
  finite_positive_integer(opt$min_copies, "--min-copies")
  finite_positive_integer(opt$min_period, "--min-period")
  finite_positive_integer(opt$max_period, "--max-period")
  finite_positive_integer(opt$min_lag_pairs, "--min-lag-pairs")
  finite_positive_integer(opt$screen_support, "--screen-support")
  finite_positive_integer(opt$high_support, "--high-support")
  finite_positive_integer(opt$moderate_support, "--moderate-support")
  finite_positive_integer(opt$max_lag, "--max-lag")
  finite_positive_integer(opt$permutations, "--permutations", allow_zero = TRUE)
  finite_positive_integer(opt$max_matrix_copies, "--max-matrix-copies")
  finite_positive_integer(opt$dpi, "--dpi")
  if (!is.finite(opt$main_identity) || opt$main_identity <= 0 || opt$main_identity > 1) {
    stop("--identity must be in (0, 1].")
  }
  if (!is.finite(opt$heatmap_min) || opt$heatmap_min < 0 || opt$heatmap_min >= 1) {
    stop("--heatmap-min must be in [0, 1).")
  }
  if (!opt$satellite %in% c("AUTO", "SAT113", "SAT483")) {
    stop("--satellite must be auto, SAT113, or SAT483.")
  }
  if (!is.null(opt$mode) && !opt$mode %in% c("auto", "manual")) {
    stop("--mode must be auto or manual.")
  }
  if (opt$max_period < opt$min_period) stop("--max-period must be >= --min-period.")
  if (!is.finite(opt$min_units) || opt$min_units < 2) stop("--min-units must be >= 2.")
  bounded <- c(
    min_pair_fraction = opt$min_pair_fraction,
    support_excess = opt$support_excess,
    screen_excess = opt$screen_excess,
    fundamental_abs_tolerance = opt$fundamental_abs_tolerance,
    fundamental_rel_tolerance = opt$fundamental_rel_tolerance,
    high_fdr = opt$high_fdr,
    high_excess = opt$high_excess,
    moderate_fdr = opt$moderate_fdr,
    moderate_excess = opt$moderate_excess
  )
  if (any(!is.finite(bounded)) || any(bounded < 0) || any(bounded > 1)) {
    stop("Scan proportions/FDR/effect thresholds must lie in [0, 1].")
  }
  if (opt$high_fdr > opt$moderate_fdr) stop("--high-fdr must not exceed --moderate-fdr.")
  opt$candidate_calls <- unique(trimws(strsplit(toupper(opt$candidate_calls), ",", fixed = TRUE)[[1L]]))
  if (opt$permutations == 0L && !opt$candidate_calls_explicit) {
    opt$candidate_calls <- c("EFFECT_ONLY_HIGH", "EFFECT_ONLY_MODERATE")
  }
  valid_calls <- c("HIGH", "MODERATE", "EFFECT_ONLY_HIGH", "EFFECT_ONLY_MODERATE")
  if (!length(opt$candidate_calls) || any(!opt$candidate_calls %in% valid_calls)) {
    stop("--candidate-calls contains unsupported values: ", paste(setdiff(opt$candidate_calls, valid_calls), collapse = ", "))
  }

  custom_present <- length(opt$blocks) > 0L || !is.null(opt$config) || opt$preset_explicit
  if (is.null(opt$mode)) opt$mode <- if (custom_present) "manual" else "auto"
  if (opt$mode == "auto" && custom_present) {
    stop("--mode auto cannot be combined with --block, --config, or --preset; use --mode manual.")
  }

  opt
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      label, " lacks required columns: ", paste(missing, collapse = ", "),
      "\nAvailable columns: ", paste(names(x), collapse = ", ")
    )
  }
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^[_ .-]+|[_ .-]+$", "", x)
  ifelse(nchar(x) > 0L, x, "block")
}

infer_satellite_family <- function(copy_lengths, input_path = "") {
  lengths <- safe_numeric(copy_lengths)
  lengths <- lengths[is.finite(lengths) & lengths > 0]
  from_length <- NA_character_
  if (length(lengths)) {
    centre <- stats::median(lengths)
    distances <- c(SAT113 = abs(centre - 113), SAT483 = abs(centre - 483))
    winner <- names(which.min(distances))[[1L]]
    # Reject unrelated repeat families rather than silently attaching a wrong
    # biological label.  The broad tolerance allows partial boundary calls.
    if (distances[[winner]] <= max(40, 0.30 * as.numeric(sub("SAT", "", winner)))) {
      from_length <- winner
    }
  }
  path_upper <- toupper(as.character(input_path))
  path_hits <- c(SAT113 = grepl("SAT113", path_upper, fixed = TRUE),
                 SAT483 = grepl("SAT483", path_upper, fixed = TRUE))
  from_path <- if (sum(path_hits) == 1L) names(path_hits)[path_hits][[1L]] else NA_character_
  if (!is.na(from_length) && !is.na(from_path) && from_length != from_path) {
    warning("Satellite name inferred from copy length (", from_length,
            ") conflicts with input path (", from_path, "); using copy length.")
  }
  if (!is.na(from_length)) return(from_length)
  if (!is.na(from_path)) return(from_path)
  stop("Could not infer SAT113 or SAT483. Supply --satellite explicitly.")
}

resolve_satellite_family <- function(requested, copy_lengths, input_path = "") {
  requested <- toupper(requested)
  if (requested == "AUTO") return(infer_satellite_family(copy_lengths, input_path))
  inferred <- tryCatch(
    infer_satellite_family(copy_lengths, input_path),
    error = function(e) NA_character_
  )
  if (!is.na(inferred) && inferred != requested) {
    warning("--satellite ", requested, " conflicts with copy-length inference ",
            inferred, ". Explicit label retained.")
  }
  requested
}

nonempty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x))) &
    !toupper(trimws(as.character(x))) %in% c("NA", "NONE", "ALL", "*")
}

resolve_input <- function(path) {
  if (dir.exists(path)) {
    return(list(root = normalizePath(path, mustWork = TRUE), temporary = NULL))
  }
  if (file.exists(path) && grepl("\\.zip$", path, ignore.case = TRUE)) {
    tmp <- tempfile(pattern = "satellite_HOR_input_")
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(path, exdir = tmp)
    return(list(root = normalizePath(tmp, mustWork = TRUE), temporary = tmp))
  }
  stop("Input must be a directory or ZIP file: ", path)
}

find_exact_file <- function(root, basename_wanted, required = TRUE) {
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  hits <- candidates[basename(candidates) %in% basename_wanted]
  if (length(hits) == 0L) {
    if (required) {
      stop("Cannot find ", paste(basename_wanted, collapse = " or "), " under: ", root)
    }
    return(NULL)
  }
  result_bonus <- ifelse(grepl("/results/", gsub("\\\\", "/", hits), fixed = TRUE), 0, 100000)
  subtype_bonus <- ifelse(grepl("/results/subtypes/", gsub("\\\\", "/", hits), fixed = TRUE), 0, 1000)
  hits[order(result_bonus + subtype_bonus + nchar(hits))][[1L]]
}

thresholds <- c(0.80, 0.85, 0.90, 0.95)
threshold_columns <- c(
  "0.8" = "variant_id_0_8",
  "0.85" = "variant_id_0_85",
  "0.9" = "variant_id_0_9",
  "0.95" = "variant_id_0_95"
)
threshold_labels <- c("0.8" = "80%", "0.85" = "85%", "0.9" = "90%", "0.95" = "95%")
threshold_colours <- c("80%" = "#CC79A7", "85%" = "#D55E00", "90%" = "#009E73", "95%" = "#0072B2")
species_colours <- c("FCA" = "#D55E00", "PPA" = "#0072B2", "PTI" = "#009E73")
heatmap_colours <- c("#FFFFFF", "#FFFFB2", "#FED976", "#FD8D3C", "#BD0026")

threshold_key <- function(x) {
  keys <- as.numeric(names(threshold_columns))
  j <- which.min(abs(keys - x))
  if (length(j) == 0L || abs(keys[[j]] - x) > 1e-8) {
    stop("Unsupported clustering identity: ", x, ". Available: ", paste(keys, collapse = ", "))
  }
  names(threshold_columns)[[j]]
}

read_and_join_assignments <- function(root) {
  copy_file <- find_exact_file(root, "copy_variant_assignments.tsv")
  subtype_file <- find_exact_file(root, "subtype_copy_assignments.tsv")
  copy_assign <- read_tsv(copy_file)
  subtype <- read_tsv(subtype_file)

  require_columns(
    copy_assign,
    c("seq_id", "species", "source_array_key", unname(threshold_columns)),
    "copy_variant_assignments.tsv"
  )
  require_columns(
    subtype,
    c("seq_id", "species", "chrom", "array_key", "region_key", "start", "end", "haplotype_id"),
    "subtype_copy_assignments.tsv"
  )
  if (anyDuplicated(copy_assign$seq_id)) stop("Duplicate seq_id values in copy_variant_assignments.tsv")
  if (anyDuplicated(subtype$seq_id)) stop("Duplicate seq_id values in subtype_copy_assignments.tsv")

  idx <- match(subtype$seq_id, copy_assign$seq_id)
  if (anyNA(idx)) {
    warning(sum(is.na(idx)), " subtype rows lack copy-variant assignments and will be removed.")
    subtype <- subtype[!is.na(idx), , drop = FALSE]
    idx <- idx[!is.na(idx)]
  }
  species_match <- as.character(copy_assign$species[idx]) == as.character(subtype$species)
  if (any(!species_match, na.rm = TRUE)) stop("Species mismatch between the two assignment tables.")

  subtype$source_array_key <- as.character(copy_assign$source_array_key[idx])
  fallback <- !nonempty(subtype$source_array_key)
  subtype$source_array_key[fallback] <- as.character(subtype$array_key[fallback])
  for (column in unname(threshold_columns)) {
    subtype[[column]] <- as.character(copy_assign[[column]][idx])
  }
  subtype$start <- safe_numeric(subtype$start)
  subtype$end <- safe_numeric(subtype$end)
  subtype <- subtype[is.finite(subtype$start) & is.finite(subtype$end), , drop = FALSE]

  message("Input files:")
  message("  Copy variants:     ", copy_file)
  message("  Ordered copies:    ", subtype_file)
  list(data = subtype, copy_file = copy_file, subtype_file = subtype_file)
}

summarize_available_blocks <- function(x) {
  key <- interaction(
    x$species, x$chrom, x$region_key, x$source_array_key,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(x)), key)
  rows <- lapply(groups, function(ii) {
    d <- x[ii, , drop = FALSE]
    data.frame(
      species = as.character(d$species[[1L]]),
      chrom = as.character(d$chrom[[1L]]),
      region_key = as.character(d$region_key[[1L]]),
      source_array_key = as.character(d$source_array_key[[1L]]),
      n_copies = nrow(d),
      block_start = min(d$start, na.rm = TRUE),
      block_end = max(d$end, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$species, out$chrom, out$block_start, out$source_array_key), , drop = FALSE]
}

empty_selection <- function(n = 0L) {
  data.frame(
    species = rep(NA_character_, n),
    chrom = rep(NA_character_, n),
    region_key = rep(NA_character_, n),
    source_array_key = rep(NA_character_, n),
    expected_period = rep(NA_integer_, n),
    copy_start = rep(NA_integer_, n),
    copy_end = rep(NA_integer_, n),
    label = rep(NA_character_, n),
    enabled = rep(TRUE, n),
    selection_source = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

normalize_selection <- function(x, source_label) {
  aliases <- list(
    source_array_key = c("array_key", "block", "block_key"),
    expected_period = c("period", "hor_period"),
    copy_start = c("start_copy"),
    copy_end = c("end_copy")
  )
  for (canonical in names(aliases)) {
    if (!canonical %in% names(x)) {
      hit <- intersect(aliases[[canonical]], names(x))
      if (length(hit) > 0L) names(x)[match(hit[[1L]], names(x))] <- canonical
    }
  }
  if (!"species" %in% names(x)) stop(source_label, " lacks required column: species")
  for (column in setdiff(names(empty_selection()), c("species", "selection_source"))) {
    if (!column %in% names(x)) x[[column]] <- NA
  }
  x$species <- trimws(as.character(x$species))
  x$chrom <- trimws(as.character(x$chrom))
  x$region_key <- trimws(as.character(x$region_key))
  x$source_array_key <- trimws(as.character(x$source_array_key))
  x$label <- trimws(as.character(x$label))
  x$expected_period <- suppressWarnings(as.integer(as.character(x$expected_period)))
  x$copy_start <- suppressWarnings(as.integer(as.character(x$copy_start)))
  x$copy_end <- suppressWarnings(as.integer(as.character(x$copy_end)))
  enabled_text <- tolower(trimws(as.character(x$enabled)))
  x$enabled <- is.na(x$enabled) | enabled_text %in% c("", "1", "true", "yes", "y", "on")
  x$selection_source <- source_label
  x <- x[x$enabled & nonempty(x$species), names(empty_selection()), drop = FALSE]
  rownames(x) <- NULL
  x
}

preset_selection <- function(name) {
  name <- tolower(name)
  if (name == "none") return(empty_selection())
  if (!name %in% c("strong", "extended", "all")) {
    stop("Unknown --preset: ", name, ". Use strong, extended, all, or none.")
  }
  x <- data.frame(
    species = c(
      rep("FCA", 6),
      rep("PPA", 8),
      "PTI"
    ),
    region_key = c(
      "ChrF1__region001", "ChrF1__region003", "ChrF1__region005",
      "ChrF1__region006", "ChrF1__region006", "ChrF1__region010",
      "ChrE1__region004", "ChrE1__region004",
      "ChrE1__region005", "ChrE1__region005", "ChrE1__region005", "ChrE1__region005",
      "ChrE1__region006", "ChrE1__region006",
      "ChrA1__region001"
    ),
    source_array_key = c(
      "ChrF1__28235", "ChrF1__28244", "ChrF1__28252",
      "ChrF1__28258", "ChrF1__28259", "ChrF1__28271",
      "ChrE1__21788", "ChrE1__21797",
      "ChrE1__21808", "ChrE1__21814", "ChrE1__21819", "ChrE1__21824",
      "ChrE1__21830", "ChrE1__21831",
      "ChrA1__694"
    ),
    expected_period = c(8, 8, 2, 8, 21, 8, 62, 65, 65, 65, 62, 65, 66, 66, 81),
    stringsAsFactors = FALSE
  )
  if (name == "strong") {
    x <- x[!x$source_array_key %in% c("ChrE1__21797", "ChrE1__21819", "ChrE1__21824"), , drop = FALSE]
  }
  normalize_selection(x, paste0("preset_", ifelse(name == "all", "extended", name)))
}

parse_inline_blocks <- function(specs) {
  if (length(specs) == 0L) return(empty_selection())
  rows <- lapply(specs, function(spec) {
    fields <- trimws(strsplit(spec, ",", fixed = TRUE)[[1L]])
    if (length(fields) < 2L || length(fields) > 6L) {
      stop("Invalid --block specification: ", spec)
    }
    fields <- c(fields, rep(NA_character_, 6L - length(fields)))
    period <- if (!is.na(fields[[3L]]) && tolower(fields[[3L]]) != "auto" && nzchar(fields[[3L]])) {
      suppressWarnings(as.integer(fields[[3L]]))
    } else {
      NA_integer_
    }
    data.frame(
      species = fields[[1L]],
      source_array_key = fields[[2L]],
      expected_period = period,
      copy_start = suppressWarnings(as.integer(fields[[4L]])),
      copy_end = suppressWarnings(as.integer(fields[[5L]])),
      label = fields[[6L]],
      stringsAsFactors = FALSE
    )
  })
  normalize_selection(do.call(rbind, rows), "inline_block")
}

read_config_selection <- function(path) {
  if (is.null(path)) return(empty_selection())
  if (!file.exists(path)) stop("Block config does not exist: ", path)
  normalize_selection(read_tsv(path), paste0("config:", basename(path)))
}

expand_selection <- function(selection, available, min_block_copies) {
  if (nrow(selection) == 0L) stop("No blocks selected.")
  expanded <- list()
  out_i <- 0L
  for (i in seq_len(nrow(selection))) {
    s <- selection[i, , drop = FALSE]
    hit <- available$species == s$species
    if (nonempty(s$chrom)) hit <- hit & available$chrom == s$chrom
    if (nonempty(s$region_key)) hit <- hit & available$region_key == s$region_key
    if (nonempty(s$source_array_key)) hit <- hit & available$source_array_key == s$source_array_key
    candidates <- available[hit, , drop = FALSE]
    if (!nonempty(s$source_array_key)) {
      candidates <- candidates[candidates$n_copies >= min_block_copies, , drop = FALSE]
    }
    if (nrow(candidates) == 0L) {
      stop(
        "Selection did not match an available block: species=", s$species,
        ", region_key=", s$region_key, ", source_array_key=", s$source_array_key
      )
    }
    candidates <- candidates[order(candidates$block_start, candidates$source_array_key), , drop = FALSE]
    for (j in seq_len(nrow(candidates))) {
      cnd <- candidates[j, , drop = FALSE]
      start_copy <- if (is.finite(s$copy_start) && s$copy_start >= 1L) s$copy_start else 1L
      end_copy <- if (is.finite(s$copy_end) && s$copy_end >= start_copy) s$copy_end else cnd$n_copies
      if (start_copy > cnd$n_copies) {
        stop("copy_start exceeds block length for ", cnd$species, " ", cnd$source_array_key)
      }
      if (end_copy > cnd$n_copies) {
        warning("Clipping copy_end from ", end_copy, " to ", cnd$n_copies, " for ", cnd$source_array_key)
        end_copy <- cnd$n_copies
      }
      label <- if (nonempty(s$label)) s$label else paste(cnd$species, cnd$chrom, cnd$region_key, cnd$source_array_key, sep = "_")
      if (nrow(candidates) > 1L && nonempty(s$label)) label <- paste(label, cnd$source_array_key, sep = "_")
      out_i <- out_i + 1L
      expanded[[out_i]] <- data.frame(
        species = cnd$species,
        chrom = cnd$chrom,
        region_key = cnd$region_key,
        source_array_key = cnd$source_array_key,
        expected_period = if (is.finite(s$expected_period) && s$expected_period >= 1L) s$expected_period else NA_integer_,
        copy_start = as.integer(start_copy),
        copy_end = as.integer(end_copy),
        label = label,
        selection_source = s$selection_source,
        total_block_copies = cnd$n_copies,
        block_start = cnd$block_start,
        block_end = cnd$block_end,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, expanded)
  dedup_key <- paste(
    out$species, out$chrom, out$region_key, out$source_array_key,
    out$copy_start, out$copy_end, sep = "|"
  )
  if (anyDuplicated(dedup_key)) {
    warning("Duplicate block selections were collapsed; the first occurrence was retained.")
    out <- out[!duplicated(dedup_key), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0L) stop("No FASTA records found in: ", path)
  ends <- c(headers[-1L] - 1L, length(lines))
  ids <- sub("^>([^[:space:]]+).*$", "\\1", lines[headers])
  seqs <- vapply(seq_along(headers), function(i) {
    if (ends[[i]] <= headers[[i]]) return("")
    toupper(paste0(lines[(headers[[i]] + 1L):ends[[i]]], collapse = ""))
  }, character(1))
  if (anyDuplicated(ids)) stop("Duplicate FASTA identifiers in: ", path)
  setNames(seqs, ids)
}

exact_identity_matrix <- function(haplotype_ids, fasta_sequences) {
  fasta_index <- match(haplotype_ids, names(fasta_sequences))
  found <- !is.na(fasta_index)
  if (sum(found) < 2L) return(NULL)
  seqs <- unname(fasta_sequences[fasta_index[found]])
  lengths <- nchar(seqs)
  if (length(unique(lengths)) != 1L || unique(lengths) < 1L) {
    warning("Selected FASTA sequences are not a same-length alignment; using cluster agreement instead.")
    return(NULL)
  }
  chars <- do.call(rbind, strsplit(seqs, "", fixed = TRUE))
  bases <- c("A", "C", "G", "T")
  numerator <- matrix(0, nrow(chars), nrow(chars))
  for (base in bases) {
    indicator <- (chars == base) * 1
    numerator <- numerator + tcrossprod(indicator)
  }
  valid <- matrix(chars %in% bases, nrow = nrow(chars), ncol = ncol(chars)) * 1
  denominator <- tcrossprod(valid)
  submatrix <- numerator / denominator
  submatrix[denominator == 0] <- NA_real_
  out <- matrix(NA_real_, length(haplotype_ids), length(haplotype_ids))
  out[found, found] <- submatrix
  list(matrix = out, coverage = mean(found), aligned_columns = ncol(chars))
}

cluster_agreement_matrix <- function(d) {
  baseline <- min(thresholds) - 0.05
  out <- matrix(baseline, nrow(d), nrow(d))
  for (threshold in thresholds) {
    key <- threshold_key(threshold)
    column <- threshold_columns[[key]]
    labels <- as.character(d[[column]])
    equal <- outer(labels, labels, "==")
    equal[is.na(equal)] <- FALSE
    out[equal] <- threshold
  }
  diag(out) <- max(thresholds)
  list(matrix = out, coverage = 1, aligned_columns = NA_integer_, baseline = baseline)
}

lag_profile_one <- function(
    labels,
    max_lag,
    threshold,
    min_lag_pairs = 1L,
    min_pair_fraction = 0) {
  labels <- as.character(labels)
  labels[is.na(labels) | !nzchar(labels) |
           toupper(labels) %in% c("NA", "NAN", "NONE", "NULL", "UNASSIGNED")] <- NA_character_
  n <- length(labels)
  valid_labels <- labels[!is.na(labels)]
  counts <- table(valid_labels)
  random_expectation <- if (length(valid_labels) > 1L) {
    sum(counts * (counts - 1)) / (length(valid_labels) * (length(valid_labels) - 1))
  } else {
    NA_real_
  }
  lags <- seq_len(min(max_lag, n - 1L))
  minimum_pairs <- max(as.integer(min_lag_pairs), as.integer(ceiling(min_pair_fraction * n)))
  pair_counts <- vapply(lags, function(lag) {
    left <- labels[seq_len(n - lag)]
    right <- labels[seq_len(n - lag) + lag]
    keep <- !is.na(left) & !is.na(right) & nzchar(left) & nzchar(right)
    sum(keep)
  }, integer(1))
  observed <- vapply(seq_along(lags), function(i) {
    lag <- lags[[i]]
    if (pair_counts[[i]] < minimum_pairs) return(NA_real_)
    left <- labels[seq_len(n - lag)]
    right <- labels[seq_len(n - lag) + lag]
    keep <- !is.na(left) & !is.na(right)
    mean(left[keep] == right[keep])
  }, numeric(1))
  data.frame(
    lag = lags,
    n_pairs = pair_counts,
    minimum_pairs = rep(minimum_pairs, length(lags)),
    meets_minimum_pairs = pair_counts >= minimum_pairs,
    threshold = threshold,
    threshold_label = threshold_labels[[threshold_key(threshold)]],
    same_cluster_fraction = observed,
    random_order_expectation = random_expectation,
    periodicity_excess = observed - random_expectation,
    stringsAsFactors = FALSE
  )
}

label_matrix_from_block <- function(d) {
  out <- as.matrix(d[, unname(threshold_columns), drop = FALSE])
  storage.mode(out) <- "character"
  missing <- is.na(out) | !nzchar(out) |
    toupper(out) %in% c("NA", "NAN", "NONE", "NULL", "UNASSIGNED")
  out[missing] <- NA_character_
  out
}

random_expectation_by_threshold <- function(label_matrix) {
  vapply(seq_len(ncol(label_matrix)), function(j) {
    values <- label_matrix[, j]
    values <- values[!is.na(values)]
    n <- length(values)
    if (n < 2L) return(NA_real_)
    counts <- table(values)
    sum(counts * (counts - 1)) / (n * (n - 1))
  }, numeric(1))
}

lag_profile_matrix <- function(label_matrix, lags, random_expectation = NULL) {
  n <- nrow(label_matrix)
  m <- ncol(label_matrix)
  if (is.null(random_expectation)) {
    random_expectation <- random_expectation_by_threshold(label_matrix)
  }
  observed <- matrix(NA_real_, nrow = length(lags), ncol = m)
  pair_counts <- matrix(0L, nrow = length(lags), ncol = m)
  for (i in seq_along(lags)) {
    lag <- lags[[i]]
    left <- label_matrix[seq_len(n - lag), , drop = FALSE]
    right <- label_matrix[seq_len(n - lag) + lag, , drop = FALSE]
    valid <- !is.na(left) & !is.na(right)
    pair_counts[i, ] <- colSums(valid)
    equal <- left == right
    equal[!valid] <- FALSE
    sums <- colSums(equal)
    keep <- pair_counts[i, ] > 0L
    observed[i, keep] <- sums[keep] / pair_counts[i, keep]
  }
  excess <- sweep(observed, 2L, random_expectation, FUN = "-")
  robust <- apply(excess, 1L, function(x) {
    if (!any(is.finite(x))) return(NA_real_)
    stats::median(x[is.finite(x)])
  })
  list(
    observed = observed,
    excess = excess,
    robust = robust,
    pair_counts = pair_counts,
    random_expectation = random_expectation
  )
}

enforce_minimum_pairs <- function(profile, minimum_pairs) {
  insufficient <- profile$pair_counts < minimum_pairs
  profile$observed[insufficient] <- NA_real_
  profile$excess[insufficient] <- NA_real_
  profile$robust <- apply(profile$excess, 1L, function(x) {
    if (!any(is.finite(x))) return(NA_real_)
    stats::median(x[is.finite(x)])
  })
  profile
}

select_fundamental_lag <- function(
    lags,
    robust_scores,
    absolute_tolerance = 0.02,
    relative_tolerance = 0.05) {
  finite_index <- which(is.finite(robust_scores))
  if (!length(finite_index)) stop("No finite lag scores are available.")
  peak_index <- finite_index[[which.max(robust_scores[finite_index])]]
  peak_lag <- as.integer(lags[[peak_index]])
  peak_score <- robust_scores[[peak_index]]
  tolerance <- max(absolute_tolerance, relative_tolerance * abs(peak_score))
  near_peak <- finite_index[robust_scores[finite_index] >= peak_score - tolerance]
  # A putative fundamental must belong to the same harmonic family as the
  # strongest peak. This prevents an unrelated low-lag fluctuation from being
  # preferred merely because it is numerically smaller and nearly as high.
  harmonic <- near_peak[peak_lag %% as.integer(lags[near_peak]) == 0L]
  candidate_index <- harmonic[[which.min(as.integer(lags[harmonic]))]]
  list(
    peak_index = as.integer(peak_index),
    peak_lag = peak_lag,
    peak_score = peak_score,
    tolerance = tolerance,
    near_peak_indices = as.integer(near_peak),
    harmonic_near_peak_indices = as.integer(harmonic),
    candidate_index = as.integer(candidate_index),
    candidate_lag = as.integer(lags[[candidate_index]])
  )
}

scan_label_matrix <- function(
    label_matrix,
    min_period = 2L,
    max_period = 120L,
    min_units = 2,
    min_lag_pairs = 8L,
    min_pair_fraction = 0.25,
    support_excess = 0.10,
    fundamental_abs_tolerance = 0.02,
    fundamental_rel_tolerance = 0.05) {
  label_matrix <- as.matrix(label_matrix)
  storage.mode(label_matrix) <- "character"
  n <- nrow(label_matrix)
  minimum_pairs <- max(as.integer(min_lag_pairs), as.integer(ceiling(min_pair_fraction * n)))
  tested_max <- min(as.integer(max_period), as.integer(floor(n / min_units)), n - minimum_pairs)
  if (n < 3L || tested_max < min_period) {
    return(list(test_status = "NOT_TESTED", not_tested_reason = "INSUFFICIENT_PERIOD_RANGE"))
  }
  lags <- seq.int(as.integer(min_period), as.integer(tested_max))
  profile <- enforce_minimum_pairs(
    lag_profile_matrix(label_matrix, lags),
    minimum_pairs
  )
  eligible <- is.finite(profile$robust)
  if (!any(eligible)) {
    return(list(test_status = "NOT_TESTED", not_tested_reason = "NO_VALID_LAG"))
  }
  eligible_scores <- profile$robust
  eligible_scores[!eligible] <- NA_real_
  selected_lag <- select_fundamental_lag(
    lags,
    eligible_scores,
    absolute_tolerance = fundamental_abs_tolerance,
    relative_tolerance = fundamental_rel_tolerance
  )
  peak_index <- selected_lag$peak_index
  peak_score <- selected_lag$peak_score
  candidate_index <- selected_lag$candidate_index
  period <- selected_lag$candidate_lag
  candidate_excess <- profile$excess[candidate_index, ]
  candidate_observed <- profile$observed[candidate_index, ]
  support <- as.integer(sum(candidate_excess >= support_excess, na.rm = TRUE))
  multiples <- which(lags %% period == 0L)
  multiple_cutoff <- max(support_excess, 0.5 * profile$robust[[candidate_index]])
  multiple_support <- as.integer(sum(profile$robust[multiples] >= multiple_cutoff, na.rm = TRUE))
  neighbour_index <- which(abs(lags - period) == 1L)
  neighbour_max <- if (length(neighbour_index) && any(is.finite(profile$robust[neighbour_index]))) {
    max(profile$robust[neighbour_index], na.rm = TRUE)
  } else {
    NA_real_
  }
  list(
    test_status = "TESTED",
    not_tested_reason = "",
    n_copies = n,
    lags = lags,
    minimum_pairs = minimum_pairs,
    observed = profile$observed,
    excess = profile$excess,
    robust = profile$robust,
    pair_counts = profile$pair_counts,
    random_expectation = profile$random_expectation,
    peak_index = peak_index,
    peak_score = peak_score,
    candidate_index = candidate_index,
    candidate_period = period,
    candidate_robust_excess = profile$robust[[candidate_index]],
    candidate_threshold_support = support,
    candidate_observed = candidate_observed,
    candidate_excess = candidate_excess,
    estimated_units = n / period,
    multiple_support = multiple_support,
    local_contrast = profile$robust[[candidate_index]] - neighbour_max
  )
}

max_stat_permutation <- function(label_matrix, scan, permutations, seed) {
  if (permutations <= 0L || scan$test_status != "TESTED") {
    return(c(candidate_point_p = NA_real_, block_maxT_p = NA_real_))
  }
  set.seed(seed)
  null_max <- numeric(permutations)
  null_candidate <- numeric(permutations)
  for (i in seq_len(permutations)) {
    shuffled <- label_matrix[sample.int(nrow(label_matrix), replace = FALSE), , drop = FALSE]
    profile <- enforce_minimum_pairs(
      lag_profile_matrix(shuffled, scan$lags, scan$random_expectation),
      scan$minimum_pairs
    )
    null_max[[i]] <- max(profile$robust, na.rm = TRUE)
    null_candidate[[i]] <- profile$robust[[scan$candidate_index]]
  }
  epsilon <- sqrt(.Machine$double.eps)
  c(
    candidate_point_p = (1 + sum(null_candidate >= scan$candidate_robust_excess - epsilon, na.rm = TRUE)) /
      (permutations + 1),
    block_maxT_p = (1 + sum(null_max >= scan$peak_score - epsilon, na.rm = TRUE)) /
      (permutations + 1)
  )
}

classify_hor_call <- function(
    q_value,
    robust_excess,
    support,
    permutations_enabled,
    high_fdr = 0.05,
    high_excess = 0.20,
    high_support = 3L,
    moderate_fdr = 0.10,
    moderate_excess = 0.10,
    moderate_support = 2L) {
  effect_high <- is.finite(robust_excess) && robust_excess >= high_excess && support >= high_support
  effect_moderate <- is.finite(robust_excess) && robust_excess >= moderate_excess && support >= moderate_support
  if (!permutations_enabled) {
    if (effect_high) return("EFFECT_ONLY_HIGH")
    if (effect_moderate) return("EFFECT_ONLY_MODERATE")
    return("NOT_DETECTED")
  }
  if (effect_high && is.finite(q_value) && q_value <= high_fdr) return("HIGH")
  if (effect_moderate && is.finite(q_value) && q_value <= moderate_fdr) return("MODERATE")
  "NOT_DETECTED"
}

lag_permutation <- function(labels, lag, n_permutations, min_valid_pairs = 1L) {
  labels <- as.character(labels)
  labels[is.na(labels) | !nzchar(labels) |
           toupper(labels) %in% c("NA", "NAN", "NONE", "NULL", "UNASSIGNED")] <- NA_character_
  n <- length(labels)
  if (!is.finite(lag) || lag < 1L || lag >= n) {
    return(c(n_valid_pairs = 0, observed = NA, null_mean = NA,
             null_low = NA, null_high = NA, p_enrichment = NA))
  }
  score <- function(x) {
    left <- x[seq_len(n - lag)]
    right <- x[seq_len(n - lag) + lag]
    valid <- !is.na(left) & !is.na(right)
    if (sum(valid) < min_valid_pairs) return(NA_real_)
    mean(left[valid] == right[valid])
  }
  observed_left <- labels[seq_len(n - lag)]
  observed_right <- labels[seq_len(n - lag) + lag]
  n_valid_pairs <- sum(!is.na(observed_left) & !is.na(observed_right))
  observed <- score(labels)
  if (n_permutations <= 0L) {
    return(c(n_valid_pairs = n_valid_pairs, observed = observed, null_mean = NA,
             null_low = NA, null_high = NA, p_enrichment = NA))
  }
  null <- replicate(n_permutations, score(sample(labels, replace = FALSE)))
  valid_null <- null[is.finite(null)]
  if (!is.finite(observed) || !length(valid_null)) {
    return(c(n_valid_pairs = n_valid_pairs, observed = observed, null_mean = NA,
             null_low = NA, null_high = NA, p_enrichment = NA))
  }
  c(
    n_valid_pairs = n_valid_pairs,
    observed = observed,
    null_mean = mean(valid_null),
    null_low = unname(stats::quantile(valid_null, 0.025)),
    null_high = unname(stats::quantile(valid_null, 0.975)),
    p_enrichment = (1 + sum(valid_null >= observed)) / (length(valid_null) + 1)
  )
}

format_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(format(x, scientific = TRUE, digits = 2))
  sprintf("%.3f", x)
}

format_percent <- function(x, digits = 1L) {
  ifelse(is.finite(x), paste0(formatC(100 * x, format = "f", digits = digits), "%"), "NA")
}

scan_all_blocks <- function(copies, available, opt, satellite_family) {
  summaries <- vector("list", nrow(available))
  lag_tables <- vector("list", nrow(available))
  message("Automatic HOR scan: ", nrow(available), " source-array blocks")
  for (i in seq_len(nrow(available))) {
    block_info <- available[i, , drop = FALSE]
    d <- copies[
      copies$species == block_info$species &
        copies$chrom == block_info$chrom &
        copies$region_key == block_info$region_key &
        copies$source_array_key == block_info$source_array_key,
      , drop = FALSE
    ]
    d <- d[order(d$start, d$end, d$seq_id), , drop = FALSE]
    n <- nrow(d)
    if (i == 1L || i %% 25L == 0L || i == nrow(available)) {
      message("  Scanning block ", i, "/", nrow(available), ": ",
              block_info$species, " ", block_info$source_array_key, " (n=", n, ")")
    }
    copy_length <- if ("phased_length" %in% names(d)) safe_numeric(d$phased_length) else d$end - d$start + 1
    copy_length <- copy_length[is.finite(copy_length) & copy_length > 0]
    median_copy_length <- if (length(copy_length)) stats::median(copy_length) else NA_real_
    base <- data.frame(
      satellite = satellite_family,
      species = as.character(block_info$species),
      chrom = as.character(block_info$chrom),
      region_key = as.character(block_info$region_key),
      source_array_key = as.character(block_info$source_array_key),
      block_start = block_info$block_start,
      block_end = block_info$block_end,
      n_copies = n,
      median_copy_length = median_copy_length,
      test_status = "NOT_TESTED",
      not_tested_reason = "",
      tested_min_period = opt$min_period,
      tested_max_period = NA_integer_,
      minimum_lag_pairs = NA_integer_,
      candidate_period = NA_integer_,
      candidate_period_bp = NA_real_,
      peak_lag = NA_integer_,
      estimated_units = NA_real_,
      complete_units = NA_integer_,
      residual_copies = NA_integer_,
      repeat_depth_flag = NA_character_,
      candidate_robust_excess = NA_real_,
      peak_robust_excess = NA_real_,
      candidate_threshold_support = NA_integer_,
      n_cluster_thresholds = length(thresholds),
      candidate_median_same_fraction = NA_real_,
      candidate_median_random_expectation = NA_real_,
      primary_identity = opt$main_identity,
      primary_cluster_count = NA_integer_,
      primary_dominant_cluster_fraction = NA_real_,
      primary_adjacent_same_fraction = NA_real_,
      primary_candidate_same_fraction = NA_real_,
      primary_candidate_random_expectation = NA_real_,
      primary_candidate_periodicity_excess = NA_real_,
      candidate_min_threshold_excess = NA_real_,
      candidate_max_threshold_excess = NA_real_,
      harmonic_multiple_support = NA_integer_,
      testable_multiple_count = NA_integer_,
      local_peak_contrast = NA_real_,
      permutation_screen_status = "NOT_SCREENED",
      permutation_elapsed_seconds = 0,
      candidate_point_p = if (opt$permutations > 0L) 1 else NA_real_,
      block_maxT_p = if (opt$permutations > 0L) 1 else NA_real_,
      block_fdr_q = NA_real_,
      hor_call = "NOT_TESTED",
      stringsAsFactors = FALSE
    )
    if (n < opt$min_copies) {
      base$not_tested_reason <- paste0("N_COPIES_LT_", opt$min_copies)
      summaries[[i]] <- base
      next
    }
    label_matrix <- label_matrix_from_block(d)
    scan <- scan_label_matrix(
      label_matrix,
      min_period = opt$min_period,
      max_period = opt$max_period,
      min_units = opt$min_units,
      min_lag_pairs = opt$min_lag_pairs,
      min_pair_fraction = opt$min_pair_fraction,
      support_excess = opt$support_excess,
      fundamental_abs_tolerance = opt$fundamental_abs_tolerance,
      fundamental_rel_tolerance = opt$fundamental_rel_tolerance
    )
    if (scan$test_status != "TESTED") {
      base$not_tested_reason <- scan$not_tested_reason
      summaries[[i]] <- base
      next
    }
    main_index <- which.min(abs(thresholds - opt$main_identity))
    base$test_status <- "TESTED"
    base$tested_max_period <- max(scan$lags)
    base$minimum_lag_pairs <- scan$minimum_pairs
    base$candidate_period <- scan$candidate_period
    base$candidate_period_bp <- scan$candidate_period * median_copy_length
    base$peak_lag <- scan$lags[[scan$peak_index]]
    base$estimated_units <- scan$estimated_units
    base$complete_units <- n %/% scan$candidate_period
    base$residual_copies <- n %% scan$candidate_period
    base$repeat_depth_flag <- if (scan$estimated_units < 3) "LIMITED_LT_3_UNITS" else "AT_LEAST_3_UNITS"
    base$candidate_robust_excess <- scan$candidate_robust_excess
    base$peak_robust_excess <- scan$peak_score
    base$candidate_threshold_support <- scan$candidate_threshold_support
    base$candidate_median_same_fraction <- stats::median(scan$candidate_observed, na.rm = TRUE)
    base$candidate_median_random_expectation <- stats::median(scan$random_expectation, na.rm = TRUE)
    primary_labels <- label_matrix[, main_index]
    primary_labels <- primary_labels[!is.na(primary_labels)]
    primary_counts <- table(primary_labels)
    base$primary_cluster_count <- length(primary_counts)
    base$primary_dominant_cluster_fraction <- if (length(primary_counts)) max(primary_counts) / sum(primary_counts) else NA_real_
    if (n >= 2L) {
      left_primary <- label_matrix[seq_len(n - 1L), main_index]
      right_primary <- label_matrix[seq_len(n - 1L) + 1L, main_index]
      valid_primary <- !is.na(left_primary) & !is.na(right_primary)
      base$primary_adjacent_same_fraction <- if (any(valid_primary)) mean(left_primary[valid_primary] == right_primary[valid_primary]) else NA_real_
    }
    base$primary_candidate_same_fraction <- scan$candidate_observed[[main_index]]
    base$primary_candidate_random_expectation <- scan$random_expectation[[main_index]]
    base$primary_candidate_periodicity_excess <- scan$candidate_excess[[main_index]]
    base$candidate_min_threshold_excess <- min(scan$candidate_excess, na.rm = TRUE)
    base$candidate_max_threshold_excess <- max(scan$candidate_excess, na.rm = TRUE)
    base$harmonic_multiple_support <- scan$multiple_support
    base$testable_multiple_count <- sum(scan$lags %% scan$candidate_period == 0L)
    base$local_peak_contrast <- scan$local_contrast

    screened <- scan$candidate_robust_excess >= opt$screen_excess &&
      scan$candidate_threshold_support >= opt$screen_support
    if (screened && opt$permutations > 0L) {
      message("    ", opt$permutations, " max-lag permutations: ",
              block_info$species, " ", block_info$source_array_key)
      permutation_timing <- system.time({
        permutation <- max_stat_permutation(
          label_matrix, scan, opt$permutations,
          opt$seed + i * 1009L
        )
      })
      base$permutation_screen_status <- "SCREENED_IN"
      base$permutation_elapsed_seconds <- unname(permutation_timing[["elapsed"]])
      base$candidate_point_p <- permutation[["candidate_point_p"]]
      base$block_maxT_p <- permutation[["block_maxT_p"]]
    } else if (opt$permutations <= 0L) {
      base$permutation_screen_status <- "PERMUTATIONS_DISABLED"
    } else {
      base$permutation_screen_status <- "SCREENED_OUT_EFFECT_OR_STABILITY"
    }
    summaries[[i]] <- base

    lag_rows <- vector("list", length(thresholds))
    support_by_lag <- rowSums(scan$excess >= opt$support_excess, na.rm = TRUE)
    median_observed <- apply(scan$observed, 1L, stats::median, na.rm = TRUE)
    for (j in seq_along(thresholds)) {
      lag_rows[[j]] <- data.frame(
        satellite = satellite_family,
        species = block_info$species,
        chrom = block_info$chrom,
        region_key = block_info$region_key,
        source_array_key = block_info$source_array_key,
        lag = scan$lags,
        threshold = thresholds[[j]],
        threshold_label = threshold_labels[[threshold_key(thresholds[[j]])]],
        n_pairs = scan$pair_counts[, j],
        same_cluster_fraction = scan$observed[, j],
        random_order_expectation = scan$random_expectation[[j]],
        periodicity_excess = scan$excess[, j],
        robust_excess_across_thresholds = scan$robust,
        median_same_fraction_across_thresholds = median_observed,
        threshold_support_count = support_by_lag,
        global_peak_lag = scan$lags == scan$lags[[scan$peak_index]],
        selected_candidate_lag = scan$lags == scan$candidate_period,
        stringsAsFactors = FALSE
      )
    }
    lag_tables[[i]] <- do.call(rbind, lag_rows)
  }
  summary <- do.call(rbind, summaries)
  tested <- summary$test_status == "TESTED"
  if (opt$permutations > 0L && any(tested)) {
    summary$block_fdr_q[tested] <- stats::p.adjust(summary$block_maxT_p[tested], method = "BH")
  }
  summary$hor_call[tested] <- vapply(which(tested), function(i) {
    classify_hor_call(
      summary$block_fdr_q[[i]],
      summary$candidate_robust_excess[[i]],
      summary$candidate_threshold_support[[i]],
      permutations_enabled = opt$permutations > 0L,
      high_fdr = opt$high_fdr,
      high_excess = opt$high_excess,
      high_support = opt$high_support,
      moderate_fdr = opt$moderate_fdr,
      moderate_excess = opt$moderate_excess,
      moderate_support = opt$moderate_support
    )
  }, character(1))
  call_order <- c(HIGH = 1L, EFFECT_ONLY_HIGH = 1L, MODERATE = 2L,
                  EFFECT_ONLY_MODERATE = 2L, NOT_DETECTED = 3L, NOT_TESTED = 4L)
  summary$.call_order <- unname(call_order[summary$hor_call])
  summary <- summary[order(summary$.call_order, -summary$candidate_robust_excess,
                           summary$species, summary$chrom, summary$block_start), , drop = FALSE]
  summary$.call_order <- NULL
  candidate_mask <- summary$hor_call %in% opt$candidate_calls
  summary$candidate_rank <- NA_integer_
  summary$candidate_rank[candidate_mask] <- seq_len(sum(candidate_mask))
  rownames(summary) <- NULL
  all_lags <- if (any(vapply(lag_tables, Negate(is.null), logical(1)))) {
    do.call(rbind, lag_tables[!vapply(lag_tables, is.null, logical(1))])
  } else {
    data.frame(
      satellite = character(), species = character(), chrom = character(),
      region_key = character(), source_array_key = character(), lag = integer(),
      threshold = numeric(), threshold_label = character(), n_pairs = integer(),
      same_cluster_fraction = numeric(), random_order_expectation = numeric(),
      periodicity_excess = numeric(), robust_excess_across_thresholds = numeric(),
      median_same_fraction_across_thresholds = numeric(), threshold_support_count = integer(),
      global_peak_lag = logical(),
      selected_candidate_lag = logical(), stringsAsFactors = FALSE
    )
  }
  list(summary = summary, candidates = summary[candidate_mask, , drop = FALSE], all_lags = all_lags)
}

make_region_summary <- function(block_summary) {
  group_key <- interaction(block_summary$species, block_summary$chrom,
                           block_summary$region_key, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(block_summary)), group_key)
  out <- lapply(groups, function(ii) {
    d <- block_summary[ii, , drop = FALSE]
    tested <- d[d$test_status == "TESTED", , drop = FALSE]
    best <- if (nrow(tested)) tested[which.max(tested$candidate_robust_excess), , drop = FALSE] else d[1L, , drop = FALSE]
    candidate <- d[d$hor_call %in% c("HIGH", "MODERATE", "EFFECT_ONLY_HIGH", "EFFECT_ONLY_MODERATE"), , drop = FALSE]
    region_call <- if (any(grepl("HIGH", candidate$hor_call, fixed = TRUE))) "HIGH" else if (nrow(candidate)) "MODERATE" else if (nrow(tested)) "NOT_DETECTED" else "NOT_TESTED"
    data.frame(
      species = d$species[[1L]], chrom = d$chrom[[1L]], region_key = d$region_key[[1L]],
      region_call = region_call, n_blocks = nrow(d), tested_block_count = nrow(tested),
      candidate_block_count = nrow(candidate),
      best_source_array_key = best$source_array_key[[1L]],
      best_candidate_period = best$candidate_period[[1L]],
      best_robust_excess = best$candidate_robust_excess[[1L]],
      best_block_fdr_q = best$block_fdr_q[[1L]],
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

make_evidence_config <- function(candidates) {
  if (!nrow(candidates)) return(empty_selection())
  data.frame(
    species = candidates$species,
    chrom = candidates$chrom,
    region_key = candidates$region_key,
    source_array_key = candidates$source_array_key,
    expected_period = as.integer(candidates$candidate_period),
    copy_start = 1L,
    copy_end = as.integer(candidates$n_copies),
    label = paste(candidates$species, candidates$region_key,
                  candidates$source_array_key, paste0("p", candidates$candidate_period), sep = "_"),
    enabled = TRUE,
    hor_call = candidates$hor_call,
    candidate_robust_excess = candidates$candidate_robust_excess,
    block_fdr_q = candidates$block_fdr_q,
    estimated_units = candidates$estimated_units,
    stringsAsFactors = FALSE
  )
}

make_scan_metadata <- function(opt, satellite_family, scan_result) {
  values <- list(
    script = SCRIPT_NAME,
    script_version = SCRIPT_VERSION,
    satellite = satellite_family,
    total_blocks = nrow(scan_result$summary),
    tested_blocks = sum(scan_result$summary$test_status == "TESTED"),
    candidate_blocks = nrow(scan_result$candidates),
    screened_permutation_blocks = sum(scan_result$summary$permutation_screen_status == "SCREENED_IN"),
    permutation_elapsed_seconds = sum(scan_result$summary$permutation_elapsed_seconds, na.rm = TRUE),
    min_copies = opt$min_copies,
    min_period = opt$min_period,
    max_period = opt$max_period,
    min_units = opt$min_units,
    min_lag_pairs = opt$min_lag_pairs,
    min_pair_fraction = opt$min_pair_fraction,
    fundamental_abs_tolerance = opt$fundamental_abs_tolerance,
    fundamental_rel_tolerance = opt$fundamental_rel_tolerance,
    cluster_thresholds = paste(thresholds, collapse = ","),
    support_excess = opt$support_excess,
    screen_excess = opt$screen_excess,
    screen_support = opt$screen_support,
    permutations = opt$permutations,
    seed = opt$seed,
    high_fdr = opt$high_fdr,
    high_excess = opt$high_excess,
    high_support = opt$high_support,
    moderate_fdr = opt$moderate_fdr,
    moderate_excess = opt$moderate_excess,
    moderate_support = opt$moderate_support,
    fundamental_selection = "smallest_near_peak_divisor_of_global_peak",
    pair_filter = "per_threshold_minimum_pairs_before_robust_median",
    candidate_calls = paste(opt$candidate_calls, collapse = ",")
  )
  data.frame(
    parameter = names(values),
    value = vapply(values, as.character, character(1)),
    stringsAsFactors = FALSE
  )
}

make_empty_plot <- function(title, message, base_theme) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title) +
    base_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}

make_scan_overview_plot <- function(scan_summary, satellite_family, opt, base_theme) {
  d <- scan_summary[scan_summary$test_status == "TESTED", , drop = FALSE]
  if (!nrow(d)) return(make_empty_plot("HOR candidate scan", "No testable blocks", base_theme))
  d$is_candidate <- d$hor_call %in% opt$candidate_calls
  d$plot_call <- ifelse(d$is_candidate, d$hor_call, "NOT_DETECTED")
  if (opt$permutations > 0L && any(is.finite(d$block_fdr_q))) {
    floor_q <- 1 / ((opt$permutations + 1) * max(1L, nrow(d)))
    d$overview_y <- -log10(pmax(d$block_fdr_q, floor_q, na.rm = TRUE))
    y_label <- expression(-log[10]("BH-adjusted q"))
    subtitle <- "Candidates combine max-lag permutation significance, effect size, and multi-threshold stability"
    fdr_guides <- data.frame(
      y = -log10(c(opt$moderate_fdr, opt$high_fdr)),
      label = c("MODERATE FDR", "HIGH FDR")
    )
  } else {
    d$overview_y <- d$estimated_units
    y_label <- "Estimated HOR units in block"
    subtitle <- "Permutations disabled: calls reflect effect size and multi-threshold stability only"
    fdr_guides <- data.frame(y = numeric(), label = character())
  }
  d$label <- ifelse(
    d$is_candidate,
    paste0(d$species, " ", d$source_array_key, " p", d$candidate_period),
    ""
  )
  shapes <- c(HIGH = 21, MODERATE = 24, EFFECT_ONLY_HIGH = 21,
              EFFECT_ONLY_MODERATE = 24, NOT_DETECTED = 16)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = candidate_robust_excess, y = overview_y)) +
    ggplot2::geom_point(
      data = d[!d$is_candidate, , drop = FALSE],
      colour = "grey75", size = 1.8, alpha = 0.65
    ) +
    ggplot2::geom_point(
      data = d[d$is_candidate, , drop = FALSE],
      ggplot2::aes(fill = species, shape = plot_call),
      colour = "black", stroke = 0.35, size = 3.0
    ) +
    ggplot2::geom_text(
      data = d[d$is_candidate, , drop = FALSE],
      ggplot2::aes(label = label),
      hjust = -0.05, vjust = 0.15, size = 2.5, check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = species_colours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
    ggplot2::labs(
      title = paste0(satellite_family, " HOR-like periodicity scan"),
      subtitle = subtitle,
      x = "Robust periodicity excess (median across clustering thresholds)",
      y = y_label,
      fill = "Species", shape = "Call"
    ) +
    base_theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::coord_cartesian(clip = "off")
  if (nrow(fdr_guides)) {
    p <- p + ggplot2::geom_hline(
      data = fdr_guides, ggplot2::aes(yintercept = y),
      inherit.aes = FALSE, colour = "grey45", linetype = 2, linewidth = 0.35
    )
  }
  p
}

analyze_block <- function(d_all, selection_row, fasta_sequences, opt, block_index, base_theme, satellite_family) {
  d <- d_all[
    d_all$species == selection_row$species &
      d_all$chrom == selection_row$chrom &
      d_all$region_key == selection_row$region_key &
      d_all$source_array_key == selection_row$source_array_key,
    , drop = FALSE
  ]
  d <- d[order(d$start, d$end, d$seq_id), , drop = FALSE]
  d$full_block_copy_index <- seq_len(nrow(d))
  d <- d[
    d$full_block_copy_index >= selection_row$copy_start &
      d$full_block_copy_index <= selection_row$copy_end,
    , drop = FALSE
  ]
  d$copy_index <- seq_len(nrow(d))
  n <- nrow(d)
  if (n < 3L) stop("Fewer than three copies after filtering: ", selection_row$source_array_key)

  main_key <- threshold_key(opt$main_identity)
  main_column <- threshold_columns[[main_key]]
  main_threshold_label <- threshold_labels[[main_key]]
  main_labels <- as.character(d[[main_column]])
  main_labels[is.na(main_labels) | !nzchar(main_labels)] <- "Unassigned"
  first_seen <- unique(main_labels)
  d$cluster_rank <- match(main_labels, first_seen)
  occurrence <- table(main_labels)
  d$cluster_occurrences <- as.numeric(occurrence[main_labels])

  supplied_period <- selection_row$expected_period
  provisional_max_lag <- min(n - 1L, max(opt$max_lag, ifelse(is.finite(supplied_period), supplied_period, 1L)))
  minimum_supported_pairs <- max(
    as.integer(opt$min_lag_pairs),
    as.integer(ceiling(opt$min_pair_fraction * n))
  )
  profiles <- do.call(rbind, lapply(thresholds, function(threshold) {
    key <- threshold_key(threshold)
    lag_profile_one(
      d[[threshold_columns[[key]]]], provisional_max_lag, threshold,
      min_lag_pairs = opt$min_lag_pairs,
      min_pair_fraction = opt$min_pair_fraction
    )
  }))
  main_profile <- profiles[abs(profiles$threshold - opt$main_identity) < 1e-8, , drop = FALSE]
  supported <- main_profile[
    main_profile$meets_minimum_pairs & is.finite(main_profile$periodicity_excess),
    , drop = FALSE
  ]
  if (nrow(supported)) {
    supported_selection <- select_fundamental_lag(
      supported$lag,
      supported$periodicity_excess,
      absolute_tolerance = opt$fundamental_abs_tolerance,
      relative_tolerance = opt$fundamental_rel_tolerance
    )
    best_row <- supported[supported_selection$candidate_index, , drop = FALSE]
  } else {
    best_row <- data.frame(
      lag = NA_integer_, same_cluster_fraction = NA_real_,
      periodicity_excess = NA_real_, stringsAsFactors = FALSE
    )
  }
  if (is.finite(supplied_period) && supplied_period >= 1L) {
    expected_period <- as.integer(supplied_period)
    period_source <- if (identical(as.character(selection_row$selection_source), "automatic_scan")) {
      "automatic_global_scan"
    } else {
      "supplied"
    }
  } else {
    robust_scan <- scan_label_matrix(
      label_matrix_from_block(d),
      min_period = opt$min_period,
      max_period = min(opt$max_period, opt$max_lag),
      min_units = opt$min_units,
      min_lag_pairs = opt$min_lag_pairs,
      min_pair_fraction = opt$min_pair_fraction,
      support_excess = opt$support_excess,
      fundamental_abs_tolerance = opt$fundamental_abs_tolerance,
      fundamental_rel_tolerance = opt$fundamental_rel_tolerance
    )
    if (robust_scan$test_status != "TESTED") {
      stop("Cannot estimate an automatic period for ", selection_row$source_array_key,
           ": ", robust_scan$not_tested_reason)
    }
    expected_period <- robust_scan$candidate_period
    period_source <- "auto_multi_threshold_peak"
  }
  if (expected_period >= n) stop("Expected period must be smaller than selected copy count for ", selection_row$source_array_key)

  if (expected_period > provisional_max_lag) {
    block_max_lag <- min(n - 1L, expected_period)
    profiles <- do.call(rbind, lapply(thresholds, function(threshold) {
      key <- threshold_key(threshold)
      lag_profile_one(
        d[[threshold_columns[[key]]]], block_max_lag, threshold,
        min_lag_pairs = opt$min_lag_pairs,
        min_pair_fraction = opt$min_pair_fraction
      )
    }))
    main_profile <- profiles[abs(profiles$threshold - opt$main_identity) < 1e-8, , drop = FALSE]
  } else {
    block_max_lag <- provisional_max_lag
  }

  periodicity_rows <- lapply(seq_along(thresholds), function(k) {
    threshold <- thresholds[[k]]
    key <- threshold_key(threshold)
    set.seed(opt$seed + block_index * 100L + k)
    evidence_permutations <- if (opt$mode == "auto") 0L else opt$permutations
    permutation <- lag_permutation(
      d[[threshold_columns[[key]]]], expected_period, evidence_permutations,
      min_valid_pairs = minimum_supported_pairs
    )
    profile_row <- profiles[abs(profiles$threshold - threshold) < 1e-8 & profiles$lag == expected_period, , drop = FALSE]
    data.frame(
      block_label = selection_row$label,
      species = selection_row$species,
      chrom = selection_row$chrom,
      region_key = selection_row$region_key,
      source_array_key = selection_row$source_array_key,
      threshold = threshold,
      threshold_label = threshold_labels[[key]],
      expected_period = expected_period,
      n_lag_pairs = if (nrow(profile_row)) profile_row$n_pairs else 0L,
      minimum_lag_pairs = minimum_supported_pairs,
      lag_pair_requirement_met = if (nrow(profile_row)) profile_row$meets_minimum_pairs else FALSE,
      same_cluster_fraction = if (nrow(profile_row)) profile_row$same_cluster_fraction else NA_real_,
      random_order_expectation = if (nrow(profile_row)) profile_row$random_order_expectation else NA_real_,
      periodicity_excess = if (nrow(profile_row)) profile_row$periodicity_excess else NA_real_,
      permutation_null_mean = permutation[["null_mean"]],
      permutation_null_low = permutation[["null_low"]],
      permutation_null_high = permutation[["null_high"]],
      permutation_p_enrichment = permutation[["p_enrichment"]],
      stringsAsFactors = FALSE
    )
  })
  periodicity <- do.call(rbind, periodicity_rows)
  main_periodicity <- periodicity[abs(periodicity$threshold - opt$main_identity) < 1e-8, , drop = FALSE]

  matrix_indices <- seq_len(min(n, opt$max_matrix_copies))
  matrix_d <- d[matrix_indices, , drop = FALSE]
  matrix_result <- NULL
  matrix_mode <- "multi-threshold cluster agreement"
  if (!is.null(fasta_sequences)) {
    matrix_result <- exact_identity_matrix(as.character(matrix_d$haplotype_id), fasta_sequences)
    if (!is.null(matrix_result)) matrix_mode <- "fixed-reference pairwise identity"
  }
  if (is.null(matrix_result)) matrix_result <- cluster_agreement_matrix(matrix_d)
  similarity_matrix <- matrix_result$matrix
  matrix_n <- nrow(similarity_matrix)

  pair_index <- which(upper.tri(similarity_matrix), arr.ind = TRUE)
  pair_data <- data.frame(
    copy_i = pair_index[, 1L],
    copy_j = pair_index[, 2L],
    lag = pair_index[, 2L] - pair_index[, 1L],
    similarity = similarity_matrix[pair_index],
    stringsAsFactors = FALSE
  )
  pair_data <- pair_data[is.finite(pair_data$similarity) & pair_data$lag <= block_max_lag, , drop = FALSE]
  pair_data$phase_group <- ifelse(
    pair_data$lag %% expected_period == 0L,
    "Expected phase (lag = k × p)",
    "Other lags"
  )
  pair_data$phase_group <- factor(pair_data$phase_group, levels = c("Expected phase (lag = k × p)", "Other lags"))
  same_values <- pair_data$similarity[pair_data$phase_group == "Expected phase (lag = k × p)"]
  other_values <- pair_data$similarity[pair_data$phase_group == "Other lags"]
  same_median <- if (length(same_values)) median(same_values, na.rm = TRUE) else NA_real_
  other_median <- if (length(other_values)) median(other_values, na.rm = TRUE) else NA_real_
  similarity_delta <- same_median - other_median

  guide_multiples <- seq(expected_period, block_max_lag, by = expected_period)
  boundary_lines <- seq(expected_period + 0.5, n - 0.5, by = expected_period)
  species_colour <- species_colours[[selection_row$species]]
  if (is.null(species_colour) || is.na(species_colour)) species_colour <- "#555555"

  p_track <- ggplot2::ggplot(d, ggplot2::aes(x = copy_index, y = cluster_rank)) +
    ggplot2::geom_vline(xintercept = boundary_lines, colour = "grey75", linetype = 3, linewidth = 0.25) +
    ggplot2::geom_line(colour = "grey78", linewidth = 0.22) +
    ggplot2::geom_point(
      ggplot2::aes(size = log2(cluster_occurrences + 1)),
      shape = 21, fill = species_colour, colour = "white", stroke = 0.15, alpha = 0.9
    ) +
    ggplot2::scale_size_continuous(range = c(0.8, 3.0), name = "log2(cluster copies + 1)") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.005, 0.005))) +
    ggplot2::labs(
      title = paste0("A  Ordered ", main_threshold_label, " cluster recurrence"),
      subtitle = "Cluster ranks are assigned by first appearance; horizontal recurrences reveal repeated monomer positions",
      x = paste0(satellite_family, " copy index in genomic order"),
      y = "Variant-cluster rank"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "right", panel.grid.minor = ggplot2::element_blank())

  matrix_frame <- expand.grid(copy_x = seq_len(matrix_n), copy_y = seq_len(matrix_n))
  matrix_frame$value <- as.vector(similarity_matrix)
  heat_limits <- if (matrix_mode == "fixed-reference pairwise identity") c(opt$heatmap_min, 1) else c(min(thresholds) - 0.05, max(thresholds))
  heat_title <- if (matrix_mode == "fixed-reference pairwise identity") {
    "B  Pairwise monomer identity"
  } else {
    "B  Multi-threshold cluster agreement"
  }
  heat_subtitle <- if (matrix_mode == "fixed-reference pairwise identity") {
    paste0("Fixed-reference alignment; coverage ", format_percent(matrix_result$coverage), "; blue guides mark lag k × p")
  } else {
    "Highest shared 80/85/90/95% cluster level; this is not exact pairwise identity"
  }
  intercepts <- unique(c(guide_multiples, -guide_multiples))
  p_matrix <- ggplot2::ggplot(matrix_frame, ggplot2::aes(x = copy_x, y = copy_y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::geom_abline(intercept = intercepts, slope = 1, colour = "#2C7FB8", linetype = 2, linewidth = 0.18, alpha = 0.7) +
    ggplot2::scale_fill_gradientn(
      colours = heatmap_colours,
      limits = heat_limits,
      oob = scales::squish,
      na.value = "grey85",
      # Both exact identity and cluster-agreement levels are proportions.
      # Use an explicit formatter so this fallback branch is independent of
      # whether ggplot2 is attached to the search path.
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(
      title = heat_title,
      subtitle = heat_subtitle,
      x = if (matrix_n < n) paste0("Copy index (first ", matrix_n, ")") else "Copy index",
      y = if (matrix_n < n) paste0("Copy index (first ", matrix_n, ")") else "Copy index",
      fill = if (matrix_mode == "fixed-reference pairwise identity") "Identity" else "Cluster level"
    ) +
    base_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")

  profile_guides <- data.frame(lag = guide_multiples)
  expected_points <- profiles[profiles$lag %in% guide_multiples, , drop = FALSE]
  p_lag <- ggplot2::ggplot(
    profiles,
    ggplot2::aes(x = lag, y = periodicity_excess, colour = threshold_label, group = threshold_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linetype = 3, linewidth = 0.35) +
    ggplot2::geom_vline(data = profile_guides, ggplot2::aes(xintercept = lag), inherit.aes = FALSE, colour = "grey55", linetype = 2, linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.72, alpha = 0.9) +
    ggplot2::geom_point(data = expected_points, size = 1.7, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = threshold_colours, breaks = names(threshold_colours), drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = pretty(c(1, block_max_lag), n = 8)) +
    ggplot2::labs(
      title = "C  Lag periodicity across clustering thresholds",
      subtitle = paste0("Candidate p = ", expected_period, "; vertical guides show p, 2p, ..."),
      x = paste0("Lag between ", satellite_family, " copies"),
      y = "Same-cluster fraction − random-order expectation",
      colour = "Clustering identity"
    ) +
    base_theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")

  if (nrow(pair_data) > 0L && length(unique(pair_data$phase_group[!is.na(pair_data$phase_group)])) == 2L) {
    set.seed(opt$seed + block_index * 1000L)
    jitter_data <- do.call(rbind, lapply(split(pair_data, pair_data$phase_group), function(z) {
      z[sample(seq_len(nrow(z)), min(nrow(z), 1500L)), , drop = FALSE]
    }))
    metric_label <- if (matrix_mode == "fixed-reference pairwise identity") "Pairwise identity" else "Highest shared cluster level"
    permutation_text <- if (opt$mode == "auto") {
      "candidate significance was assessed by the global max-lag/BH scan"
    } else {
      paste0(main_threshold_label, " lag-p permutation P = ",
             format_p(main_periodicity$permutation_p_enrichment))
    }
    subtitle_text <- paste0(
      "Median same phase = ", if (matrix_mode == "fixed-reference pairwise identity") format_percent(same_median) else sprintf("%.2f", same_median),
      "; other = ", if (matrix_mode == "fixed-reference pairwise identity") format_percent(other_median) else sprintf("%.2f", other_median),
      "; Δ = ", if (matrix_mode == "fixed-reference pairwise identity") format_percent(similarity_delta) else sprintf("%.2f", similarity_delta),
      "; ", permutation_text
    )
    p_phase <- ggplot2::ggplot(pair_data, ggplot2::aes(x = phase_group, y = similarity, fill = phase_group)) +
      ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.4) +
      ggplot2::geom_jitter(data = jitter_data, width = 0.12, height = 0, alpha = 0.12, size = 0.45, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c("Expected phase (lag = k × p)" = species_colour, "Other lags" = "grey72"), guide = "none") +
      ggplot2::labs(
        title = "D  Monomers at the same expected HOR phase are more similar",
        subtitle = subtitle_text,
        x = NULL,
        y = metric_label
      ) +
      base_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 12, hjust = 1), panel.grid.minor = ggplot2::element_blank())
    if (matrix_mode == "fixed-reference pairwise identity") {
      p_phase <- p_phase + ggplot2::coord_cartesian(ylim = c(opt$heatmap_min, 1)) +
        ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"))
    } else {
      cluster_breaks <- c(min(thresholds) - 0.05, thresholds)
      p_phase <- p_phase + ggplot2::scale_y_continuous(
        breaks = cluster_breaks,
        labels = c("<80%", "80%", "85%", "90%", "95%"),
        limits = range(cluster_breaks)
      )
    }
  } else {
    p_phase <- make_empty_plot("D  Expected-phase similarity", "Insufficient finite pairwise values", base_theme)
  }

  period_word <- paste0(expected_period, "-mer candidate")
  page_title <- paste0(
    selection_row$species, " ", selection_row$chrom, " ", selection_row$region_key,
    " | block ", selection_row$source_array_key, " | ", period_word
  )
  page_subtitle <- paste0(
    "Selected copies ", selection_row$copy_start, "–", selection_row$copy_end,
    " (n = ", n, "; ", sprintf("%.2f", n / expected_period), " candidate units); period source: ", period_source
  )
  page_caption <- if (matrix_mode == "fixed-reference pairwise identity") {
    "A recurring ordered cluster pattern, off-diagonal identity bands at lag p, lag-p enrichment across thresholds, and greater same-phase identity jointly support an HOR. Read mapping across ≥2 HOR units remains the final assembly-level validation."
  } else {
    "Compact result input lacks the fixed-reference haplotype alignment, so panels B and D use cluster agreement rather than exact identity. Supply --alignment OUTDIR/prepared/exact_haplotypes.reference.aln.fa for the stronger sequence-identity version."
  }
  design <- "\nAA\nBC\nBD\n"
  page <- p_track + p_matrix + p_lag + p_phase +
    patchwork::plot_layout(design = design, widths = c(1.03, 1.25), heights = c(0.70, 1, 1)) +
    patchwork::plot_annotation(
      title = page_title,
      subtitle = page_subtitle,
      caption = page_caption,
      theme = ggplot2::theme(
        text = ggplot2::element_text(family = opt$font_family),
        plot.title = ggplot2::element_text(face = "bold", size = 15),
        plot.subtitle = ggplot2::element_text(size = 10),
        plot.caption = ggplot2::element_text(size = 8, colour = "grey30")
      )
    )

  lag_out <- profiles
  lag_out$satellite <- satellite_family
  lag_out$block_label <- selection_row$label
  lag_out$species <- selection_row$species
  lag_out$chrom <- selection_row$chrom
  lag_out$region_key <- selection_row$region_key
  lag_out$source_array_key <- selection_row$source_array_key
  lag_out$expected_period <- expected_period

  selected_copy_out <- d[, unique(c(
    "seq_id", "species", "chrom", "region_key", "source_array_key", "array_key",
    "start", "end", "haplotype_id", "copy_index", "full_block_copy_index",
    unname(threshold_columns), "cluster_rank", "cluster_occurrences"
  )), drop = FALSE]
  selected_copy_out$block_label <- selection_row$label
  selected_copy_out$expected_period <- expected_period
  selected_copy_out$satellite <- satellite_family

  summary <- data.frame(
    satellite = satellite_family,
    block_label = selection_row$label,
    species = selection_row$species,
    chrom = selection_row$chrom,
    region_key = selection_row$region_key,
    source_array_key = selection_row$source_array_key,
    selection_source = selection_row$selection_source,
    selected_copy_start = selection_row$copy_start,
    selected_copy_end = selection_row$copy_end,
    n_copies = n,
    expected_period = expected_period,
    period_source = period_source,
    candidate_unit_count = n / expected_period,
    main_identity = opt$main_identity,
    best_supported_lag = best_row$lag,
    best_supported_same_cluster_fraction = best_row$same_cluster_fraction,
    best_supported_periodicity_excess = best_row$periodicity_excess,
    expected_lag_same_cluster_fraction = main_periodicity$same_cluster_fraction,
    expected_lag_random_order_expectation = main_periodicity$random_order_expectation,
    expected_lag_periodicity_excess = main_periodicity$periodicity_excess,
    expected_lag_permutation_p = main_periodicity$permutation_p_enrichment,
    minimum_pairs_for_best_lag = minimum_supported_pairs,
    matrix_mode = matrix_mode,
    matrix_n_copies = matrix_n,
    alignment_haplotype_coverage = matrix_result$coverage,
    aligned_reference_columns = matrix_result$aligned_columns,
    same_phase_pair_count = length(same_values),
    other_lag_pair_count = length(other_values),
    same_phase_median_similarity = same_median,
    other_lag_median_similarity = other_median,
    median_similarity_delta = similarity_delta,
    stringsAsFactors = FALSE
  )

  list(
    page = page,
    summary = summary,
    periodicity = periodicity,
    lag_profile = lag_out,
    selected_copies = selected_copy_out,
    file_label = safe_filename(selection_row$label)
  )
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  resolved <- resolve_input(opt$input)
  if (!is.null(resolved$temporary)) {
    on.exit(unlink(resolved$temporary, recursive = TRUE, force = TRUE), add = TRUE)
  }
  joined <- read_and_join_assignments(resolved$root)
  copies <- joined$data
  length_values <- if ("local_consensus_length" %in% names(copies)) {
    copies$local_consensus_length
  } else if ("phased_length" %in% names(copies)) {
    copies$phased_length
  } else {
    copies$end - copies$start + 1
  }
  satellite_family <- resolve_satellite_family(opt$satellite, length_values, opt$input)
  if (is.null(opt$prefix)) opt$prefix <- satellite_family
  if (is.null(opt$outdir)) opt$outdir <- paste0(satellite_family, "_HOR_analysis")
  message("Satellite family: ", satellite_family)
  message("Analysis mode:    ", opt$mode)
  available <- summarize_available_blocks(copies)
  available$satellite <- satellite_family
  available <- available[, c("satellite", setdiff(names(available), "satellite")), drop = FALSE]

  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  available_file <- file.path(opt$outdir, paste0(opt$prefix, ".available_blocks.tsv"))
  write_tsv(available, available_file)
  if (opt$list_blocks) {
    message("Available blocks written: ", available_file)
    return(invisible(NULL))
  }

  scan_result <- NULL
  if (opt$mode == "auto") {
    scan_result <- scan_all_blocks(copies, available, opt, satellite_family)
    scan_all_file <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.all_blocks.tsv"))
    scan_candidate_file <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.candidates.tsv"))
    scan_region_file <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.regions.tsv"))
    scan_lag_file <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.all_lags.tsv"))
    scan_config_file <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.evidence_config.tsv"))
    scan_metadata_file <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.metadata.tsv"))
    write_tsv(scan_result$summary, scan_all_file)
    write_tsv(scan_result$candidates, scan_candidate_file)
    write_tsv(make_region_summary(scan_result$summary), scan_region_file)
    write_tsv(scan_result$all_lags, scan_lag_file)
    evidence_config <- make_evidence_config(scan_result$candidates)
    write_tsv(evidence_config, scan_config_file)
    write_tsv(make_scan_metadata(opt, satellite_family, scan_result), scan_metadata_file)
    message("Candidate scan complete: ", nrow(scan_result$candidates), " candidate blocks among ",
            sum(scan_result$summary$test_status == "TESTED"), " tested blocks")
    if (!opt$write_pdf && !opt$write_png) {
      message("Plots disabled; scan tables written to: ", opt$outdir)
      return(invisible(scan_result))
    }
    selection <- normalize_selection(evidence_config, "automatic_scan")
  } else {
    custom_present <- length(opt$blocks) > 0L || !is.null(opt$config)
    if (!is.null(opt$preset) && satellite_family != "SAT113" && opt$preset != "none") {
      stop("The strong/extended legacy presets apply only to SAT113. Use --block or --config for SAT483 manual mode.")
    }
    preset_name <- if (!is.null(opt$preset)) {
      opt$preset
    } else if (custom_present) {
      "none"
    } else if (satellite_family == "SAT113") {
      "strong"
    } else {
      stop("SAT483 manual mode requires --block or --config; use --mode auto for discovery.")
    }
    preset <- preset_selection(preset_name)
    inline <- parse_inline_blocks(opt$blocks)
    config <- read_config_selection(opt$config)
    selection <- rbind(preset, inline, config)
  }
  if (!nrow(selection)) {
    message("No candidate blocks met the requested criteria; scan tables were still written.")
    return(invisible(scan_result))
  }
  resolved_selection <- expand_selection(selection, available, opt$min_block_copies)
  selection_file <- file.path(opt$outdir, paste0(opt$prefix, ".resolved_blocks.tsv"))
  write_tsv(resolved_selection, selection_file)

  missing_packages <- c("ggplot2", "patchwork", "scales")[
    !vapply(c("ggplot2", "patchwork", "scales"), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing R packages: ", paste(missing_packages, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
    )
  }

  alignment_file <- NULL
  if (!is.null(opt$alignment) && tolower(opt$alignment) != "none") {
    if (!file.exists(opt$alignment)) stop("Alignment FASTA does not exist: ", opt$alignment)
    alignment_file <- normalizePath(opt$alignment, mustWork = TRUE)
  } else if (is.null(opt$alignment)) {
    alignment_file <- find_exact_file(
      resolved$root,
      c("exact_haplotypes.reference.aln.fa", "exact_haplotypes.reference_phase.aln.fa"),
      required = FALSE
    )
  }
  fasta_sequences <- NULL
  if (!is.null(alignment_file)) {
    message("  Haplotype alignment: ", alignment_file)
    fasta_sequences <- read_fasta(alignment_file)
  } else {
    message("  Haplotype alignment: not found; using multi-threshold cluster agreement")
  }

  base_theme <- ggplot2::theme_bw(base_size = 9.5, base_family = opt$font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "grey25"),
      axis.title = ggplot2::element_text(size = 8.8),
      axis.text = ggplot2::element_text(size = 7.7),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7.4)
    )

  if (opt$mode == "auto") {
    overview <- make_scan_overview_plot(scan_result$summary, satellite_family, opt, base_theme)
    if (opt$write_pdf) {
      ggplot2::ggsave(file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.overview.pdf")),
                      overview, width = 11.5, height = 7.5, units = "in", useDingbats = FALSE)
    }
    if (opt$write_png) {
      ggplot2::ggsave(file.path(opt$outdir, paste0(opt$prefix, ".HOR_scan.overview.png")),
                      overview, width = 11.5, height = 7.5, units = "in", dpi = opt$dpi, limitsize = FALSE)
    }
  }

  pdf_path <- file.path(opt$outdir, paste0(opt$prefix, ".HOR_evidence.pdf"))
  pdf_open <- FALSE
  if (opt$write_pdf) {
    grDevices::pdf(pdf_path, width = 14, height = 11, onefile = TRUE, useDingbats = FALSE)
    pdf_open <- TRUE
    on.exit({
      if (isTRUE(pdf_open)) try(grDevices::dev.off(), silent = TRUE)
    }, add = TRUE)
  }
  summaries <- list()
  periodicities <- list()
  lag_profiles <- list()
  selected_copies <- list()

  for (i in seq_len(nrow(resolved_selection))) {
    row <- resolved_selection[i, , drop = FALSE]
    message(
      "[", i, "/", nrow(resolved_selection), "] ", row$species, " ",
      row$source_array_key, " copies ", row$copy_start, "-", row$copy_end,
      if (is.finite(row$expected_period)) paste0(" p=", row$expected_period) else " p=auto"
    )
    result <- analyze_block(copies, row, fasta_sequences, opt, i, base_theme, satellite_family)
    if (pdf_open) print(result$page)
    if (opt$write_png) {
      png_path <- file.path(opt$outdir, paste0(opt$prefix, ".", result$file_label, ".HOR_evidence.png"))
      ggplot2::ggsave(
        png_path,
        result$page,
        width = 14,
        height = 11,
        units = "in",
        dpi = opt$dpi,
        limitsize = FALSE
      )
    }
    summaries[[i]] <- result$summary
    periodicities[[i]] <- result$periodicity
    lag_profiles[[i]] <- result$lag_profile
    selected_copies[[i]] <- result$selected_copies
  }
  if (pdf_open) {
    grDevices::dev.off()
    pdf_open <- FALSE
  }

  summary_file <- file.path(opt$outdir, paste0(opt$prefix, ".block_summary.tsv"))
  periodicity_file <- file.path(opt$outdir, paste0(opt$prefix, ".periodicity_at_expected_lag.tsv"))
  lag_file <- file.path(opt$outdir, paste0(opt$prefix, ".lag_profiles.tsv"))
  copies_file <- file.path(opt$outdir, paste0(opt$prefix, ".selected_copy_assignments.tsv"))
  write_tsv(do.call(rbind, summaries), summary_file)
  write_tsv(do.call(rbind, periodicities), periodicity_file)
  write_tsv(do.call(rbind, lag_profiles), lag_file)
  write_tsv(do.call(rbind, selected_copies), copies_file)

  message("Done.")
  if (opt$write_pdf) message("  Multipage evidence PDF: ", pdf_path)
  if (opt$write_png) message("  Per-block PNG files:     ", opt$outdir)
  message("  Block summary:           ", summary_file)
  message("  Threshold sensitivity:   ", periodicity_file)
  message("  Lag profiles:            ", lag_file)
  message("  Resolved block config:   ", selection_file)
}

if (!identical(Sys.getenv("SATELLITE_HOR_SOURCE_ONLY"), "1")) main()
