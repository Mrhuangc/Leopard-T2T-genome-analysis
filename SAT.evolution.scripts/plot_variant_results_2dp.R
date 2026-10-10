#!/usr/bin/env Rscript

# Publication-oriented plotting of cross-species satellite-variant analyses.
#
# This script replaces the Python/Matplotlib plotting step. It reads the TSV
# files produced by run_sat113_variant_origin.sh / variant_origin v2.3 and writes:
#   1. R equivalents of the original Python figures; and
#   2. the reorganized main figure and Supplementary Figures S1-S6.
#
# Required R packages:
#   ggplot2, patchwork, scales
#
# Simplest examples:
#   Rscript plot_SAT113_variant_results.R SAT483.variant_origin.zip
#   Rscript plot_SAT113_variant_results.R SAT113.variant_origin/results
#
# The satellite name is inferred from the input path (for example, SAT113 or
# SAT483). Use --satellite to override the inferred name when the input path is
# generic.

suppressPackageStartupMessages({
  required_packages <- c("ggplot2", "patchwork", "scales")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing R package(s): ", paste(missing_packages, collapse = ", "),
      "\nInstall them with:\n",
      "  install.packages(c(",
      paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
      "))",
      call. = FALSE
    )
  }
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

# Shared palette for every quantitative heatmap. Values increase from white
# through yellow to red, so empty/low-abundance cells remain visually quiet and
# high-abundance/high-frequency cells are emphasized consistently.
HEATMAP_COLORS <- c("#FFFFFF", "#FFFFB2", "#FED976", "#FD8D3C", "#BD0026")
HEATMAP_NA_COLOR <- "#F2F2F2"

usage <- function(status = 0L) {
  txt <- paste0(
    "Usage:\n",
    "  Rscript plot_SAT113_variant_results.R RUN_OUTPUT_DIR_OR_ZIP [options]\n",
    "  Rscript plot_SAT113_variant_results.R --input RUN_OUTPUT_DIR_OR_ZIP [options]\n\n",
    "Input (choose one form):\n",
    "  PATH                        Positional run output directory/root or zip file.\n",
    "  --input PATH                Same as positional PATH.\n",
    "  --primary-dir PATH          Backward-compatible alias for --input.\n\n",
    "Optional:\n",
    "  --sensitivity-dir PATH      Combined sensitivity result directory/root or zip.\n",
    "  --satellite NAME            Satellite label; inferred from input path when omitted.\n",
    "  --outdir DIR                Output directory [<SATELLITE>_R_figures].\n",
    "  --main-threshold FLOAT      Variant-cluster identity [0.85].\n",
    "  --dpi INT                   PNG resolution [400].\n",
    "  --formats LIST              Comma-separated pdf,png,svg [pdf,png].\n",
    "  --top-variants INT          Variant clusters in heatmap [40].\n",
    "  --top-regions-per-species INT  Regions per species in heatmap [25].\n",
    "  --layout MODE               all, manuscript, or legacy/original [all].\n",
    "  --primary-label TEXT        Short legend label [Primary].\n",
    "  --sensitivity-label TEXT    Short legend label [Sensitivity].\n",
    "  -h, --help                  Show this help.\n\n",
    "The result path may point to the run root, its results/ directory, or a zip.\n",
    "With one input, all figures supported by that run are redrawn in R.\n",
    "Use --layout legacy (or original) for only the R replacements of figures\n",
    "created by run_sat113_variant_origin.sh. Supp. Figures S1, S5, and the\n",
    "two-run S6 require --sensitivity-dir.\n"
  )
  cat(txt)
  quit(save = "no", status = status)
}

parse_cli <- function(args) {
  cfg <- list(
    primary_dir = NULL,
    sensitivity_dir = NULL,
    satellite = NULL,
    outdir = NULL,
    main_threshold = 0.85,
    dpi = 400L,
    formats = c("pdf", "png"),
    top_variants = 40L,
    top_regions_per_species = 25L,
    layout = "all",
    primary_label = "Primary",
    sensitivity_label = "Sensitivity"
  )

  if (length(args) == 0L) usage(1L)
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (token %in% c("-h", "--help")) usage(0L)
    if (!startsWith(token, "--")) {
      positional <- c(positional, token)
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop("Missing value after ", token, call. = FALSE)
    }
    key <- gsub("-", "_", substring(token, 3L), fixed = TRUE)
    value <- args[[i + 1L]]
    if (key %in% c("input", "primary_dir", "result_dir", "run_output")) {
      cfg$primary_dir <- value
    } else if (key == "sensitivity_dir") {
      cfg$sensitivity_dir <- value
    } else if (key %in% c("satellite", "satellite_name", "repeat_name")) {
      cfg$satellite <- value
    } else if (key == "outdir") {
      cfg$outdir <- value
    } else if (key == "main_threshold") {
      cfg$main_threshold <- as.numeric(value)
    } else if (key == "dpi") {
      cfg$dpi <- as.integer(value)
    } else if (key == "formats") {
      cfg$formats <- unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
    } else if (key == "top_variants") {
      cfg$top_variants <- as.integer(value)
    } else if (key == "top_regions_per_species") {
      cfg$top_regions_per_species <- as.integer(value)
    } else if (key == "layout") {
      cfg$layout <- tolower(value)
    } else if (key == "primary_label") {
      cfg$primary_label <- value
    } else if (key == "sensitivity_label") {
      cfg$sensitivity_label <- value
    } else {
      stop("Unknown option: ", token, call. = FALSE)
    }
    i <- i + 2L
  }

  if (length(positional) > 1L) {
    stop("Only one positional run-output path is allowed.", call. = FALSE)
  }
  if (length(positional) == 1L) {
    if (!is.null(cfg$primary_dir) && !identical(cfg$primary_dir, positional[[1L]])) {
      stop("Provide the run output either positionally or with --input/--primary-dir, not both.", call. = FALSE)
    }
    cfg$primary_dir <- positional[[1L]]
  }
  if (is.null(cfg$primary_dir)) {
    stop("A run-output path is required (positional PATH, --input, or --primary-dir).", call. = FALSE)
  }
  if (!is.null(cfg$satellite) && !nzchar(trimws(cfg$satellite))) {
    stop("--satellite must not be empty.", call. = FALSE)
  }
  if (!is.finite(cfg$main_threshold) || cfg$main_threshold <= 0 || cfg$main_threshold > 1) {
    stop("--main-threshold must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(cfg$dpi) || cfg$dpi < 72L) stop("--dpi must be >= 72.", call. = FALSE)
  if (cfg$layout %in% c("original", "pipeline")) cfg$layout <- "legacy"
  if (!cfg$layout %in% c("all", "manuscript", "legacy")) {
    stop("--layout must be all, manuscript, or legacy/original.", call. = FALSE)
  }
  allowed_formats <- c("pdf", "png", "svg")
  if (length(cfg$formats) == 0L || any(!cfg$formats %in% allowed_formats)) {
    stop("--formats may contain only pdf,png,svg.", call. = FALSE)
  }
  cfg
}

cfg <- parse_cli(commandArgs(trailingOnly = TRUE))

extract_satellite_tokens <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(character())
  hits <- regmatches(toupper(x), gregexpr("SAT[0-9]+", toupper(x), perl = TRUE))[[1L]]
  unique(hits[hits != ""])
}

infer_satellite_reference_length <- function(path) {
  conclusion <- NULL
  if (file.exists(path) && grepl("\\.zip$", path, ignore.case = TRUE)) {
    listing <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) NULL)
    if (!is.null(listing) && "Name" %in% names(listing)) {
      entries <- listing$Name[
        grepl("(^|/)subtype_validation_conclusion\\.tsv$", listing$Name)
      ]
      entries <- entries[!grepl("__MACOSX", entries, fixed = TRUE)]
      if (length(entries) > 0L) {
        conclusion <- tryCatch(
          read.delim(
            unz(path, entries[[1L]]),
            header = TRUE,
            sep = "\t",
            quote = "",
            comment.char = "",
            check.names = FALSE,
            stringsAsFactors = FALSE
          ),
          error = function(e) NULL
        )
      }
    }
  } else if (dir.exists(path)) {
    hits <- list.files(
      path,
      pattern = "^subtype_validation_conclusion\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
    hits <- hits[!grepl("__MACOSX", hits, fixed = TRUE)]
    if (length(hits) > 0L) {
      conclusion <- tryCatch(
        read.delim(
          hits[[1L]],
          header = TRUE,
          sep = "\t",
          quote = "",
          comment.char = "",
          check.names = FALSE,
          stringsAsFactors = FALSE
        ),
        error = function(e) NULL
      )
    }
  }
  if (
    is.null(conclusion) ||
      !all(c("metric", "value") %in% names(conclusion)) ||
      !any(conclusion$metric == "fixed_reference_length")
  ) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(
    conclusion$value[conclusion$metric == "fixed_reference_length"][[1L]]
  ))
  if (!is.finite(value) || value < 10 || value > 100000 || abs(value - round(value)) > 1e-8) {
    return(NA_real_)
  }
  value
}

detect_satellite_name <- function(path) {
  # Prefer the user-visible input path. This correctly handles legacy SAT483
  # archives whose internal Python figure names were still hard-coded as SAT113.
  path_candidates <- c(path, normalizePath(path, mustWork = FALSE))
  from_path <- unique(unlist(lapply(path_candidates, extract_satellite_tokens), use.names = FALSE))
  reference_length <- infer_satellite_reference_length(path)
  if (length(from_path) > 0L) {
    path_length <- suppressWarnings(as.numeric(sub("^SAT", "", from_path[[1L]])))
    tolerance <- if (is.finite(reference_length)) max(2, 0.02 * reference_length) else NA_real_
    if (is.finite(reference_length) && is.finite(path_length) && abs(path_length - reference_length) > tolerance) {
      warning(
        "Satellite label inferred from path (", from_path[[1L]],
        ") differs from subtype fixed_reference_length (", reference_length,
        "); using the path label. Use --satellite to override."
      )
    }
    return(from_path[[1L]])
  }

  internal_names <- character()
  if (file.exists(path) && grepl("\\.zip$", path, ignore.case = TRUE)) {
    listing <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) NULL)
    if (!is.null(listing) && "Name" %in% names(listing)) internal_names <- listing$Name
  } else if (dir.exists(path)) {
    internal_names <- list.files(path, recursive = TRUE, full.names = FALSE)
  }
  internal_tokens <- unique(unlist(lapply(internal_names, extract_satellite_tokens), use.names = FALSE))
  if (length(internal_tokens) == 1L && is.finite(reference_length)) {
    internal_length <- suppressWarnings(as.numeric(sub("^SAT", "", internal_tokens[[1L]])))
    tolerance <- max(2, 0.02 * reference_length)
    if (is.finite(internal_length) && abs(internal_length - reference_length) <= tolerance) {
      return(internal_tokens[[1L]])
    }
    inferred <- paste0("SAT", as.integer(round(reference_length)))
    warning(
      "Internal figure label ", internal_tokens[[1L]],
      " is inconsistent with subtype fixed_reference_length ", reference_length,
      "; using ", inferred, ". Use --satellite to override."
    )
    return(inferred)
  }
  if (is.finite(reference_length)) {
    return(paste0("SAT", as.integer(round(reference_length))))
  }
  if (length(internal_tokens) == 1L) return(internal_tokens[[1L]])
  if (length(internal_tokens) > 1L) {
    warning(
      "Multiple satellite labels were found under the input (",
      paste(internal_tokens, collapse = ", "),
      "); using ", internal_tokens[[1L]], ". Use --satellite to override."
    )
    return(internal_tokens[[1L]])
  }
  NA_character_
}

normalize_satellite_name <- function(x) {
  value <- toupper(trimws(as.character(x)))
  value <- gsub("[[:space:]]+", "", value)
  if (!nzchar(value)) stop("Satellite label is empty.", call. = FALSE)
  value
}

if (is.null(cfg$satellite)) {
  cfg$satellite <- detect_satellite_name(cfg$primary_dir)
}
if (is.na(cfg$satellite) || !nzchar(cfg$satellite)) {
  cfg$satellite <- "SAT113"
  warning(
    "Could not infer a satellite label from the input path; using SAT113. ",
    "Use --satellite NAME to set it explicitly."
  )
}
cfg$satellite <- normalize_satellite_name(cfg$satellite)
satellite_label <- cfg$satellite
satellite_stem <- gsub("[^A-Z0-9._-]+", "_", satellite_label)
satellite_stem <- gsub("^_+|_+$", "", satellite_stem)
if (!nzchar(satellite_stem)) stop("Satellite label cannot be converted to a safe file prefix.", call. = FALSE)
if (is.null(cfg$outdir)) cfg$outdir <- paste0(satellite_stem, "_R_figures")

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required input file not found: ", path, call. = FALSE)
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "nan")
  )
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

threshold_tag <- function(x) {
  out <- format(x, scientific = FALSE, trim = TRUE, digits = 12)
  out <- sub("0+$", "", out)
  sub("\\.$", "", out)
}

resolve_result_dir <- function(path, label) {
  path <- path.expand(trimws(path))
  if (!file.exists(path) && !dir.exists(path)) {
    stop(
      label, " path does not exist: ", path,
      "\nCurrent working directory: ", getwd(),
      call. = FALSE
    )
  }

  search_root <- path
  if (dir.exists(path)) {
    search_root <- path
  } else if (file.exists(path) && grepl("\\.zip$", path, ignore.case = TRUE)) {
    unzip_root <- file.path(
      tempdir(),
      paste0(satellite_stem, "_R_", label, "_", as.integer(Sys.time()), "_", sample.int(1e6, 1L))
    )
    dir.create(unzip_root, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(path, exdir = unzip_root)
    search_root <- unzip_root
  } else {
    stop(label, " must be a directory or .zip file: ", path, call. = FALSE)
  }

  direct_candidates <- c(search_root, file.path(search_root, "results"))
  for (candidate in direct_candidates) {
    if (file.exists(file.path(candidate, "threshold_sensitivity.tsv"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  hits <- list.files(
    search_root,
    pattern = "^threshold_sensitivity\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  hits <- hits[!grepl("__MACOSX", hits, fixed = TRUE)]
  if (length(hits) == 0L) {
    stop("Cannot locate results/threshold_sensitivity.tsv under ", path, call. = FALSE)
  }
  if (length(hits) > 1L) {
    depth <- lengths(strsplit(normalizePath(hits), .Platform$file.sep, fixed = TRUE))
    hits <- hits[order(depth, nchar(hits))]
    warning("Multiple result tables found under ", path, "; using ", dirname(hits[[1L]]))
  }
  normalizePath(dirname(hits[[1L]]), mustWork = TRUE)
}

find_cluster_file <- function(result_dir, threshold) {
  direct <- file.path(
    result_dir,
    paste0("variant_cluster_summary.identity_", threshold_tag(threshold), ".tsv")
  )
  if (file.exists(direct)) return(direct)
  candidates <- list.files(
    result_dir,
    pattern = "^variant_cluster_summary\\.identity_.*\\.tsv$",
    full.names = TRUE
  )
  values <- suppressWarnings(as.numeric(sub(
    "^variant_cluster_summary\\.identity_(.*)\\.tsv$", "\\1", basename(candidates)
  )))
  keep <- which(is.finite(values) & abs(values - threshold) < 1e-9)
  if (length(keep) == 0L) {
    stop("No variant-cluster summary found for identity ", threshold, " in ", result_dir, call. = FALSE)
  }
  candidates[[keep[[1L]]]]
}

metric_value <- function(metric_table, metric, default = NA_real_) {
  hit <- metric_table$value[metric_table$metric == metric]
  if (length(hit) == 0L) return(default)
  suppressWarnings(as.numeric(hit[[1L]]))
}

conclusion_value <- function(conclusion, metric, default = NA_character_) {
  if (is.null(conclusion) || nrow(conclusion) == 0L) return(default)
  hit <- conclusion$value[conclusion$metric == metric]
  if (length(hit) == 0L) default else as.character(hit[[1L]])
}

load_analysis <- function(path, label, threshold) {
  result_dir <- resolve_result_dir(path, label)
  subtype_dir <- file.path(result_dir, "subtypes")

  dat <- list(
    label = label,
    result_dir = result_dir,
    threshold = read_tsv(file.path(result_dir, "threshold_sensitivity.tsv")),
    arrays = read_tsv(file.path(result_dir, "array_variant_summary.tsv")),
    pairs = read_tsv(file.path(result_dir, "array_pair_similarity.tsv")),
    linkage = read_tsv(file.path(result_dir, "site_pair_linkage.tsv")),
    linkage_summary = read_tsv(file.path(result_dir, "site_pair_linkage.summary.tsv")),
    metrics = read_tsv(file.path(result_dir, "hypothesis_diagnostic_metrics.tsv")),
    pca = read_tsv(file.path(result_dir, "haplotype_pca.tsv")),
    species_pairs = read_tsv(file.path(result_dir, "species_pair_origin_metrics.tsv")),
    rarefaction = read_tsv(file.path(result_dir, "rarefaction_species_sharing.tsv")),
    clusters = read_tsv(find_cluster_file(result_dir, threshold)),
    variant_composition = read_tsv(file.path(result_dir, "array_variant_composition.long.tsv")),
    core_columns = read_tsv(file.path(result_dir, "alignment_core_columns.tsv"))
  )

  subtype_files <- c(
    model = "subtype_model_selection.tsv",
    cv = "subtype_cross_validation.tsv",
    summary = "subtype_summary.tsv",
    diagnostics = "subtype_diagnostic_sites.tsv",
    composition = "subtype_array_composition.tsv",
    conclusion = "subtype_validation_conclusion.tsv",
    copy_assignments = "subtype_copy_assignments.tsv"
  )
  dat$subtypes <- list()
  for (name in names(subtype_files)) {
    f <- file.path(subtype_dir, subtype_files[[name]])
    dat$subtypes[[name]] <- if (file.exists(f)) read_tsv(f) else NULL
  }
  dat$recommended_k <- suppressWarnings(as.integer(conclusion_value(
    dat$subtypes$conclusion, "recommended_k", NA_character_
  )))
  if (!is.finite(dat$recommended_k) && !is.null(dat$subtypes$model)) {
    recommended <- as_flag(dat$subtypes$model$recommended)
    if (any(recommended)) dat$recommended_k <- dat$subtypes$model$k[which(recommended)[[1L]]]
  }
  dat
}

primary <- load_analysis(cfg$primary_dir, cfg$primary_label, cfg$main_threshold)
sensitivity <- if (!is.null(cfg$sensitivity_dir)) {
  load_analysis(cfg$sensitivity_dir, cfg$sensitivity_label, cfg$main_threshold)
} else {
  NULL
}

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
cfg$outdir <- normalizePath(cfg$outdir, mustWork = TRUE)

# One authoritative palette for every geometry that encodes species identity.
# Keep species geoms opaque: using different alpha values blends the same hex
# color with white and makes panels appear to use different saturation levels.
species_colors <- c(FCA = "#D55E00", PPA = "#0072B2", PTI = "#009E73")
species_geom_alpha <- 1
all_group_color <- "#7A7A7A"
analysis_colors <- setNames(c("#3B6FB6", "#C65D35"), c(cfg$primary_label, cfg$sensitivity_label))

# Publication number formatting:
# continuous quantities are shown with exactly two decimal places;
# counts remain integers.
label_continuous_2dp <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 2))
}

label_percent_2dp <- function(x) {
  ifelse(is.na(x), "NA", paste0(formatC(100 * x, format = "f", digits = 2), "%"))
}

label_percent_value_2dp <- function(x) {
  ifelse(is.na(x), "NA", paste0(formatC(x, format = "f", digits = 2), "%"))
}

label_count_integer <- function(x) {
  ifelse(is.na(x), "NA", scales::comma(round(x), accuracy = 1))
}
sharing_colors <- c(
  "Multi-region cluster" = "#0072B2",
  "Cross-species cluster" = "#D55E00",
  "Species-specific cluster" = "#009E73"
)
metric_colors <- c(
  "Consensus divergence" = "#7A7A7A",
  "Within-species pi" = "#56B4E9",
  "Between-species Dxy" = "#D55E00",
  "Net divergence Da" = "#009E73"
)

theme_sat113 <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.30, color = "black"),
      panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.28),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.4, color = "#444444"),
      plot.caption = element_text(size = base_size - 1.2, color = "#555555", hjust = 0),
      plot.tag = element_text(face = "bold", size = base_size + 3),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      strip.text = element_text(face = "bold", color = "black")
    )
}

empty_plot <- function(title, message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, color = "#777777", size = 3.2) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0))
}

saved_manifest <- data.frame(
  file = character(),
  description = character(),
  stringsAsFactors = FALSE
)

save_figure <- function(plot, stem, width, height, description) {
  for (fmt in cfg$formats) {
    path <- file.path(cfg$outdir, paste0(stem, ".", fmt))
    if (fmt == "pdf") {
      grDevices::pdf(path, width = width, height = height, useDingbats = FALSE, family = "sans")
      print(plot)
      grDevices::dev.off()
    } else if (fmt == "png") {
      png_type <- if (capabilities("cairo")) "cairo" else getOption("bitmapType", "Xlib")
      grDevices::png(
        path,
        width = width,
        height = height,
        units = "in",
        res = cfg$dpi,
        type = png_type,
        bg = "white"
      )
      print(plot)
      grDevices::dev.off()
    } else if (fmt == "svg") {
      grDevices::svg(path, width = width, height = height, family = "sans", bg = "white")
      print(plot)
      grDevices::dev.off()
    }
    saved_manifest <<- rbind(
      saved_manifest,
      data.frame(file = normalizePath(path), description = description, stringsAsFactors = FALSE)
    )
  }
}

species_scale_values <- function(x) {
  observed <- unique(as.character(x))
  missing <- setdiff(observed, names(species_colors))
  if (length(missing) == 0L) return(species_colors)
  extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
  c(species_colors, setNames(extra, missing))
}

pretty_class <- function(x) {
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("^shared ", "Shared: ", x)
  x <- gsub(" specific$", "-specific", x)
  x <- gsub("^rare unclassified$", "Rare/unclassified", x)
  x
}

plot_threshold_sensitivity <- function(dat, main_threshold) {
  d <- dat$threshold
  long <- rbind(
    data.frame(
      identity = 100 * as.numeric(d$identity_threshold),
      series = "Multi-region cluster",
      represented = 100 * as.numeric(d$fraction_copies_cross_array_cluster)
    ),
    data.frame(
      identity = 100 * as.numeric(d$identity_threshold),
      series = "Cross-species cluster",
      represented = 100 * as.numeric(d$fraction_copies_cross_species_cluster)
    ),
    data.frame(
      identity = 100 * as.numeric(d$identity_threshold),
      series = "Species-specific cluster",
      represented = 100 * as.numeric(d$fraction_copies_species_specific_cluster)
    )
  )
  long$series <- factor(long$series, levels = names(sharing_colors))
  ggplot(long, aes(identity, represented, color = series, group = series)) +
    geom_vline(xintercept = 100 * main_threshold, linetype = 2, linewidth = 0.45, color = "#555555") +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    scale_color_manual(values = sharing_colors, drop = FALSE) +
    scale_x_continuous(
      breaks = sort(unique(long$identity)),
      labels = label_continuous_2dp
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = label_continuous_2dp,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Variant-clustering identity (%)",
      y = "Copies represented (%)",
      title = "Threshold sensitivity"
    ) +
    theme_sat113() +
    theme(legend.position = c(0.03, 0.05), legend.justification = c(0, 0))
}

plot_exact_vs_cluster <- function(dat, main_threshold) {
  m <- dat$metrics
  cluster_label <- sprintf("Variant cluster (%d%%)", round(100 * main_threshold))
  d <- data.frame(
    scope = factor(
      rep(c("Across regions", "Across species"), each = 2L),
      levels = c("Across regions", "Across species")
    ),
    definition = factor(
      rep(c("Exact haplotype", cluster_label), 2L),
      levels = c("Exact haplotype", cluster_label)
    ),
    represented = 100 * c(
      metric_value(m, "fraction_copies_in_cross_array_exact_haplotypes"),
      metric_value(m, "fraction_copies_in_cross_array_variant_clusters"),
      metric_value(m, "fraction_copies_in_cross_species_exact_haplotypes"),
      metric_value(m, "fraction_copies_in_cross_species_variant_clusters")
    )
  )
  ggplot(d, aes(scope, represented, fill = definition)) +
    geom_col(position = position_dodge(width = 0.76), width = 0.68, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = sprintf("%.2f%%", represented)),
      position = position_dodge(width = 0.76),
      vjust = -0.35,
      size = 2.7
    ) +
    scale_fill_manual(values = setNames(c("#8A8A8A", "#56B4E9"), c("Exact haplotype", cluster_label))) +
    scale_y_continuous(
      labels = label_continuous_2dp,
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      x = NULL,
      y = "Copies represented (%)",
      title = "Exact versus similarity-based sharing"
    ) +
    theme_sat113() +
    theme(legend.position = "top")
}

region_pair_metric_data <- function(dat, main_threshold) {
  d <- dat$pairs[as_flag(dat$pairs$same_species), , drop = FALSE]
  required <- c(
    "consensus_divergence_pct",
    "copy_pool_dxy_pct",
    "allele_frequency_mean_jsd",
    "variant_composition_jsd"
  )
  missing <- setdiff(required, names(d))
  if (length(missing) > 0L) {
    stop(
      "array_pair_similarity.tsv is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_metric <- sprintf("%d%%-cluster composition", round(100 * main_threshold))
  specifications <- data.frame(
    column = required,
    family = c("Divergence (%)", "Divergence (%)", "JSD (0-1)", "JSD (0-1)"),
    metric = c("Region consensus", "Copy-pool Dxy", "Site-frequency JSD", cluster_metric),
    stringsAsFactors = FALSE
  )
  long <- do.call(rbind, lapply(seq_len(nrow(specifications)), function(i) {
    data.frame(
      family = specifications$family[[i]],
      metric = specifications$metric[[i]],
      value = suppressWarnings(as.numeric(d[[specifications$column[[i]]]])),
      stringsAsFactors = FALSE
    )
  }))
  long <- long[is.finite(long$value), , drop = FALSE]
  long$family <- factor(long$family, levels = c("Divergence (%)", "JSD (0-1)"))
  long$metric <- factor(
    long$metric,
    levels = c("Region consensus", "Copy-pool Dxy", "Site-frequency JSD", cluster_metric)
  )
  long
}

summarize_region_pair_metrics <- function(dat, main_threshold) {
  long <- region_pair_metric_data(dat, main_threshold)
  groups <- split(
    long,
    interaction(long$family, long$metric, drop = TRUE, lex.order = TRUE)
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- x$value[is.finite(x$value)]
    data.frame(
      family = as.character(x$family[[1L]]),
      metric = as.character(x$metric[[1L]]),
      n_region_pairs = length(values),
      q25 = unname(stats::quantile(values, 0.25, na.rm = TRUE)),
      median = stats::median(values, na.rm = TRUE),
      q75 = unname(stats::quantile(values, 0.75, na.rm = TRUE)),
      fraction_approximately_one = if (identical(
        as.character(x$metric[[1L]]),
        sprintf("%d%%-cluster composition", round(100 * main_threshold))
      )) {
        mean(values >= 1 - 1e-9)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

plot_consensus_vs_jsd <- function(dat, main_threshold) {
  d <- dat$pairs[as_flag(dat$pairs$same_species), , drop = FALSE]
  if (nrow(d) == 0L) {
    return(empty_plot("Region-pair turnover", "No same-species region pairs"))
  }

  long <- region_pair_metric_data(dat, main_threshold)
  summary <- summarize_region_pair_metrics(dat, main_threshold)
  summary$family <- factor(summary$family, levels = levels(long$family))
  summary$metric <- factor(summary$metric, levels = levels(long$metric))
  summary$median_label <- ifelse(
    summary$family == "Divergence (%)",
    sprintf("%.2f%%", summary$median),
    sprintf("%.2f", summary$median)
  )
  cluster_metric <- sprintf("%d%%-cluster composition", round(100 * main_threshold))
  cluster_fraction <- summary$fraction_approximately_one[
    as.character(summary$metric) == cluster_metric
  ]
  cluster_fraction <- if (length(cluster_fraction) == 0L) NA_real_ else cluster_fraction[[1L]]
  n_pairs <- sum(as_flag(dat$pairs$same_species))

  metric_palette <- c(
    "Region consensus" = "#7A7A7A",
    "Copy-pool Dxy" = "#D55E00",
    "Site-frequency JSD" = "#0072B2",
    setNames("#CC79A7", cluster_metric)
  )
  boundary_annotation <- if (is.finite(cluster_fraction)) {
    data.frame(
      family = factor("JSD (0-1)", levels = levels(long$family)),
      x = Inf,
      y = 0.12,
      label = sprintf("%.2f%% at JSD approximately 1", 100 * cluster_fraction),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      family = factor(character(), levels = levels(long$family)),
      x = numeric(),
      y = numeric(),
      label = character(),
      stringsAsFactors = FALSE
    )
  }

  ggplot(long, aes(value, color = metric, group = metric)) +
    stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.78) +
    geom_vline(
      data = summary,
      aes(xintercept = median, color = metric),
      linetype = 3,
      linewidth = 0.42,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = summary,
      aes(x = median, y = 0.50, label = median_label, color = metric),
      fill = "white",
      label.size = 0.16,
      size = 2.35,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = boundary_annotation,
      aes(x = x, y = y, label = label),
      hjust = 1.04,
      color = "#8B3F72",
      fill = "white",
      label.size = 0.16,
      size = 2.25,
      inherit.aes = FALSE
    ) +
    facet_wrap(vars(family), scales = "free_x", nrow = 1) +
    scale_color_manual(values = metric_palette, drop = FALSE) +
    scale_x_continuous(
      labels = label_continuous_2dp,
      expand = expansion(mult = c(0.03, 0.14))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = label_percent_2dp,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Metric value",
      y = "Cumulative fraction of region pairs",
      title = "Consensus summaries mask copy-pool turnover",
      subtitle = sprintf("Within-species region pairs (n = %s); labels mark medians", scales::comma(n_pairs))
    ) +
    theme_sat113(base_size = 8.7) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 6.9),
      legend.key.width = grid::unit(0.85, "lines"),
      legend.spacing.x = grid::unit(0.18, "lines"),
      strip.text = element_text(size = 7.5),
      plot.subtitle = element_text(size = 7.3)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

plot_linkage <- function(dat) {
  d <- dat$linkage
  d$NMI <- suppressWarnings(as.numeric(d$categorical_normalized_mutual_information))
  d <- d[is.finite(d$NMI), , drop = FALSE]
  if (nrow(d) == 0L) return(empty_plot("Within-monomer linkage", "No finite site-pair NMI values"))
  wanted <- c("ALL", "FCA", "PPA", "PTI")
  d$group <- factor(d$group, levels = c(intersect(wanted, unique(d$group)), setdiff(unique(d$group), wanted)))
  group_colors <- c(ALL = all_group_color, species_scale_values(setdiff(as.character(d$group), "ALL")))
  ggplot(d, aes(group, NMI, fill = group, color = group)) +
    geom_violin(
      linewidth = 0.35,
      trim = TRUE,
      scale = "width",
      alpha = species_geom_alpha
    ) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.32) +
    geom_hline(yintercept = 0.2, linetype = 2, linewidth = 0.45, color = "#555555") +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(
      x = NULL,
      y = "Normalized mutual information",
      title = "Mutation co-occurrence within monomers"
    ) +
    theme_sat113() +
    theme(legend.position = "none")
}

plot_spatial <- function(dat) {
  d <- dat$arrays
  keep <- is.finite(as.numeric(d$adjacent_same_expected_random_order)) &
    is.finite(as.numeric(d$adjacent_same_observed))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0L) return(empty_plot("Spatial organization", "No spatially testable regions"))
  d$point_size <- sqrt(pmax(as.numeric(d$n_copies), 1))
  upper <- max(c(d$adjacent_same_expected_random_order, d$adjacent_same_observed), na.rm = TRUE) * 1.04
  upper <- max(upper, 0.05)
  colors <- species_scale_values(d$species)
  ggplot(d, aes(adjacent_same_expected_random_order, adjacent_same_observed, color = species)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, color = "#444444") +
    geom_point(aes(size = point_size), alpha = species_geom_alpha, stroke = 0.20) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(1.1, 4.5), guide = "none") +
    scale_x_continuous(labels = label_continuous_2dp) +
    scale_y_continuous(labels = label_continuous_2dp) +
    coord_equal(xlim = c(0, upper), ylim = c(0, upper), expand = FALSE) +
    labs(
      x = "Expected adjacent-same fraction after shuffling",
      y = "Observed adjacent-same fraction",
      title = "Spatial blocks versus interspersion"
    ) +
    theme_sat113() +
    theme(legend.position = "top")
}

plot_array_diversity <- function(dat) {
  d <- dat$arrays
  keep <- is.finite(as.numeric(d$within_array_pi_pct)) &
    is.finite(as.numeric(d$fraction_copies_in_cross_species_clusters))
  d <- d[keep, , drop = FALSE]
  colors <- species_scale_values(d$species)
  ggplot(d, aes(within_array_pi_pct, fraction_copies_in_cross_species_clusters, color = species)) +
    geom_point(size = 1.7, alpha = species_geom_alpha) +
    scale_color_manual(values = colors) +
    scale_x_continuous(labels = label_continuous_2dp) +
    scale_y_continuous(labels = label_percent_2dp) +
    labs(
      x = "Within-region copy diversity, pi (%)",
      y = "Copies in cross-species clusters",
      title = "Region diversity and the shared variant pool"
    ) +
    theme_sat113() +
    theme(legend.position = "top")
}

plot_haplotype_pca <- function(dat) {
  d <- dat$pca
  if ("n_copies" %in% names(d)) d <- d[as.numeric(d$n_copies) >= 2, , drop = FALSE]
  if (nrow(d) > 6000L) {
    d <- d[order(as.numeric(d$n_copies), decreasing = TRUE)[seq_len(6000L)], , drop = FALSE]
  }
  if (nrow(d) == 0L) return(empty_plot("Exact monomer haplotypes", "No haplotypes pass the display filter"))
  d$size_value <- log2(pmax(as.numeric(d$n_copies), 1) + 1)
  shared <- d[d$sharing_class == "shared_exact_haplotype", , drop = FALSE]
  colors <- species_scale_values(d$dominant_species)
  pc1_ratio <- suppressWarnings(as.numeric(d$PC1_explained_variance_ratio[[1L]]))
  pc2_ratio <- suppressWarnings(as.numeric(d$PC2_explained_variance_ratio[[1L]]))
  xlab <- if (is.finite(pc1_ratio)) sprintf("PC1 (%.2f%%)", 100 * pc1_ratio) else "PC1"
  ylab <- if (is.finite(pc2_ratio)) sprintf("PC2 (%.2f%%)", 100 * pc2_ratio) else "PC2"

  ggplot(d, aes(PC1, PC2)) +
    geom_point(
      aes(color = dominant_species, size = size_value),
      alpha = species_geom_alpha,
      stroke = 0
    ) +
    geom_point(
      data = shared,
      aes(fill = dominant_species, size = size_value),
      shape = 21,
      color = "black",
      alpha = species_geom_alpha,
      stroke = 0.32
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_size_continuous(range = c(0.45, 3.4), guide = "none") +
    scale_x_continuous(labels = label_continuous_2dp) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = xlab, y = ylab, title = "Reference-phased exact monomer haplotypes") +
    theme_sat113() +
    theme(legend.position = "top")
}

plot_cluster_origin <- function(dat, main_threshold) {
  d <- aggregate(n_copies ~ classification, data = dat$clusters, FUN = sum)
  d <- d[order(d$n_copies), , drop = FALSE]
  d$classification <- factor(d$classification, levels = d$classification)
  d$display <- pretty_class(as.character(d$classification))
  d$display <- factor(d$display, levels = d$display)
  d$fill_group <- ifelse(
    startsWith(as.character(d$classification), "shared"),
    "Shared",
    ifelse(endsWith(as.character(d$classification), "specific"), "Species-specific", "Other")
  )
  ggplot(d, aes(n_copies, display, fill = fill_group)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(Shared = "#56B4E9", `Species-specific` = "#E69F00", Other = "#8A8A8A")) +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(
      x = "Copies",
      y = NULL,
      title = sprintf("Variant-cluster origin classes (%d%%)", round(100 * main_threshold))
    ) +
    theme_sat113() +
    theme(legend.position = "top")
}

species_pair_long <- function(dat, analysis = NULL) {
  d <- dat$species_pairs
  rows <- lapply(seq_len(nrow(d)), function(i) {
    pair <- paste0(d$species_1[[i]], "-", d$species_2[[i]])
    data.frame(
      pair = pair,
      metric = c("Consensus divergence", "Within-species pi", "Between-species Dxy", "Net divergence Da"),
      value = c(
        as.numeric(d$consensus_divergence_pct[[i]]),
        as.numeric(d$mean_within_species_pi_pct[[i]]),
        as.numeric(d$copy_pool_dxy_pct[[i]]),
        as.numeric(d$net_divergence_da_pct[[i]])
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (!is.null(analysis)) out$analysis <- analysis
  out
}

plot_species_metrics <- function(dat) {
  d <- species_pair_long(dat)
  d$metric <- factor(d$metric, levels = names(metric_colors))
  ggplot(d, aes(pair, value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.72) +
    scale_fill_manual(values = metric_colors, drop = FALSE) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(
      x = NULL,
      y = "Substitution divergence (%)",
      title = "Consensus versus full copy pools"
    ) +
    theme_sat113() +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "top",
      legend.text = element_text(size = 7.5)
    )
}

plot_rarefaction <- function(dat, main_threshold) {
  d <- dat$rarefaction[
    abs(as.numeric(dat$rarefaction$identity_threshold) - main_threshold) < 1e-9,
    , drop = FALSE
  ]
  cols <- grep("^jaccard_", names(d), value = TRUE)
  if (nrow(d) == 0L || length(cols) == 0L) {
    return(empty_plot("Rarefied sharing", "No rarefaction results at the requested threshold"))
  }
  long <- do.call(rbind, lapply(cols, function(column) {
    data.frame(
      replicate = d$replicate,
      species_pair = gsub("_", "-", sub("^jaccard_", "", column)),
      Jaccard = as.numeric(d[[column]]),
      stringsAsFactors = FALSE
    )
  }))
  ggplot(long, aes(species_pair, Jaccard)) +
    geom_boxplot(width = 0.55, fill = "#009E73", color = "#145B45", outlier.shape = NA, alpha = 0.78) +
    geom_jitter(width = 0.14, height = 0, size = 0.65, alpha = 0.16, color = "black") +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(
      x = NULL,
      y = "Cluster-presence Jaccard",
      title = "Rarefied cross-species sharing"
    ) +
    theme_sat113() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

top_n_per_group <- function(d, group, value, n) {
  split_d <- split(d, d[[group]])
  out <- lapply(split_d, function(x) {
    x <- x[order(x[[value]], decreasing = TRUE), , drop = FALSE]
    x[seq_len(min(n, nrow(x))), , drop = FALSE]
  })
  do.call(rbind, out)
}

plot_top_variant_heatmap <- function(dat, top_variants = 40L, top_regions_per_species = 25L) {
  clusters <- dat$clusters[order(as.numeric(dat$clusters$n_copies), decreasing = TRUE), , drop = FALSE]
  top_ids <- head(as.character(clusters$variant_cluster), top_variants)

  comp <- dat$variant_composition
  totals <- aggregate(n_copies ~ array_uid + species, data = comp, FUN = sum)
  selected_regions <- top_n_per_group(
    totals, group = "species", value = "n_copies", n = top_regions_per_species
  )
  species_order <- intersect(c("FCA", "PPA", "PTI"), unique(selected_regions$species))
  selected_regions$species_rank <- match(selected_regions$species, species_order)
  selected_regions <- selected_regions[order(selected_regions$species_rank, -selected_regions$n_copies), , drop = FALSE]
  region_ids <- as.character(selected_regions$array_uid)

  subcomp <- comp[
    comp$variant_cluster %in% top_ids & comp$array_uid %in% region_ids,
    c("array_uid", "species", "variant_cluster", "fraction"),
    drop = FALSE
  ]
  grid <- expand.grid(
    array_uid = region_ids,
    variant_cluster = top_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- merge(grid, subcomp, by = c("array_uid", "variant_cluster"), all.x = TRUE, sort = FALSE)
  grid$fraction[is.na(grid$fraction)] <- 0
  region_meta <- selected_regions[, c("array_uid", "species"), drop = FALSE]
  grid <- merge(grid, region_meta, by = "array_uid", all.x = TRUE, suffixes = c("", ".meta"), sort = FALSE)
  if ("species.meta" %in% names(grid)) {
    grid$species <- ifelse(is.na(grid$species), grid$species.meta, grid$species)
  }
  grid$array_uid <- factor(grid$array_uid, levels = region_ids)
  grid$variant_cluster <- factor(grid$variant_cluster, levels = rev(top_ids))
  grid$log_fraction <- log10(as.numeric(grid$fraction) + 1e-4)

  ggplot(grid, aes(array_uid, variant_cluster, fill = log_fraction)) +
    geom_tile(linewidth = 0.05, color = "#151515") +
    facet_grid(cols = vars(species), scales = "free_x", space = "free_x") +
    scale_fill_gradientn(
      colors = HEATMAP_COLORS,
      limits = c(-4, 0.00005),
      breaks = c(-4, -3, -2, -1, 0),
      labels = c("0", "0.1%", "1%", "10%", "100%"),
      na.value = HEATMAP_NA_COLOR,
      name = "Region fraction"
    ) +
    scale_x_discrete(labels = function(x) gsub("__", ":", sub("^[^|]+\\|", "", x), fixed = TRUE)) +
    labs(
      x = sprintf("Representative %s-rich region", satellite_label),
      y = sprintf("Top %d variant clusters", length(top_ids)),
      title = sprintf("Abundant %s variant clusters across regions", satellite_label)
    ) +
    theme_minimal(base_size = 8.5, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5, color = "black"),
      axis.text.y = element_text(size = 6.5, color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = 10.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      legend.position = "right"
    )
}

component_number <- function(x) {
  suppressWarnings(as.integer(sub(".*_", "", as.character(x))))
}

plot_model_bic <- function(dat, title = "Model complexity") {
  m <- dat$subtypes$model
  if (is.null(m)) return(empty_plot(title, "No model-selection table"))
  selected <- m[as.numeric(m$k) == dat$recommended_k, , drop = FALSE]
  ggplot(m, aes(as.numeric(k), as.numeric(delta_bic_vs_k1))) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.40, color = "#555555") +
    geom_line(color = "#0072B2", linewidth = 0.78) +
    geom_point(color = "#0072B2", size = 1.8) +
    geom_vline(xintercept = dat$recommended_k, linetype = 3, linewidth = 0.65, color = "#D55E00") +
    geom_point(data = selected, shape = 21, size = 3.1, fill = "#D55E00", color = "black", stroke = 0.35) +
    scale_x_continuous(breaks = sort(unique(as.numeric(m$k)))) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = "Number of categorical components (K)", y = "Delta BIC vs K=1", title = title) +
    theme_sat113()
}

plot_model_cv <- function(dat, title = "Generalization to held-out regions") {
  m <- dat$subtypes$model
  if (is.null(m)) return(empty_plot(title, "No model-selection table"))
  selected <- m[as.numeric(m$k) == dat$recommended_k, , drop = FALSE]
  ggplot(m, aes(as.numeric(k), as.numeric(cv_loglik_per_monomer_mean))) +
    geom_errorbar(
      aes(
        ymin = as.numeric(cv_loglik_per_monomer_mean) - as.numeric(cv_loglik_per_monomer_se),
        ymax = as.numeric(cv_loglik_per_monomer_mean) + as.numeric(cv_loglik_per_monomer_se)
      ),
      width = 0.16,
      linewidth = 0.40,
      color = "#007A5E"
    ) +
    geom_line(color = "#009E73", linewidth = 0.78) +
    geom_point(color = "#009E73", size = 1.8) +
    geom_vline(xintercept = dat$recommended_k, linetype = 3, linewidth = 0.65, color = "#D55E00") +
    geom_point(data = selected, shape = 21, size = 3.1, fill = "#D55E00", color = "black", stroke = 0.35) +
    scale_x_continuous(breaks = sort(unique(as.numeric(m$k)))) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = "K", y = "Held-out log likelihood per monomer", title = title) +
    theme_sat113()
}

plot_component_species <- function(dat) {
  d <- dat$subtypes$summary
  if (is.null(d)) return(empty_plot("Species composition", "No component summary"))
  species_cols <- grep("^n_", names(d), value = TRUE)
  species_cols <- setdiff(species_cols, c("n_copies", "n_haplotypes", "n_units"))
  long <- do.call(rbind, lapply(species_cols, function(column) {
    data.frame(
      component = component_number(d$subtype),
      species = sub("^n_", "", column),
      n_copies = as.numeric(d[[column]]),
      stringsAsFactors = FALSE
    )
  }))
  long$component <- factor(long$component, levels = sort(unique(long$component)))
  colors <- species_scale_values(long$species)
  ggplot(long, aes(component, n_copies, fill = species)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
      x = "Model-derived component",
      y = "Hard-assigned copies",
      title = "Species composition of model-derived components"
    ) +
    theme_sat113() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
}

plot_diagnostic_sites <- function(dat, top_n = 25L) {
  d <- dat$subtypes$diagnostics
  k <- dat$recommended_k
  if (is.null(d) || !is.finite(k)) return(empty_plot("Multi-site component pattern", "No diagnostic-site table"))
  d <- head(d, min(top_n, nrow(d)))
  rows <- vector("list", k)
  for (component in seq_len(k)) {
    freq_col <- paste0("subtype_", component, "_consensus_frequency")
    base_col <- paste0("subtype_", component, "_consensus_base")
    if (!all(c(freq_col, base_col) %in% names(d))) next
    rows[[component]] <- data.frame(
      component = component,
      position = as.character(d$alignment_position_1based),
      frequency = as.numeric(d[[freq_col]]),
      base = as.character(d[[base_col]]),
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  long$position <- factor(long$position, levels = as.character(d$alignment_position_1based))
  long$component <- factor(long$component, levels = rev(seq_len(k)), labels = paste("Component", rev(seq_len(k))))
  long$text_color <- ifelse(long$frequency >= 0.82, "white", "black")
  ggplot(long, aes(position, component, fill = frequency)) +
    geom_tile(color = "white", linewidth = 0.18) +
    geom_text(aes(label = base, color = text_color), size = 2.35, fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = HEATMAP_COLORS,
      limits = c(0.25, 1),
      oob = scales::squish,
      na.value = HEATMAP_NA_COLOR,
      labels = label_continuous_2dp,
      name = "Major-base frequency"
    ) +
    labs(
      x = "Top diagnostic alignment positions",
      y = NULL,
      title = "Multi-site frequency pattern"
    ) +
    theme_minimal(base_size = 8.5, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6.5, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      plot.title = element_text(face = "bold", size = 10.3),
      legend.position = "right"
    )
}

region_order <- function(composition) {
  meta <- unique(composition[, c("unit_uid", "species", "chrom", "unit_key"), drop = FALSE])
  preferred <- c("FCA", "PPA", "PTI")
  meta$species_rank <- match(meta$species, c(intersect(preferred, meta$species), setdiff(unique(meta$species), preferred)))
  meta$chrom_prefix <- sub("[0-9]+.*$", "", meta$chrom)
  chrom_num_text <- sub("^[^0-9]*([0-9]+).*$", "\\1", meta$chrom)
  meta$chrom_num <- suppressWarnings(as.numeric(chrom_num_text))
  meta$chrom_num[!is.finite(meta$chrom_num)] <- Inf
  meta <- meta[order(meta$species_rank, meta$chrom_prefix, meta$chrom_num, meta$chrom, meta$unit_key), , drop = FALSE]
  meta
}

plot_component_array_heatmap <- function(dat, title = NULL, show_y_labels = TRUE) {
  d <- dat$subtypes$composition
  k <- dat$recommended_k
  if (is.null(d) || !is.finite(k)) return(empty_plot(title %||% "Region composition", "No component-composition table"))
  meta <- region_order(d)
  components <- sort(unique(component_number(d$subtype)))
  grid <- expand.grid(
    unit_uid = meta$unit_uid,
    component = components,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  d2 <- data.frame(
    unit_uid = d$unit_uid,
    component = component_number(d$subtype),
    fraction = as.numeric(d$fraction),
    stringsAsFactors = FALSE
  )
  grid <- merge(grid, d2, by = c("unit_uid", "component"), all.x = TRUE, sort = FALSE)
  grid$fraction[is.na(grid$fraction)] <- 0
  grid <- merge(grid, meta[, c("unit_uid", "species"), drop = FALSE], by = "unit_uid", all.x = TRUE, sort = FALSE)
  grid$unit_uid <- factor(grid$unit_uid, levels = rev(meta$unit_uid))
  grid$component <- factor(grid$component, levels = components, labels = paste0("C", components))
  y_labels <- if (show_y_labels) {
    function(x) gsub("__", ":", sub("^[^|]+\\|", "", x), fixed = TRUE)
  } else {
    NULL
  }
  ggplot(grid, aes(component, unit_uid, fill = fraction)) +
    geom_tile(color = "white", linewidth = 0.10) +
    facet_grid(rows = vars(species), scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors = HEATMAP_COLORS,
      limits = c(0, 1),
      na.value = HEATMAP_NA_COLOR,
      labels = label_continuous_2dp,
      name = "Region fraction"
    ) +
    scale_y_discrete(labels = y_labels) +
    labs(
      x = sprintf("Model-derived component (diagnostic K=%d)", k),
      y = sprintf("%s-rich region", satellite_label),
      title = title %||% sprintf("Region-level composition (K=%d)", k)
    ) +
    theme_minimal(base_size = 7.5, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 6.5),
      axis.text.y = element_text(color = "black", size = if (show_y_labels) 4.5 else 0),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      strip.text.y = element_text(face = "bold", angle = 0),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      legend.position = "right"
    )
}

analysis_summary <- function(dat, label) {
  p_cluster <- suppressWarnings(as.numeric(dat$arrays$permutation_p_spatial_clustering))
  p_alt <- suppressWarnings(as.numeric(dat$arrays$permutation_p_spatial_alternation))
  q_cluster <- rep(NA_real_, length(p_cluster))
  q_alt <- rep(NA_real_, length(p_alt))
  q_cluster[is.finite(p_cluster)] <- p.adjust(p_cluster[is.finite(p_cluster)], method = "BH")
  q_alt[is.finite(p_alt)] <- p.adjust(p_alt[is.finite(p_alt)], method = "BH")
  n_spatial_q <- sum((q_cluster < 0.05) | (q_alt < 0.05), na.rm = TRUE)
  data.frame(
    analysis = label,
    n_copies = metric_value(dat$metrics, "n_copies"),
    n_core_positions = sum(as_flag(dat$core_columns$is_core_column), na.rm = TRUE),
    cross_species_exact = metric_value(dat$metrics, "fraction_copies_in_cross_species_exact_haplotypes"),
    cross_species_cluster = metric_value(dat$metrics, "fraction_copies_in_cross_species_variant_clusters"),
    median_nmi = metric_value(dat$metrics, "global_median_site_pair_nmi"),
    q95_nmi = metric_value(dat$metrics, "global_q95_site_pair_nmi"),
    n_spatial_q_lt_0_05 = n_spatial_q,
    stringsAsFactors = FALSE
  )
}

plot_sensitivity_inputs <- function(primary, sensitivity) {
  d <- rbind(
    analysis_summary(primary, cfg$primary_label),
    analysis_summary(sensitivity, cfg$sensitivity_label)
  )
  long <- rbind(
    data.frame(analysis = d$analysis, metric = "Retained monomer copies", value = d$n_copies),
    data.frame(analysis = d$analysis, metric = "Retained core positions", value = d$n_core_positions)
  )
  long$analysis <- factor(long$analysis, levels = c(cfg$primary_label, cfg$sensitivity_label))
  ggplot(long, aes(analysis, value, fill = analysis)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = label_count_integer(value)), vjust = -0.35, size = 2.8) +
    facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
    scale_fill_manual(values = analysis_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14)), labels = scales::label_comma()) +
    labs(x = NULL, y = NULL, title = "Retained data after filtering") +
    theme_sat113() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
}

plot_sensitivity_sharing <- function(primary, sensitivity) {
  summary <- rbind(
    analysis_summary(primary, cfg$primary_label),
    analysis_summary(sensitivity, cfg$sensitivity_label)
  )
  d <- rbind(
    data.frame(analysis = summary$analysis, definition = "Exact haplotype", value = summary$cross_species_exact),
    data.frame(
      analysis = summary$analysis,
      definition = sprintf("Variant cluster (%d%%)", round(100 * cfg$main_threshold)),
      value = summary$cross_species_cluster
    )
  )
  d$analysis <- factor(d$analysis, levels = c(cfg$primary_label, cfg$sensitivity_label))
  ggplot(d, aes(definition, value, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.74), width = 0.66) +
    geom_text(
      aes(label = scales::percent(value, accuracy = 0.01)),
      position = position_dodge(width = 0.74),
      vjust = -0.35,
      size = 2.7
    ) +
    scale_fill_manual(values = analysis_colors) +
    scale_y_continuous(labels = label_percent_2dp, expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Cross-species copies", title = "Cross-species sharing") +
    theme_sat113() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))
}

plot_sensitivity_divergence <- function(primary, sensitivity) {
  p <- species_pair_long(primary, cfg$primary_label)
  s <- species_pair_long(sensitivity, cfg$sensitivity_label)
  d <- rbind(p, s)
  d <- d[d$metric %in% c("Within-species pi", "Between-species Dxy", "Net divergence Da"), , drop = FALSE]
  d$metric <- factor(
    d$metric,
    levels = c("Within-species pi", "Between-species Dxy", "Net divergence Da"),
    labels = c("Mean within-species pi", "Between-species Dxy", "Net divergence Da")
  )
  d$analysis <- factor(d$analysis, levels = c(cfg$primary_label, cfg$sensitivity_label))
  ggplot(d, aes(pair, value, color = analysis, group = analysis)) +
    geom_point(position = position_dodge(width = 0.32), size = 2.0) +
    facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
    scale_color_manual(values = analysis_colors) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = NULL, y = "Divergence (%)", title = "Copy-pool divergence") +
    theme_sat113() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1))
}

plot_sensitivity_linkage <- function(primary, sensitivity) {
  summary <- rbind(
    analysis_summary(primary, cfg$primary_label),
    analysis_summary(sensitivity, cfg$sensitivity_label)
  )
  d <- rbind(
    data.frame(analysis = summary$analysis, statistic = "Median", NMI = summary$median_nmi),
    data.frame(analysis = summary$analysis, statistic = "95th percentile", NMI = summary$q95_nmi)
  )
  d$analysis <- factor(d$analysis, levels = c(cfg$primary_label, cfg$sensitivity_label))
  spatial_note <- paste0(
    "BH-significant spatial regions: ",
    paste(paste0(summary$analysis, "=", summary$n_spatial_q_lt_0_05), collapse = "; ")
  )
  ggplot(d, aes(statistic, NMI, color = analysis, group = analysis)) +
    geom_line(linewidth = 0.60, position = position_dodge(width = 0.18)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.18)) +
    scale_color_manual(values = analysis_colors) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(
      x = NULL,
      y = "Site-pair NMI",
      title = "Weak linkage is insensitive to filtering",
      subtitle = spatial_note
    ) +
    theme_sat113() +
    theme(legend.position = "top")
}

choose2_sum <- function(x) sum(x * (x - 1) / 2)

adjusted_rand_index <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- as.character(x[keep])
  y <- as.character(y[keep])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  tab <- table(x, y)
  nij <- choose2_sum(as.numeric(tab))
  ai <- choose2_sum(rowSums(tab))
  bj <- choose2_sum(colSums(tab))
  total <- n * (n - 1) / 2
  expected <- ai * bj / total
  maximum <- 0.5 * (ai + bj)
  if (abs(maximum - expected) < .Machine$double.eps) return(1)
  (nij - expected) / (maximum - expected)
}

component_assignment_agreement <- function(primary, sensitivity) {
  p <- primary$subtypes$copy_assignments
  s <- sensitivity$subtypes$copy_assignments
  if (is.null(p) || is.null(s)) return(list(ari = NA_real_, n = 0L))
  key_columns <- c("species", "chrom", "start", "end", "strand")
  if (!all(key_columns %in% names(p)) || !all(key_columns %in% names(s))) {
    key_columns <- intersect(c("species", "copy_id"), intersect(names(p), names(s)))
  }
  if (length(key_columns) == 0L) return(list(ari = NA_real_, n = 0L))
  p$key <- do.call(paste, c(p[key_columns], sep = "|"))
  s$key <- do.call(paste, c(s[key_columns], sep = "|"))
  merged <- merge(
    p[, c("key", "subtype")],
    s[, c("key", "subtype")],
    by = "key",
    suffixes = c(".primary", ".sensitivity")
  )
  list(
    ari = adjusted_rand_index(merged$subtype.primary, merged$subtype.sensitivity),
    n = nrow(merged)
  )
}

combined_model_table <- function(primary, sensitivity) {
  p <- primary$subtypes$model
  s <- sensitivity$subtypes$model
  if (is.null(p) || is.null(s)) return(NULL)
  p$analysis <- cfg$primary_label
  s$analysis <- cfg$sensitivity_label
  out <- rbind(p, s)
  out$analysis <- factor(out$analysis, levels = c(cfg$primary_label, cfg$sensitivity_label))
  out
}

selected_k_table <- function(primary, sensitivity) {
  data.frame(
    analysis = factor(
      c(cfg$primary_label, cfg$sensitivity_label),
      levels = c(cfg$primary_label, cfg$sensitivity_label)
    ),
    k = c(primary$recommended_k, sensitivity$recommended_k)
  )
}

plot_model_comparison <- function(primary, sensitivity) {
  d <- combined_model_table(primary, sensitivity)
  if (is.null(d)) return(empty_plot("Model sensitivity", "Model tables are unavailable"))
  selected <- selected_k_table(primary, sensitivity)
  agreement <- component_assignment_agreement(primary, sensitivity)
  caption <- if (is.finite(agreement$ari)) {
    sprintf(
      "Agreement between selected solutions on %s shared copies: ARI = %.2f. Components are descriptive model-derived patterns, not confirmed biological subtypes.",
      scales::comma(agreement$n), agreement$ari
    )
  } else {
    "Components are descriptive model-derived patterns, not confirmed biological subtypes."
  }

  p_bic <- ggplot(d, aes(as.numeric(k), as.numeric(delta_bic_vs_k1), color = analysis)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, color = "#555555") +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.55) +
    geom_vline(data = selected, aes(xintercept = k, color = analysis), linetype = 3, linewidth = 0.52, show.legend = FALSE) +
    scale_color_manual(values = analysis_colors) +
    scale_x_continuous(breaks = sort(unique(as.numeric(d$k)))) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = "K", y = "Delta BIC vs K=1", title = "Model complexity") +
    theme_sat113() +
    theme(legend.position = "top")

  p_cv <- ggplot(d, aes(as.numeric(k), as.numeric(cv_loglik_per_monomer_mean), color = analysis)) +
    geom_errorbar(
      aes(
        ymin = as.numeric(cv_loglik_per_monomer_mean) - as.numeric(cv_loglik_per_monomer_se),
        ymax = as.numeric(cv_loglik_per_monomer_mean) + as.numeric(cv_loglik_per_monomer_se)
      ),
      width = 0.12,
      linewidth = 0.30
    ) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.55) +
    geom_vline(data = selected, aes(xintercept = k, color = analysis), linetype = 3, linewidth = 0.52, show.legend = FALSE) +
    scale_color_manual(values = analysis_colors) +
    scale_x_continuous(breaks = sort(unique(as.numeric(d$k)))) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = "K", y = "Held-out log likelihood per monomer", title = "Held-out-region cross-validation") +
    theme_sat113() +
    theme(legend.position = "top")

  p_stability <- ggplot(d, aes(as.numeric(k), as.numeric(stability_mean_ari), color = analysis)) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.55) +
    geom_vline(data = selected, aes(xintercept = k, color = analysis), linetype = 3, linewidth = 0.52, show.legend = FALSE) +
    scale_color_manual(values = analysis_colors) +
    scale_x_continuous(breaks = sort(unique(as.numeric(d$k)))) +
    scale_y_continuous(labels = label_continuous_2dp) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "K", y = "Mean ARI among initializations", title = "Initialization stability") +
    theme_sat113() +
    theme(legend.position = "top")

  thresholds <- data.frame(
    analysis = factor(
      c(cfg$primary_label, cfg$sensitivity_label),
      levels = c(cfg$primary_label, cfg$sensitivity_label)
    ),
    threshold = c(
      suppressWarnings(as.numeric(conclusion_value(
        primary$subtypes$conclusion, "minimum_component_weight_threshold", NA_character_
      ))),
      suppressWarnings(as.numeric(conclusion_value(
        sensitivity$subtypes$conclusion, "minimum_component_weight_threshold", NA_character_
      )))
    )
  )
  p_component <- ggplot(d, aes(as.numeric(k), as.numeric(min_hard_component_weight), color = analysis)) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.55) +
    geom_hline(
      data = thresholds[is.finite(thresholds$threshold), , drop = FALSE],
      aes(yintercept = threshold, color = analysis),
      linetype = 2,
      linewidth = 0.50,
      show.legend = FALSE
    ) +
    geom_vline(data = selected, aes(xintercept = k, color = analysis), linetype = 3, linewidth = 0.52, show.legend = FALSE) +
    scale_color_manual(values = analysis_colors) +
    scale_x_continuous(breaks = sort(unique(as.numeric(d$k)))) +
    scale_y_continuous(labels = label_continuous_2dp) +
    labs(x = "K", y = "Smallest component weight", title = "Minimum component-size criterion") +
    theme_sat113() +
    theme(legend.position = "top")

  wrap_plots(p_bic, p_cv, p_stability, p_component, ncol = 2) +
    plot_annotation(
      title = "Sensitivity of the diagnostic categorical-mixture solution",
      caption = caption,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.caption = element_text(size = 8.5, hjust = 0),
        plot.tag = element_text(face = "bold", size = 13)
      )
    )
}

make_legacy_figures <- function(dat) {
  diagnostics <- wrap_plots(
    plot_threshold_sensitivity(dat, cfg$main_threshold),
    plot_exact_vs_cluster(dat, cfg$main_threshold),
    plot_consensus_vs_jsd(dat, cfg$main_threshold),
    plot_linkage(dat),
    plot_spatial(dat),
    plot_array_diversity(dat),
    ncol = 3
  ) +
    plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 13)))
  save_figure(
    diagnostics,
    paste0(satellite_stem, "_variant_diagnostics_R"),
    15.5,
    9.2,
    sprintf(
      "R replacement for the original Python %s_variant_diagnostics figure.",
      satellite_label
    )
  )

  origin <- wrap_plots(
    plot_haplotype_pca(dat),
    plot_cluster_origin(dat, cfg$main_threshold),
    plot_species_metrics(dat),
    plot_rarefaction(dat, cfg$main_threshold),
    ncol = 2
  ) +
    plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 13)))
  save_figure(
    origin,
    paste0(satellite_stem, "_species_origin_R"),
    13.5,
    10.2,
    sprintf(
      "R replacement for the original Python %s_species_origin figure.",
      satellite_label
    )
  )

  variant_heatmap <- plot_top_variant_heatmap(
    dat,
    top_variants = cfg$top_variants,
    top_regions_per_species = cfg$top_regions_per_species
  )
  save_figure(
    variant_heatmap,
    paste0(satellite_stem, "_top_variant_array_heatmap_R"),
    18,
    9,
    "R replacement for the original Python top-variant region heatmap."
  )

  if (!is.null(dat$subtypes$model)) {
    validation <- wrap_plots(
      plot_model_bic(dat),
      plot_model_cv(dat),
      plot_component_species(dat),
      plot_diagnostic_sites(dat),
      ncol = 2
    ) +
      plot_annotation(
        title = sprintf("Diagnostic categorical-mixture analysis (K=%d)", dat$recommended_k),
        subtitle = "Model-derived components are descriptive patterns, not confirmed biological subtypes.",
        tag_levels = "A",
        theme = theme(
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "#555555"),
          plot.tag = element_text(face = "bold", size = 13)
        )
      )
    save_figure(
      validation,
      paste0(satellite_stem, "_model_component_validation_R"),
      13.2,
      9.4,
      "R replacement for the original subtype-validation figure, using model-derived-component terminology."
    )

    component_heatmap <- plot_component_array_heatmap(dat)
    save_figure(
      component_heatmap,
      paste0(satellite_stem, "_model_component_array_composition_R"),
      8.5,
      17,
      "R replacement for the original subtype-array-composition heatmap."
    )
  }
}

make_manuscript_figures <- function(primary, sensitivity = NULL) {
  main_figure <- wrap_plots(
    plot_haplotype_pca(primary),
    plot_species_metrics(primary),
    plot_exact_vs_cluster(primary, cfg$main_threshold),
    plot_consensus_vs_jsd(primary, cfg$main_threshold),
    plot_linkage(primary),
    plot_spatial(primary),
    ncol = 3
  ) +
    plot_annotation(
      title = sprintf(
        "Cross-species conservation and region-level turnover of %s monomer variation",
        satellite_label
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.tag = element_text(face = "bold", size = 13)
      )
    )
  save_figure(
    main_figure,
    paste0("Figure_main_", satellite_stem, "_variation_R"),
    15.8,
    9.4,
    paste0(
      "Recommended six-panel main figure based on the strict primary analysis; ",
      "panel D compares region-pair divergence and composition distributions."
    )
  )

  origin_aux <- wrap_plots(
    plot_cluster_origin(primary, cfg$main_threshold),
    plot_rarefaction(primary, cfg$main_threshold),
    ncol = 2
  ) +
    plot_annotation(
      title = sprintf(
        "Cross-species distribution of reference-phased %s variation",
        satellite_label
      ),
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 12), plot.tag = element_text(face = "bold", size = 13))
    )
  save_figure(
    origin_aux,
    paste0("Figure_S2_", satellite_stem, "_cross_species_auxiliary_R"),
    11.5,
    5.2,
    "Supplementary Figure S2: cluster-origin classes and rarefied cross-species sharing."
  )

  threshold_aux <- wrap_plots(
    plot_threshold_sensitivity(primary, cfg$main_threshold),
    plot_array_diversity(primary),
    ncol = 2
  ) +
    plot_annotation(
      title = "Clustering-threshold and region-diversity diagnostics",
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 12), plot.tag = element_text(face = "bold", size = 13))
    )
  save_figure(
    threshold_aux,
    paste0("Figure_S3_", satellite_stem, "_threshold_and_array_diversity_R"),
    11.5,
    5.2,
    "Supplementary Figure S3: clustering-threshold sensitivity and within-region diversity."
  )

  variant_heatmap <- plot_top_variant_heatmap(
    primary,
    top_variants = cfg$top_variants,
    top_regions_per_species = cfg$top_regions_per_species
  )
  save_figure(
      variant_heatmap,
      paste0("Figure_S4_", satellite_stem, "_top_variant_array_heatmap_R"),
      18,
      9,
      sprintf(
        "Supplementary Figure S4: abundant %d%%-identity variant clusters across representative regions.",
        round(100 * cfg$main_threshold)
      )
  )

  if (!is.null(sensitivity)) {
    sensitivity_figure <- wrap_plots(
      plot_sensitivity_inputs(primary, sensitivity),
      plot_sensitivity_sharing(primary, sensitivity),
      plot_sensitivity_divergence(primary, sensitivity),
      plot_sensitivity_linkage(primary, sensitivity),
      ncol = 2
    ) +
      plot_annotation(
        title = sprintf(
          "Robustness of %s conclusions to copy-length and alignment-edge filtering",
          satellite_label
        ),
        tag_levels = "A",
        theme = theme(plot.title = element_text(face = "bold", size = 12.5), plot.tag = element_text(face = "bold", size = 13))
      )
    save_figure(
      sensitivity_figure,
      paste0("Figure_S1_", satellite_stem, "_parameter_sensitivity_R"),
      13.2,
      9.0,
      "Supplementary Figure S1: primary versus combined filtering-sensitivity analysis."
    )

    model_figure <- plot_model_comparison(primary, sensitivity)
    save_figure(
      model_figure,
      paste0("Figure_S5_", satellite_stem, "_model_sensitivity_R"),
      13.2,
      9.2,
      "Supplementary Figure S5: K, cross-validation, initialization stability, component size, and cross-run ARI."
    )

    p_primary <- plot_component_array_heatmap(
      primary,
      title = sprintf("%s: diagnostic K=%d", cfg$primary_label, primary$recommended_k),
      show_y_labels = TRUE
    )
    p_sensitivity <- plot_component_array_heatmap(
      sensitivity,
      title = sprintf("%s: diagnostic K=%d", cfg$sensitivity_label, sensitivity$recommended_k),
      show_y_labels = FALSE
    )
    component_comparison <- wrap_plots(p_primary, p_sensitivity, ncol = 2, guides = "collect") +
      plot_annotation(
        title = "Filtering sensitivity of model-derived component composition",
        caption = "Component numbers are run-specific labels and should not be matched across panels by number.",
        tag_levels = "A",
        theme = theme(
          plot.title = element_text(face = "bold", size = 12.5),
          plot.caption = element_text(size = 8.5, hjust = 0),
          plot.tag = element_text(face = "bold", size = 13)
        )
      ) & theme(legend.position = "right")
    save_figure(
      component_comparison,
      paste0("Figure_S6_", satellite_stem, "_component_array_composition_R"),
      15.5,
      18,
      sprintf(
        "Supplementary Figure S6: side-by-side K=%d/K=%d region composition heatmaps.",
        primary$recommended_k,
        sensitivity$recommended_k
      )
    )
  } else if (!is.null(primary$subtypes$composition)) {
    component_heatmap <- plot_component_array_heatmap(primary)
    save_figure(
      component_heatmap,
      paste0("Figure_S6_", satellite_stem, "_component_array_composition_R"),
      8.5,
      17,
      "Supplementary Figure S6: model-derived component composition for the primary run only."
    )
  }
}

message("Primary results: ", primary$result_dir)
if (!is.null(sensitivity)) message("Sensitivity results: ", sensitivity$result_dir)
message("Writing R figures to: ", cfg$outdir)

if (cfg$layout %in% c("all", "legacy")) make_legacy_figures(primary)
if (cfg$layout %in% c("all", "manuscript")) make_manuscript_figures(primary, sensitivity)

region_pair_summary_path <- file.path(
  cfg$outdir,
  paste0(satellite_stem, "_region_pair_distribution_summary.tsv")
)
region_pair_summary_out <- summarize_region_pair_metrics(primary, cfg$main_threshold)

if ("n_region_pairs" %in% names(region_pair_summary_out)) {
  region_pair_summary_out$n_region_pairs <- label_count_integer(
    region_pair_summary_out$n_region_pairs
  )
}

continuous_summary_cols <- intersect(
  c("q25", "median", "q75", "fraction_approximately_one"),
  names(region_pair_summary_out)
)
region_pair_summary_out[continuous_summary_cols] <- lapply(
  region_pair_summary_out[continuous_summary_cols],
  label_continuous_2dp
)

write.table(
  region_pair_summary_out,
  region_pair_summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

palette_path <- file.path(cfg$outdir, paste0(satellite_stem, "_species_palette.tsv"))
write.table(
  data.frame(
    species = names(species_colors),
    color_hex = unname(species_colors),
    alpha = formatC(species_geom_alpha, format = "f", digits = 2),
    stringsAsFactors = FALSE
  ),
  palette_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest_path <- file.path(cfg$outdir, paste0(satellite_stem, "_R_figure_manifest.tsv"))
write.table(saved_manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Completed: ", nrow(saved_manifest), " figure files")
message("Region-pair summary: ", region_pair_summary_path)
message("Species palette: ", palette_path)
message("Manifest: ", manifest_path)
