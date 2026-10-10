#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript build_fig5b_subfamily_count_data.R <input_snapshot_dir> <output_dir>")
}

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

counts_path <- file.path(input_dir, "species_OR_subfamily_counts.tsv")
metadata_path <- file.path(input_dir, "species_metadata.tsv")
completion_path <- file.path(input_dir, "ASSIGNMENT_COMPLETE.json")
stopifnot(file.exists(counts_path), file.exists(metadata_path), file.exists(completion_path))

counts <- read.delim(counts_path, sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                     stringsAsFactors = FALSE, na.strings = c("NA", ""))
metadata <- read.delim(metadata_path, sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = c("NA", ""))
completion <- jsonlite::read_json(completion_path, simplifyVector = TRUE)

required_counts <- c("analysis_unit", "or_family", "or_subfamily",
                     "assignment_confidence", "functional_loci")
required_meta <- c("analysis_unit", "scientific_name")
if (!all(required_counts %in% names(counts))) stop("Counts table is missing required columns")
if (!all(required_meta %in% names(metadata))) stop("Metadata table is missing required columns")

counts$functional_loci <- as.integer(counts$functional_loci)
if (anyNA(counts$functional_loci) || any(counts$functional_loci < 0L)) {
  stop("functional_loci must contain non-negative integers")
}

unit_order <- c("Fca", "Lge", "Pbe", "Ple", "Ppa", "Ppa_or", "Pti")
if (!setequal(unique(counts$analysis_unit), unit_order)) stop("Expected exactly seven analysis units")
if (!setequal(unique(metadata$analysis_unit), unit_order)) stop("Metadata does not match seven analysis units")

allowed_confidence <- c("high", "medium", "low", "unresolved")
unexpected_confidence <- setdiff(unique(counts$assignment_confidence), allowed_confidence)
if (length(unexpected_confidence) > 0L) stop("Unexpected confidence category")

valid_named <- counts %>%
  filter(!is.na(or_subfamily), nzchar(or_subfamily),
         !tolower(or_subfamily) %in% c("unresolved", "unknown", "na"))

aggregate_matrix <- function(confidence_levels) {
  valid_named %>%
    filter(assignment_confidence %in% confidence_levels) %>%
    group_by(analysis_unit, or_family, or_subfamily) %>%
    summarise(functional_loci = sum(functional_loci), .groups = "drop") %>%
    complete(
      analysis_unit = unit_order,
      nesting(or_family, or_subfamily),
      fill = list(functional_loci = 0L)
    ) %>%
    mutate(analysis_unit = factor(analysis_unit, levels = unit_order))
}

main_long <- aggregate_matrix(c("high", "medium"))
sensitivity_long <- aggregate_matrix(c("high", "medium", "low"))

main_metrics <- main_long %>%
  group_by(or_family, or_subfamily) %>%
  summarise(
    minimum = min(functional_loci),
    maximum = max(functional_loci),
    count_range = maximum - minimum,
    total = sum(functional_loci),
    median_count = median(functional_loci),
    .groups = "drop"
  ) %>%
  arrange(desc(count_range), desc(maximum), desc(total), or_subfamily)

focal_subfamilies <- c("OR6C", "OR9K", "OR2AP")
selected_metrics <- slice_head(main_metrics, n = 30L) %>%
  arrange(desc(count_range), desc(maximum), desc(total), or_subfamily)

sensitivity_metrics <- sensitivity_long %>%
  group_by(or_family, or_subfamily) %>%
  summarise(
    minimum = min(functional_loci),
    maximum = max(functional_loci),
    count_range = maximum - minimum,
    total = sum(functional_loci),
    median_count = median(functional_loci),
    .groups = "drop"
  )

sensitivity_selected_metrics <- bind_rows(
  selected_metrics,
  filter(sensitivity_metrics, or_subfamily %in% focal_subfamilies)
) %>%
  distinct(or_family, or_subfamily, .keep_all = TRUE) %>%
  arrange(desc(count_range), desc(maximum), desc(total), or_subfamily)

selected_order <- rev(selected_metrics$or_subfamily)
sensitivity_selected_order <- rev(sensitivity_selected_metrics$or_subfamily)

prepare_display <- function(long_table, selected_table, row_order) {
  medians <- long_table %>%
    group_by(or_subfamily) %>%
    summarise(row_median = median(functional_loci), .groups = "drop")
  long_table %>%
    filter(or_subfamily %in% selected_table$or_subfamily) %>%
    left_join(medians, by = "or_subfamily") %>%
    mutate(
      deviation_from_median = functional_loci - row_median,
      or_subfamily = factor(or_subfamily, levels = row_order),
      analysis_unit = factor(analysis_unit, levels = unit_order)
    )
}

main_display <- prepare_display(main_long, selected_metrics, selected_order)
sensitivity_display <- prepare_display(
  sensitivity_long, sensitivity_selected_metrics, sensitivity_selected_order
)

max_abs <- max(abs(c(main_display$deviation_from_median,
                     sensitivity_display$deviation_from_median)), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

species_labels <- c(
  Fca = expression(italic("F. catus")),
  Lge = expression(italic("L. geoffroyi")),
  Pbe = expression(italic("P. bengalensis")),
  Ple = expression(italic("P. leo")),
  Ppa = expression(italic("P. pardus")),
  Ppa_or = expression(italic("P. p. orientalis")),
  Pti = expression(italic("P. tigris"))
)

theme_or <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7),
      axis.text.x = element_text(size = 6.5, angle = 40, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.8, colour = "black"),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6.5, colour = "#555555"),
      plot.caption = element_text(size = 5.5, colour = "#555555", hjust = 0),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
}

make_heatmap <- function(display_data, title, subtitle) {
  display_data <- display_data %>%
    mutate(label_colour = ifelse(abs(deviation_from_median) >= 0.62 * max_abs,
                                 "white", "black"))
  ggplot(display_data, aes(x = analysis_unit, y = or_subfamily,
                           fill = deviation_from_median)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = functional_loci, colour = label_colour),
              size = 1.75, family = "Arial", show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = "#3B6FB6", mid = "#F7F7F7", high = "#C54E4E",
      midpoint = 0, limits = c(-max_abs, max_abs),
      name = "Difference from\n7-unit median\n(loci)"
    ) +
    scale_x_discrete(labels = species_labels, drop = FALSE) +
    labs(
      x = NULL, y = "OR subfamily", title = title, subtitle = subtitle,
      caption = paste0(
        "Numbers are functional-locus counts. Colour is centred within each subfamily; ",
        "it does not encode statistical significance."
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_or()
}

main_plot <- make_heatmap(
  main_display,
  "Functional OR subfamily composition across seven felid assemblies",
  "High- and medium-confidence assignments; 30 subfamilies with the largest assembly-relative count ranges"
)

sensitivity_plot <- make_heatmap(
  sensitivity_display,
  "Sensitivity analysis of functional OR subfamily composition",
  "High-, medium- and low-confidence assignments; primary display set plus low-confidence OR9K and OR2AP"
)

save_plot <- function(plot, stem, width_mm = 183, height_mm = 170, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
                       family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = dpi)
  print(plot)
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
                units = "in", res = 300)
  print(plot)
  dev.off()
}

main_stem <- file.path(output_dir, "OR_subfamily_composition.primary")
sensitivity_stem <- file.path(output_dir, "OR_subfamily_composition.sensitivity")
save_plot(main_plot, main_stem)
save_plot(sensitivity_plot, sensitivity_stem)

write.table(main_long, file.path(output_dir, "OR_subfamily_counts.primary_complete.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sensitivity_long, file.path(output_dir, "OR_subfamily_counts.sensitivity_complete.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(main_display, file.path(output_dir, "OR_subfamily_heatmap.primary.source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sensitivity_display,
            file.path(output_dir, "OR_subfamily_heatmap.sensitivity.source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(selected_metrics,
            file.path(output_dir, "OR_subfamily_heatmap.primary_display_selection.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sensitivity_selected_metrics,
            file.path(output_dir, "OR_subfamily_heatmap.sensitivity_display_selection.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

all_subfamilies <- sort(unique(c(
  as.character(main_long$or_subfamily),
  as.character(sensitivity_long$or_subfamily)
)))

comparison_long <- expand_grid(
  analysis_unit = unit_order,
  or_subfamily = all_subfamilies
) %>%
  left_join(
    main_long %>%
      transmute(
        analysis_unit = as.character(analysis_unit),
        or_subfamily = as.character(or_subfamily),
        primary_count = functional_loci
      ),
    by = c("analysis_unit", "or_subfamily")
  ) %>%
  left_join(
    sensitivity_long %>%
      transmute(
        analysis_unit = as.character(analysis_unit),
        or_subfamily = as.character(or_subfamily),
        sensitivity_count = functional_loci
      ),
    by = c("analysis_unit", "or_subfamily")
  ) %>%
  mutate(
    primary_count = replace_na(primary_count, 0L),
    sensitivity_count = replace_na(sensitivity_count, 0L)
  )

sensitivity_summary <- comparison_long %>%
  group_by(analysis_unit) %>%
  summarise(
    primary_loci = sum(primary_count),
    sensitivity_loci = sum(sensitivity_count),
    low_confidence_loci_added = sensitivity_loci - primary_loci,
    spearman_rho_all_176 = suppressWarnings(
      cor(primary_count, sensitivity_count, method = "spearman")
    ),
    spearman_rho_primary_82 = suppressWarnings(
      cor(
        primary_count[or_subfamily %in% unique(as.character(main_long$or_subfamily))],
        sensitivity_count[or_subfamily %in% unique(as.character(main_long$or_subfamily))],
        method = "spearman"
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(analysis_unit = factor(analysis_unit, levels = unit_order)) %>%
  arrange(analysis_unit) %>%
  mutate(analysis_unit = as.character(analysis_unit))
write.table(sensitivity_summary,
            file.path(output_dir, "OR_subfamily_sensitivity.summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

observed_named_total <- sum(valid_named$functional_loci)
primary_total <- sum(main_long$functional_loci)
sensitivity_total <- sum(sensitivity_long$functional_loci)
functional_total <- as.integer(completion$counts$functional_felid_loci)
unresolved_or_unnamed_total <- functional_total - observed_named_total

qc <- list(
  schema_version = "1.0",
  status = "pass",
  observational_unit = "genome-assembly analysis unit",
  statistical_scope = "descriptive; loci are not treated as independent replicates",
  inclusion_rules = list(
    primary = c("high", "medium"),
    sensitivity = c("high", "medium", "low"),
    unresolved = "retained as unresolved and excluded from named-subfamily heatmaps",
    primary_display = "top 30 subfamilies by across-unit count range",
    sensitivity_display = "primary display set plus low-confidence focal OR9K and OR2AP"
  ),
  counts = list(
    analysis_units = length(unit_order),
    functional_loci_upstream = functional_total,
    primary_named_loci = primary_total,
    sensitivity_named_loci = sensitivity_total,
    named_loci_all_confidences = observed_named_total,
    unresolved_or_unnamed_loci = unresolved_or_unnamed_total,
    primary_subfamilies = n_distinct(main_long$or_subfamily),
    sensitivity_subfamilies = n_distinct(sensitivity_long$or_subfamily),
    primary_displayed_subfamilies = nrow(selected_metrics),
    sensitivity_displayed_subfamilies = nrow(sensitivity_selected_metrics)
  ),
  checks = list(
    upstream_status_pass = identical(completion$status, "pass"),
    seven_analysis_units_present = setequal(unique(as.character(main_long$analysis_unit)), unit_order),
    nonnegative_primary_counts = all(main_long$functional_loci >= 0L),
    sensitivity_contains_primary = sensitivity_total >= primary_total,
    OR6C_displayed_in_primary = "OR6C" %in% selected_metrics$or_subfamily,
    all_focal_subfamilies_displayed_in_sensitivity = all(
      focal_subfamilies %in% sensitivity_selected_metrics$or_subfamily
    ),
    primary_source_rows_match_display = nrow(main_display) == nrow(selected_metrics) * length(unit_order),
    sensitivity_source_rows_match_display = nrow(sensitivity_display) ==
      nrow(sensitivity_selected_metrics) * length(unit_order),
    sensitivity_summary_primary_total_matches = sum(sensitivity_summary$primary_loci) == primary_total,
    sensitivity_summary_full_total_matches = sum(sensitivity_summary$sensitivity_loci) == sensitivity_total,
    sensitivity_summary_added_total_matches =
      sum(sensitivity_summary$low_confidence_loci_added) == sensitivity_total - primary_total
  )
)
if (!all(unlist(qc$checks))) qc$status <- "fail"
jsonlite::write_json(qc, file.path(output_dir, "OR_subfamily_comparison.QC.json"),
                     pretty = TRUE, auto_unbox = TRUE)

capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("OR_SUBFAMILY_COMPARISON_COMPLETED\n")
cat(sprintf("primary_loci=%d\n", primary_total))
cat(sprintf("sensitivity_loci=%d\n", sensitivity_total))
cat(sprintf("primary_displayed_subfamilies=%d\n", nrow(selected_metrics)))
cat(sprintf("sensitivity_displayed_subfamilies=%d\n", nrow(sensitivity_selected_metrics)))
cat(sprintf("output_dir=%s\n", output_dir))
