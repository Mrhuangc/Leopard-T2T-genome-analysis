#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(readr)
  library(scales)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_fig5d_chrd1_or_arrays.R <input_dir> <output_dir>")
}

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

locus_file <- file.path(
  input_dir,
  "OR_ChrD1_macro_OR_locus_tracks_broken_axis.source_data.tsv"
)
qc_file <- file.path(
  input_dir,
  "OR_ChrD1_macro_OR_locus_tracks_broken_axis.QC.json"
)
if (!file.exists(locus_file) || !file.exists(qc_file)) {
  stop("Required audited source-data or QC file is missing from: ", input_dir)
}

loci <- read_tsv(locus_file, show_col_types = FALSE)
source_qc <- read_json(qc_file, simplifyVector = TRUE)

required_columns <- c(
  "catalog_id", "analysis_unit", "classification", "array_identity",
  "oriented_mid_mb", "plot_x", "array_start_mb", "array_end_mb"
)
if (!all(required_columns %in% names(loci))) {
  stop("Source data lack: ", paste(setdiff(required_columns, names(loci)), collapse = ", "))
}

unit_order <- c("Fca", "Lge", "Pbe", "Ple", "Ppa", "Ppa_or", "Pti")
species_labels <- c(
  Fca = "F. catus",
  Lge = "L. geoffroyi",
  Pbe = "P. bengalensis",
  Ple = "P. leo",
  Ppa = "P. p. japonensis",
  Ppa_or = "P. p. orientalis",
  Pti = "P. tigris"
)
array_order <- c("gamma-only array", "alpha-dominated mixed array")

plot_data <- loci |>
  mutate(
    analysis_unit = factor(analysis_unit, levels = rev(unit_order)),
    classification = factor(
      tolower(classification),
      levels = c("functional", "pseudogene", "truncated")
    ),
    array_identity = factor(array_identity, levels = array_order)
  )

if (nrow(plot_data) != 2777L || n_distinct(plot_data$catalog_id) != 2777L) {
  stop("Expected 2,777 unique loci; observed ", nrow(plot_data), " rows and ",
       n_distinct(plot_data$catalog_id), " unique IDs")
}
if (anyNA(plot_data$analysis_unit) || anyNA(plot_data$classification) ||
    anyNA(plot_data$array_identity)) {
  stop("At least one record has an unexpected species, class, or array label")
}

gap_start <- as.numeric(source_qc$broken_axis$omitted_from_mb)
gap_end <- as.numeric(source_qc$broken_axis$omitted_to_mb)
compressed_width <- as.numeric(source_qc$broken_axis$displayed_break_width_mb)
omitted_width <- gap_end - gap_start
if (!is.finite(omitted_width) || omitted_width <= 0 || compressed_width <= 0) {
  stop("Invalid broken-axis parameters in source QC")
}

fold_coordinate <- function(x) {
  ifelse(
    x <= gap_start,
    x,
    ifelse(
      x >= gap_end,
      gap_start + compressed_width + (x - gap_end),
      gap_start + compressed_width * (x - gap_start) / omitted_width
    )
  )
}

x_min_original <- min(plot_data$array_start_mb)
x_max_original <- max(plot_data$array_end_mb)
x_min_plot <- fold_coordinate(x_min_original)
x_max_plot <- fold_coordinate(x_max_original)
break_x <- gap_start + compressed_width / 2

track_data <- plot_data |>
  count(analysis_unit, name = "total_loci") |>
  mutate(
    x = x_min_plot,
    xend = x_max_plot,
    count_x = x_max_plot + 0.26,
    count_label = paste0("n = ", comma(total_loci))
  )

array_ranges <- plot_data |>
  group_by(array_identity) |>
  summarise(
    xmin_original = min(array_start_mb),
    xmax_original = max(array_end_mb),
    .groups = "drop"
  ) |>
  mutate(
    xmin = fold_coordinate(xmin_original),
    xmax = fold_coordinate(xmax_original),
    fill = c("#E9F5F2", "#EEF2F8"),
    label = c("Gamma-only OR array", "Alpha-dominated mixed OR array"),
    label_colour = c("#27665F", "#465D80")
  )

left_break_values <- pretty(c(x_min_original, gap_start), n = 4)
left_break_values <- left_break_values[
  left_break_values >= x_min_original & left_break_values <= gap_start
]
right_break_values <- pretty(c(gap_end, x_max_original), n = 4)
right_break_values <- right_break_values[
  right_break_values >= gap_end & right_break_values <= x_max_original
]
axis_original_values <- unique(c(left_break_values, right_break_values))
axis_plot_values <- fold_coordinate(axis_original_values)

class_palette <- c(
  functional = "#007F79",
  pseudogene = "#E99518",
  truncated = "#B23A73"
)

p <- ggplot() +
  geom_rect(
    data = array_ranges,
    aes(xmin = xmin, xmax = xmax, ymin = 0.48, ymax = 7.48, fill = fill),
    colour = NA,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  geom_segment(
    data = track_data,
    aes(x = x, xend = xend, y = analysis_unit, yend = analysis_unit),
    linewidth = 0.31,
    colour = "#B8B8B8"
  ) +
  geom_segment(
    data = plot_data,
    aes(
      x = plot_x,
      xend = plot_x,
      y = as.numeric(analysis_unit) - 0.30,
      yend = as.numeric(analysis_unit) + 0.30,
      colour = classification
    ),
    linewidth = 0.55,
    alpha = 0.96,
    lineend = "butt"
  ) +
  geom_label(
    data = track_data,
    aes(x = break_x, y = analysis_unit, label = "//"),
    family = "Arial",
    fontface = "bold",
    size = 2.75,
    linewidth = 0,
    label.padding = grid::unit(0.015, "lines"),
    fill = "white",
    colour = "#4B4B4B"
  ) +
  geom_text(
    data = track_data,
    aes(x = count_x, y = analysis_unit, label = count_label),
    family = "Arial",
    hjust = 0,
    size = 2.55,
    colour = "#333333"
  ) +
  annotate(
    "text",
    x = mean(c(array_ranges$xmin[[1]], array_ranges$xmax[[1]])),
    y = 7.67,
    label = array_ranges$label[[1]],
    family = "Arial",
    fontface = "bold",
    size = 2.95,
    colour = array_ranges$label_colour[[1]]
  ) +
  annotate(
    "text",
    x = mean(c(array_ranges$xmin[[2]], array_ranges$xmax[[2]])),
    y = 7.67,
    label = array_ranges$label[[2]],
    family = "Arial",
    fontface = "bold",
    size = 2.95,
    colour = array_ranges$label_colour[[2]]
  ) +
  scale_colour_manual(
    values = class_palette,
    breaks = c("functional", "pseudogene", "truncated"),
    labels = c("Functional", "Pseudogene", "Truncated"),
    name = "OR locus class"
  ) +
  scale_y_discrete(
    limits = rev(unit_order),
    labels = species_labels[rev(unit_order)],
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = axis_plot_values,
    labels = number(axis_original_values, accuracy = 1),
    expand = expansion(mult = c(0.01, 0))
  ) +
  coord_cartesian(
    xlim = c(x_min_plot, x_max_plot + 0.73),
    ylim = c(0.45, 7.84),
    clip = "off"
  ) +
  labs(
    x = "Orientation-normalized position along ChrD1 (Mb)",
    y = NULL
  ) +
  theme_classic(base_size = 8.4, base_family = "Arial") +
  theme(
    axis.text.y = element_text(face = "italic", size = 8.1, colour = "#222222"),
    axis.text.x = element_text(size = 7.6, colour = "#222222"),
    axis.title.x = element_text(size = 8.2, margin = margin(t = 5)),
    axis.line = element_line(linewidth = 0.35, colour = "#1F1F1F"),
    axis.ticks = element_line(linewidth = 0.32, colour = "#1F1F1F"),
    axis.ticks.length = grid::unit(1.5, "mm"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8.0),
    legend.text = element_text(size = 7.8),
    legend.key.width = grid::unit(6.5, "mm"),
    legend.key.height = grid::unit(2.6, "mm"),
    legend.spacing.x = grid::unit(1.3, "mm"),
    legend.margin = margin(0, 0, 2, 0),
    panel.grid = element_blank(),
    plot.margin = margin(2, 27, 3, 3)
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.0)))

prefix <- file.path(output_dir, "OR_ChrD1_panelD_wide")
svg_file <- paste0(prefix, ".svg")
pdf_file <- paste0(prefix, ".pdf")
png_file <- paste0(prefix, ".png")
tiff_file <- paste0(prefix, ".tiff")

svglite::svglite(svg_file, width = 183 / 25.4, height = 72 / 25.4, bg = "white")
print(p)
grDevices::dev.off()

grDevices::cairo_pdf(
  pdf_file,
  width = 183 / 25.4,
  height = 72 / 25.4,
  family = "Arial",
  bg = "white"
)
print(p)
grDevices::dev.off()

ggsave(
  filename = png_file, plot = p,
  width = 183, height = 72, units = "mm", dpi = 600, bg = "white"
)
ggsave(
  filename = tiff_file, plot = p,
  width = 183, height = 72, units = "mm", dpi = 600,
  compression = "lzw", bg = "white"
)

output_qc <- list(
  schema_version = "1.0",
  status = "pass",
  figure = "OR_ChrD1_panelD_wide",
  source_figure = "OR_ChrD1_macro_OR_locus_tracks_broken_axis",
  data_policy = "all audited loci shown; no sampling, binning, or smoothing",
  dimensions_mm = list(width = 183, height = 72),
  counts = list(
    species = n_distinct(plot_data$analysis_unit),
    arrays = n_distinct(plot_data$array_identity),
    loci = nrow(plot_data),
    functional = sum(plot_data$classification == "functional"),
    pseudogene = sum(plot_data$classification == "pseudogene"),
    truncated = sum(plot_data$classification == "truncated")
  ),
  broken_axis = source_qc$broken_axis,
  checks = list(
    all_loci_unique = n_distinct(plot_data$catalog_id) == nrow(plot_data),
    no_locus_hidden_by_break = !any(
      plot_data$oriented_mid_mb > gap_start & plot_data$oriented_mid_mb < gap_end
    ),
    seven_species_present = n_distinct(plot_data$analysis_unit) == 7L,
    two_arrays_present = n_distinct(plot_data$array_identity) == 2L,
    northern_chinese_leopard_label_updated = species_labels[["Ppa"]] == "P. p. japonensis"
  ),
  palette = as.list(class_palette)
)
write_json(
  output_qc,
  paste0(prefix, ".QC.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_tsv(
  track_data |>
    mutate(
      analysis_unit = as.character(analysis_unit),
      species_label = unname(species_labels[analysis_unit])
    ),
  paste0(prefix, ".counts.tsv")
)

cat("FIGURE_COMPLETED\n")
cat("loci=", nrow(plot_data), "\n", sep = "")
cat("dimensions_mm=183x72\n")
cat("output_prefix=", prefix, "\n", sep = "")
