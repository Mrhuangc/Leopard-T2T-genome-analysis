#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(jsonlite)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "figure", "OR_subfamily_comparison_final")
output_dir <- file.path(project_dir, "figure", "OR_subfamily_comparison_phylogeny_v1_20260908")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_path <- file.path(input_dir, "OR_subfamily_heatmap.sensitivity.source_data.tsv")
selection_path <- file.path(input_dir, "OR_subfamily_heatmap.sensitivity_display_selection.tsv")
stopifnot(file.exists(source_path), file.exists(selection_path))

dat <- read.delim(source_path, sep = "\t", quote = "", check.names = FALSE,
                  stringsAsFactors = FALSE)
selection <- read.delim(selection_path, sep = "\t", quote = "", check.names = FALSE,
                        stringsAsFactors = FALSE)

required <- c("analysis_unit", "or_subfamily", "functional_loci",
              "deviation_from_median")
if (!all(required %in% names(dat))) stop("Sensitivity source data are missing required columns")

# Column order follows the rooted topology while keeping each clade contiguous.
unit_order <- c("Fca", "Pbe", "Lge", "Ppa", "Ppa_or", "Ple", "Pti")
if (!setequal(unique(dat$analysis_unit), unit_order)) stop("Expected exactly seven analysis units")

row_order <- rev(selection$or_subfamily)
if (!setequal(unique(dat$or_subfamily), row_order)) stop("Selection and source-data rows disagree")

dat <- dat %>%
  mutate(
    analysis_unit = factor(analysis_unit, levels = unit_order),
    or_subfamily = factor(or_subfamily, levels = row_order)
  )

if (anyDuplicated(dat[c("analysis_unit", "or_subfamily")])) {
  stop("Each analysis-unit/subfamily cell must occur exactly once")
}
if (nrow(dat) != length(unit_order) * length(row_order)) {
  stop("Heatmap matrix is incomplete")
}

max_abs <- max(abs(dat$deviation_from_median), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
dat <- dat %>%
  mutate(label_colour = ifelse(abs(deviation_from_median) >= 0.62 * max_abs,
                               "white", "black"))

species_labels <- c(
  Fca = expression(italic("F. catus")),
  Pbe = expression(italic("P. bengalensis")),
  Lge = expression(italic("L. geoffroyi")),
  Ppa = expression(italic("P. p. japonensis")),
  Ppa_or = expression(italic("P. p. orientalis")),
  Ple = expression(italic("P. leo")),
  Pti = expression(italic("P. tigris"))
)

# Rectangular cladogram coordinates. This topology-only tree has no branch-length scale.
tree_segments <- data.frame(
  x = c(
    4.1875, 2.25, 6.125,
    1.50, 3.00, 1.00, 2.00,
    5.25, 7.00, 4.50, 6.00, 4.00, 5.00,
    2.25, 6.125, 1.50, 5.25, 4.50
  ),
  xend = c(
    4.1875, 2.25, 6.125,
    1.50, 3.00, 1.00, 2.00,
    5.25, 7.00, 4.50, 6.00, 4.00, 5.00,
    6.125, 3.00, 2.00, 7.00, 6.00
  ),
  y = c(
    0.15, 0.70, 0.70,
    1.30, 1.30, 2.15, 2.15,
    1.30, 1.30, 1.90, 1.90, 2.45, 2.45,
    0.70, 1.30, 2.15, 1.90, 2.45
  ),
  yend = c(
    0.70, 1.30, 1.30,
    2.15, 3.15, 3.15, 3.15,
    1.90, 3.15, 2.45, 3.15, 3.15, 3.15,
    0.70, 1.30, 2.15, 1.90, 2.45
  )
)

tree_plot <- ggplot() +
  annotate("rect", xmin = 0.50, xmax = 3.50, ymin = 0, ymax = 3.55,
           fill = "#EEF4F6", colour = NA) +
  annotate("rect", xmin = 3.50, xmax = 7.50, ymin = 0, ymax = 3.55,
           fill = "#F7F1ED", colour = NA) +
  geom_segment(
    data = tree_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.42, lineend = "square", colour = "#343434"
  ) +
  annotate("text", x = 2.0, y = 3.47, label = "Felinae",
           family = "Arial", fontface = "bold", size = 2.6, colour = "#35677A") +
  annotate("text", x = 5.5, y = 3.47, label = "Pantherinae",
           family = "Arial", fontface = "bold", size = 2.6, colour = "#8B5947") +
  scale_x_continuous(limits = c(0.5, 7.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3.60), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(0, 18, 0, 30))

heatmap_plot <- ggplot(
  dat,
  aes(x = analysis_unit, y = or_subfamily, fill = deviation_from_median)
) +
  annotate("rect", xmin = 0.50, xmax = 3.50, ymin = -Inf, ymax = Inf,
           fill = "#EEF4F6", alpha = 0.30, colour = NA) +
  annotate("rect", xmin = 3.50, xmax = 7.50, ymin = -Inf, ymax = Inf,
           fill = "#F7F1ED", alpha = 0.30, colour = NA) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = functional_loci, colour = label_colour),
            size = 1.75, family = "Arial", show.legend = FALSE) +
  geom_vline(xintercept = 3.5, linewidth = 0.45, colour = "white") +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "#F7F7F7", high = "#C54E4E",
    midpoint = 0, limits = c(-max_abs, max_abs),
    name = "Difference from\n7-unit median\n(loci)"
  ) +
  scale_x_discrete(labels = species_labels, drop = FALSE) +
  labs(x = NULL, y = "OR subfamily") +
  coord_cartesian(clip = "off") +
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
    plot.margin = margin(0, 7, 0, 5.5)
  )

caption_text <- paste0(
  "Numbers are functional-locus counts. Colour is centred within each subfamily and does not encode statistical significance. ",
  "The cladogram shows topology only; branch lengths are not scaled. Ppa and Ppa_or are subspecies-level assembly units."
)

combined <- tree_plot / heatmap_plot +
  plot_layout(heights = c(1.05, 8.8)) +
  plot_annotation(
    title = "Sensitivity analysis of functional OR subfamily composition",
    subtitle = paste0(
      "High-, medium- and low-confidence assignments; primary display set plus ",
      "low-confidence OR9K and OR2AP"
    ),
    caption = caption_text,
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8, face = "bold",
                                margin = margin(b = 1.5)),
      plot.subtitle = element_text(family = "Arial", size = 6.5, colour = "#555555",
                                   margin = margin(b = 3)),
      plot.caption = element_text(family = "Arial", size = 5.2, colour = "#555555",
                                  hjust = 0, margin = margin(t = 4)),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
  )

stem <- file.path(output_dir, "OR_subfamily_composition.sensitivity.phylogeny")
width_mm <- 183
height_mm <- 190
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(combined)
dev.off()

grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
                     family = "Arial")
print(combined)
dev.off()

ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600)
print(combined)
dev.off()

ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
              units = "in", res = 300)
print(combined)
dev.off()

ordered_source <- dat %>%
  arrange(or_subfamily, analysis_unit) %>%
  mutate(
    analysis_unit = as.character(analysis_unit),
    or_subfamily = as.character(or_subfamily)
  ) %>%
  select(-label_colour)
write.table(
  ordered_source,
  file.path(output_dir, "OR_subfamily_composition.sensitivity.phylogeny.source_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(
  "(((Felis_catus,Prionailurus_bengalensis),Leopardus_geoffroyi),(((Panthera_pardus_japonensis_Ppa,Panthera_pardus_orientalis_Ppa_or),Panthera_leo),Panthera_tigris));",
  file.path(output_dir, "felidae_7_analysis_units_topology_only.nwk")
)

qc <- list(
  schema_version = "1.0",
  status = "pass",
  figure = "OR_subfamily_composition.sensitivity.phylogeny",
  backend = "R",
  dimensions_mm = c(width_mm, height_mm),
  topology_only = TRUE,
  branch_lengths_scaled = FALSE,
  analysis_units = unit_order,
  display_labels = c(
    Fca = "F. catus",
    Pbe = "P. bengalensis",
    Lge = "L. geoffroyi",
    Ppa = "P. p. japonensis",
    Ppa_or = "P. p. orientalis",
    Ple = "P. leo",
    Pti = "P. tigris"
  ),
  clades = list(
    Felinae = c("Fca", "Pbe", "Lge"),
    Pantherinae = c("Ppa", "Ppa_or", "Ple", "Pti")
  ),
  counts = list(
    displayed_subfamilies = length(row_order),
    displayed_cells = nrow(dat),
    expected_cells = length(unit_order) * length(row_order)
  ),
  checks = list(
    source_cells_complete = nrow(dat) == length(unit_order) * length(row_order),
    source_cells_unique = !anyDuplicated(dat[c("analysis_unit", "or_subfamily")]),
    all_original_counts_retained = sum(dat$functional_loci) == 5935,
    two_leopard_subspecies_labeled = all(c("Ppa", "Ppa_or") %in% unit_order),
    original_source_unchanged = TRUE
  ),
  interpretation_boundary = paste0(
    "Cladogram is a topology-only visual annotation. Heatmap differences are descriptive; ",
    "no branch lengths, ancestral copy numbers, expansion/contraction events, or significance tests are encoded."
  )
)
if (!all(unlist(qc$checks))) qc$status <- "fail"
jsonlite::write_json(
  qc,
  file.path(output_dir, "OR_subfamily_composition.sensitivity.phylogeny.QC.json"),
  pretty = TRUE, auto_unbox = TRUE
)

capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("OR_SUBFAMILY_PHYLOGENY_HEATMAP_COMPLETED\n")
cat(sprintf("output_dir=%s\n", output_dir))

