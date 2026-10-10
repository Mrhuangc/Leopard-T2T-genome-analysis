#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_or_subfamily_absolute_copy_number_heatmap.R <results_dir> <figure_dir>")
}

results_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- normalizePath(figure_dir, winslash = "/", mustWork = TRUE)

source_path <- file.path(results_dir, "absolute_heatmap_top30.source_data.tsv")
metadata_path <- file.path(results_dir, "species_metadata.tsv")
stopifnot(file.exists(source_path), file.exists(metadata_path))

dat <- read.delim(source_path, sep = "\t", quote = "", check.names = FALSE,
                  stringsAsFactors = FALSE, na.strings = character())
metadata <- read.delim(metadata_path, sep = "\t", quote = "", check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = character())
required <- c("display_rank", "or_subfamily", "broad_clade", "analysis_unit", "absolute_copy_number")
if (!all(required %in% names(dat))) stop("Heatmap source data lacks required columns")

unit_order <- c("Fca", "Lge", "Pbe", "Ple", "Ppa", "Ppa_or", "Pti")
if (!setequal(unique(dat$analysis_unit), unit_order)) stop("Expected exactly seven analysis units")
if (nrow(dat) != 30L * length(unit_order)) stop("Expected 30 x 7 heatmap cells")
if (anyDuplicated(dat[c("or_subfamily", "analysis_unit")])) stop("Duplicate heatmap cells")

dat$absolute_copy_number <- as.integer(dat$absolute_copy_number)
if (anyNA(dat$absolute_copy_number) || any(dat$absolute_copy_number < 0L)) {
  stop("Absolute copy numbers must be non-negative integers")
}
row_order <- unique(dat$or_subfamily[order(dat$display_rank)])
dat$or_subfamily <- factor(dat$or_subfamily, levels = rev(row_order))
dat$analysis_unit <- factor(dat$analysis_unit, levels = unit_order)
max_count <- max(dat$absolute_copy_number)
dat$label_colour <- ifelse(dat$absolute_copy_number >= 0.58 * max_count, "white", "black")

species_labels <- c(
  Fca = expression(italic("F. catus")),
  Lge = expression(italic("L. geoffroyi")),
  Pbe = expression(italic("P. bengalensis")),
  Ple = expression(italic("P. leo")),
  Ppa = expression(italic("P. p. japonensis")),
  Ppa_or = expression(italic("P. p. orientalis")),
  Pti = expression(italic("P. tigris"))
)

plot <- ggplot(dat, aes(x = analysis_unit, y = or_subfamily, fill = absolute_copy_number)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = absolute_copy_number, colour = label_colour),
            family = "Arial", size = 2.0, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits = c(0, max_count),
    name = "Functional OR\nloci"
  ) +
  scale_x_discrete(labels = species_labels, drop = FALSE) +
  labs(
    x = NULL,
    y = "Standard OR subfamily",
    title = "Absolute functional OR copy number across seven felid assemblies",
    subtitle = "High- and medium-confidence assignments; 30 subfamilies with the largest across-unit count ranges",
    caption = paste0(
      "Numbers and colour both encode absolute functional-locus counts. Differences are descriptive; ",
      "they are not, by themselves, evidence of expansion or contraction."
    )
  ) +
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
    plot.title = element_text(size = 8.5, face = "bold"),
    plot.subtitle = element_text(size = 6.5, colour = "#444444"),
    plot.caption = element_text(size = 5.5, colour = "#555555", hjust = 0),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

width_in <- 183 / 25.4
height_in <- 176 / 25.4
stem <- file.path(figure_dir, "OR_subfamily_absolute_copy_number.high_medium.top30")

svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(plot)
dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(plot)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, compression = "lzw")
print(plot)
dev.off()
ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
              units = "in", res = 300)
print(plot)
dev.off()

capture.output(sessionInfo(), file = file.path(figure_dir, "R_sessionInfo.txt"))
outputs <- paste0(stem, c(".svg", ".pdf", ".tiff", ".png"))
if (!all(file.exists(outputs)) || any(file.info(outputs)$size <= 0)) stop("One or more figure exports failed")

cat("ABSOLUTE_COPY_NUMBER_HEATMAP_COMPLETE\n")
cat(sprintf("displayed_subfamilies=%d\n", length(row_order)))
cat(sprintf("heatmap_cells=%d\n", nrow(dat)))
cat(sprintf("maximum_displayed_copy_number=%d\n", max_count))
cat(sprintf("figure_dir=%s\n", figure_dir))
