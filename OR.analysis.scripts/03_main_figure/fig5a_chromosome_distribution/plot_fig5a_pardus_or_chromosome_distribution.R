#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(ragg)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_fig5a_pardus_or_chromosome_distribution.R <project_root> <output_dir>")
}

project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE,
             na.strings = c("", "NA"))
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

source_dir <- file.path(project_root, "figure", "Ppa_OR_cluster_landscape_main_v1")
chromosome_path <- file.path(source_dir, "Ppa_chromosome_lengths.source_data.tsv")
locus_path <- file.path(source_dir, "Ppa_OR_loci_genomewide.source_data.tsv")
cluster_path <- file.path(source_dir, "Ppa_OR_clusters.source_data.tsv")

required <- c(chromosome_path, locus_path, cluster_path)
bad <- required[!file.exists(required) | file.info(required)$size == 0]
if (length(bad)) stop("Missing or empty required input: ", paste(basename(bad), collapse = ", "))

chromosomes <- read_tsv(chromosome_path) %>%
  transmute(seqid, chromosome_bp = as.numeric(chromosome_bp),
            length_mb = as.numeric(length_mb))
loci <- read_tsv(locus_path) %>%
  mutate(start_1 = as.numeric(start_1), end_1 = as.numeric(end_1),
         midpoint_mb = as.numeric(midpoint_mb),
         classification = tolower(classification))
clusters <- read_tsv(cluster_path) %>%
  mutate(start_1 = as.numeric(start_1), end_1 = as.numeric(end_1),
         start_mb = start_1 / 1e6, end_mb = end_1 / 1e6,
         midpoint_mb = (start_1 + end_1) / 2e6,
         member_count = as.integer(member_count))

chromosome_order <- c(
  "ChrA1", "ChrA2", "ChrA3",
  "ChrB1", "ChrB2", "ChrB3", "ChrB4",
  "ChrC1", "ChrC2",
  "ChrD1", "ChrD2", "ChrD3", "ChrD4",
  "ChrE1", "ChrE2", "ChrE3",
  "ChrF2", "ChrF3", "ChrX", "ChrY"
)

if (nrow(chromosomes) != 20L || !setequal(chromosomes$seqid, chromosome_order)) {
  stop("Expected exactly the 20 declared Ppa chromosomes")
}
if (nrow(loci) != 1042L || dplyr::n_distinct(loci$global_id) != 1042L) {
  stop("Expected 1,042 unique Ppa OR loci")
}
allowed_classes <- c("functional", "pseudogene", "truncated")
if (!setequal(unique(loci$classification), allowed_classes)) {
  stop("Unexpected or missing OR classification")
}
if (any(!loci$seqid %in% chromosome_order)) stop("OR locus found outside declared chromosomes")
if (any(loci$midpoint_mb < 0, na.rm = TRUE)) stop("Negative locus coordinate")

chromosomes <- chromosomes %>%
  mutate(order_index = match(seqid, chromosome_order),
         y = length(chromosome_order) - order_index + 1L) %>%
  arrange(order_index)

loci <- loci %>%
  mutate(classification = factor(classification, levels = allowed_classes)) %>%
  left_join(chromosomes %>% select(seqid, chromosome_bp, length_mb, y), by = "seqid")

if (any(is.na(loci$y))) stop("Chromosome join failed for one or more OR loci")
if (any(loci$end_1 > loci$chromosome_bp, na.rm = TRUE)) stop("OR locus exceeds chromosome length")

clusters <- clusters %>%
  left_join(chromosomes %>% select(sseqid = seqid, y), by = "sseqid") %>%
  mutate(
    cluster_class = case_when(
      cluster_id == "Ppa_ORCL000020" ~ "ChrD1 gamma array",
      cluster_id == "Ppa_ORCL000021" ~ "ChrD1 alpha-rich array",
      cluster_id == "Ppa_ORCL000019" ~ "ChrB4 OR6C-rich array",
      member_count >= 50L ~ "Other large OR cluster",
      TRUE ~ "Not highlighted"
    )
  )

large_clusters <- clusters %>% filter(member_count >= 50L)
focus_clusters <- clusters %>%
  filter(cluster_id %in% c("Ppa_ORCL000019", "Ppa_ORCL000020", "Ppa_ORCL000021")) %>%
  mutate(
    annotation = case_when(
      cluster_id == "Ppa_ORCL000020" ~ "185",
      cluster_id == "Ppa_ORCL000021" ~ "214",
      cluster_id == "Ppa_ORCL000019" ~ "OR6C-rich, 65"
    ),
    annotation_x = case_when(
      cluster_id == "Ppa_ORCL000019" ~ end_mb + 2.5,
      TRUE ~ midpoint_mb
    ),
    annotation_y = if_else(cluster_id == "Ppa_ORCL000019", y + 0.18, y + 0.56),
    annotation_hjust = if_else(cluster_id == "Ppa_ORCL000019", 0, 0.5)
  )
if (nrow(focus_clusters) != 3L) stop("Expected three focus clusters")

class_palette <- c(
  functional = "#007F79",
  pseudogene = "#E99518",
  truncated = "#B23A73"
)
cluster_palette <- c(
  "ChrD1 gamma array" = "#AEBFBC",
  "ChrD1 alpha-rich array" = "#AEBFBC",
  "ChrB4 OR6C-rich array" = "#AEBFBC",
  "Other large OR cluster" = "#AEBFBC"
)

plot_data_loci <- loci %>%
  mutate(
    track_ymin = case_when(
      classification == "functional" ~ y + 0.04,
      classification == "pseudogene" ~ y - 0.42,
      classification == "truncated" ~ y - 0.48
    ),
    track_ymax = case_when(
      classification == "functional" ~ y + 0.42,
      classification == "pseudogene" ~ y - 0.04,
      classification == "truncated" ~ y + 0.48
    )
  )

legend_key <- data.frame(
  classification = factor(allowed_classes, levels = allowed_classes),
  label = c("Functional", "Pseudogene", "Truncated"),
  x = 151,
  y = c(5.25, 4.45, 3.65)
)

p <- ggplot() +
  geom_segment(
    data = chromosomes,
    aes(x = 0, xend = length_mb, y = y, yend = y),
    linewidth = 3.15, colour = "#D8D8D5", lineend = "round"
  ) +
  geom_rect(
    data = large_clusters,
    aes(xmin = start_mb, xmax = end_mb, ymin = y - 0.34, ymax = y + 0.34,
        fill = cluster_class),
    alpha = 0.26, colour = NA
  ) +
  geom_segment(
    data = plot_data_loci,
    aes(x = midpoint_mb, xend = midpoint_mb,
        y = track_ymin, yend = track_ymax,
        colour = classification),
    linewidth = 0.48, alpha = 0.98, lineend = "butt"
  ) +
  geom_text(
    data = focus_clusters,
    aes(x = annotation_x, y = annotation_y, label = annotation, hjust = annotation_hjust),
    size = 2.35, fontface = "bold", colour = "#292929"
  ) +
  geom_segment(
    data = legend_key,
    aes(x = x, xend = x, y = y - 0.22, yend = y + 0.22,
        colour = classification),
    linewidth = 1.05, lineend = "butt", inherit.aes = FALSE
  ) +
  geom_text(
    data = legend_key,
    aes(x = x + 3.2, y = y, label = label),
    hjust = 0, vjust = 0.5, size = 2.35, colour = "#292929",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = class_palette, limits = allowed_classes, drop = FALSE,
    labels = c("Functional", "Pseudogene", "Truncated"), name = NULL,
    guide = "none"
  ) +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(chromosome_order), labels = rev(chromosome_order),
    expand = expansion(add = c(0.65, 0.80))
  ) +
  scale_x_continuous(
    name = "Chromosome position (Mb)",
    breaks = seq(0, 250, 50), limits = c(0, 250),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_classic(base_size = 7.0, base_family = "Arial") +
  theme(
    axis.title.y = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.5, colour = "#292929", margin = margin(r = 2.5)),
    axis.text.x = element_text(size = 6.2, colour = "#292929"),
    axis.title.x = element_text(size = 7.0, colour = "#292929", margin = margin(t = 4)),
    axis.line.x = element_line(linewidth = 0.38, colour = "#292929"),
    axis.ticks.x = element_line(linewidth = 0.34, colour = "#292929"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(4, 6, 2, 4)
  )

prefix <- file.path(output_dir, "Ppa_OR_chromosome_distribution_separated_tracks")
width_in <- 183 / 25.4
height_in <- 112 / 25.4

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in)
print(p)
dev.off()

grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in,
                     family = "Arial")
print(p)
dev.off()

ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, compression = "lzw", background = "white")
print(p)
dev.off()

ragg::agg_png(paste0(prefix, ".png"), width = width_in, height = height_in,
              units = "in", res = 300, background = "white")
print(p)
dev.off()

write_tsv(chromosomes %>% select(seqid, chromosome_bp, length_mb, order_index),
          file.path(output_dir, "Ppa_chromosome_lengths.source_data.tsv"))
write_tsv(loci %>% select(global_id, catalog_id, seqid, start_1, end_1, midpoint_mb,
                          classification, broad_clade, or_subfamily,
                          operational_cluster_id, operational_cluster_size),
          file.path(output_dir, "Ppa_OR_loci_genomewide.source_data.tsv"))
write_tsv(clusters %>% select(cluster_id, sseqid, start_1, end_1, span_bp,
                              member_count, functional_count, pseudogene_count,
                              truncated_count, clade_composition, cluster_class),
          file.path(output_dir, "Ppa_OR_clusters.source_data.tsv"))

class_counts <- table(loci$classification)
qc <- list(
  status = "pass",
  figure = "Ppa_OR_chromosome_distribution_separated_tracks",
  backend = "R",
  dimensions_mm = c(width = 183, height = 112),
  chromosome_bar_linewidth = 3.15,
  coordinate_contract = "true chromosome and OR midpoint positions on a shared Mb axis",
  status_track_encoding = list(functional = "above chromosome midline",
                               pseudogene = "below chromosome midline",
                               truncated = "full-height central tick"),
  legend_placement = "inside lower-right blank genomic area",
  chromosomes = nrow(chromosomes),
  chromosome_order = chromosome_order,
  total_OR_loci = nrow(loci),
  unique_global_ids = dplyr::n_distinct(loci$global_id),
  classification_counts = as.list(setNames(as.integer(class_counts), names(class_counts))),
  operational_clusters = nrow(clusters),
  large_clusters_ge_50 = nrow(large_clusters),
  focus_cluster_labels = as.list(setNames(focus_clusters$annotation, focus_clusters$cluster_id)),
  data_policy = "all Ppa OR loci shown; no sampling or status filtering",
  interpretation_boundary = "descriptive distribution only; no expansion, contraction, fixation, or adaptation claim"
)
write_json(qc, file.path(output_dir, "Ppa_OR_chromosome_distribution_separated_tracks.QC.json"),
           pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

cat("FIGURE_COMPLETED\n")
cat("output_prefix=", prefix, "\n", sep = "")
cat("Ppa_OR_loci=", nrow(loci), "\n", sep = "")
cat("functional=", class_counts[["functional"]], "\n", sep = "")
cat("pseudogene=", class_counts[["pseudogene"]], "\n", sep = "")
cat("truncated=", class_counts[["truncated"]], "\n", sep = "")
