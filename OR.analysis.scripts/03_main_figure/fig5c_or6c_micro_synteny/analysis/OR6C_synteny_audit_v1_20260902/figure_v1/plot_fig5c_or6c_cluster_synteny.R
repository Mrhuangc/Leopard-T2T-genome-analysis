suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript plot_fig5c_or6c_cluster_synteny.R <project_root> <out_check> <output_dir>")
}

project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_check <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
}

interval_path <- file.path(out_check, "extracted_intervals.tsv")
paf_path <- file.path(out_check, "Ppa_vs_Ppa_or.OR6C.asm5.detailed.v2.paf")
event_path <- file.path(project_root, "analysis", "OR6C_synteny_audit_v1_20260902",
                        "derived_v2", "OR6C.structural_events.min50bp.tsv")
locus_path <- file.path(project_root, "analysis", "OR_subfamily_comparative_v1_20260902",
                        "locus_OR_subfamily_genomic_context.tsv")
ppa_catalog_path <- file.path(project_root, "analysis", "OR_multispecies_20260825", "source",
                              "Ppa", "11_complete_OR_catalog", "Ppa.complete_OR_catalog.tsv")
ppa_or_catalog_path <- file.path(project_root, "analysis", "OR_multispecies_20260825", "source",
                                 "Ppa_or", "11_complete_OR_catalog", "Ppa_or.complete_OR_catalog.tsv")

required_files <- c(interval_path, paf_path, event_path, locus_path,
                    ppa_catalog_path, ppa_or_catalog_path)
missing_files <- required_files[!file.exists(required_files) | file.info(required_files)$size == 0]
if (length(missing_files) > 0L) {
  stop("Required input missing or empty: ", paste(basename(missing_files), collapse = ", "))
}

intervals <- read_tsv(interval_path) %>%
  mutate(
    analysis_unit = as.character(analysis_unit),
    extracted_start_1 = as.numeric(extracted_start_1),
    extracted_end_1 = as.numeric(extracted_end_1),
    extracted_bp = extracted_end_1 - extracted_start_1 + 1,
    species_label = if_else(analysis_unit == "Ppa", "P. p. japonensis", "P. p. orientalis"),
    track_y = if_else(analysis_unit == "Ppa", 1.25, 0.15)
  )
if (!setequal(intervals$analysis_unit, c("Ppa", "Ppa_or")) || nrow(intervals) != 2L) {
  stop("Interval table must contain exactly one Ppa and one Ppa_or row")
}

locus_map <- read_tsv(locus_path) %>%
  filter(analysis_unit %in% c("Ppa", "Ppa_or")) %>%
  select(analysis_unit, catalog_id, final_family, final_subfamily,
         assignment_confidence, representative_id)
if (anyDuplicated(locus_map[c("analysis_unit", "catalog_id")])) {
  stop("Functional subfamily map contains duplicate analysis_unit/catalog_id keys")
}

ppa_catalog <- read_tsv(ppa_catalog_path) %>% mutate(analysis_unit = "Ppa")
ppa_or_catalog <- read_tsv(ppa_or_catalog_path) %>% mutate(analysis_unit = "Ppa_or")
catalog <- bind_rows(ppa_catalog, ppa_or_catalog) %>%
  mutate(start_1 = as.numeric(start_1), end_1 = as.numeric(end_1)) %>%
  inner_join(intervals %>% select(analysis_unit, sseqid = contig, extracted_start_1,
                                  extracted_end_1, track_y, species_label),
             by = c("analysis_unit", "sseqid")) %>%
  filter(start_1 <= extracted_end_1, end_1 >= extracted_start_1) %>%
  left_join(locus_map, by = c("analysis_unit", "catalog_id")) %>%
  mutate(
    relative_mid_mb = (((start_1 + end_1) / 2) - extracted_start_1) / 1e6,
    display_class = case_when(
      classification == "functional" & final_subfamily == "OR6C" ~ "OR6C",
      classification == "functional" & final_subfamily == "OR9K" ~ "OR9K",
      classification == "functional" & final_subfamily == "OR2AP" ~ "OR2AP",
      classification == "functional" ~ "Other functional OR",
      classification == "pseudogene" ~ "Pseudogene",
      classification == "truncated" ~ "Truncated",
      TRUE ~ "Other catalog locus"
    ),
    locus_linetype = if_else(classification == "functional", "solid", "dashed")
  )

expected_cluster_counts <- c(Ppa = 65L, Ppa_or = 61L)
observed_counts <- table(catalog$analysis_unit)
if (any(observed_counts[names(expected_cluster_counts)] != expected_cluster_counts)) {
  stop("Extracted-region catalog counts disagree with audited cluster totals")
}

paf_lines <- readLines(paf_path, warn = FALSE)
paf_fields <- strsplit(paf_lines, "\t", fixed = TRUE)
if (length(paf_fields) == 0L || any(lengths(paf_fields) < 12L)) {
  stop("Detailed PAF is empty or malformed")
}
paf <- bind_rows(lapply(seq_along(paf_fields), function(i) {
  f <- paf_fields[[i]]
  data.frame(
    block_id = i,
    qstart0 = as.numeric(f[[3]]), qend0 = as.numeric(f[[4]]),
    strand = f[[5]], tstart0 = as.numeric(f[[8]]), tend0 = as.numeric(f[[9]]),
    mapq = as.numeric(f[[12]]), stringsAsFactors = FALSE
  )
})) %>%
  mutate(block_class = if_else(strand == "+", "Collinear alignment", "Reverse repeat match"))

ribbon_data <- bind_rows(lapply(seq_len(nrow(paf)), function(i) {
  x <- paf[i, ]
  target_order <- if (x$strand == "+") c(x$tend0, x$tstart0) else c(x$tstart0, x$tend0)
  data.frame(
    block_id = x$block_id,
    block_class = x$block_class,
    x = c(x$qstart0, x$qend0, target_order) / 1e6,
    y = c(1.17, 1.17, 0.23, 0.23)
  )
}))

events <- read_tsv(event_path) %>%
  mutate(displayed_in_figure = pmax(Ppa_length_bp, Ppa_or_length_bp, na.rm = TRUE) >= 1000)
major_events <- events %>%
  filter(displayed_in_figure) %>%
  mutate(
    event_label = case_when(
      event_type == "Ppa_insertion" ~ paste0("Ppa insertion +", round(Ppa_length_bp / 1000, 1), " kb"),
      event_type == "Ppa_or_insertion" ~ paste0("Ppa_or insertion +", round(Ppa_or_length_bp / 1000, 1), " kb"),
      event_type == "unaligned_chain_gap" ~ "Complex replacement",
      TRUE ~ event_type
    )
  )

event_rects <- bind_rows(
  major_events %>% filter(!is.na(Ppa_start_1), Ppa_length_bp > 0) %>%
    transmute(event_type, event_label, analysis_unit = "Ppa",
              xmin = (Ppa_start_1 - intervals$extracted_start_1[intervals$analysis_unit == "Ppa"]) / 1e6,
              xmax = (Ppa_end_1 - intervals$extracted_start_1[intervals$analysis_unit == "Ppa"] + 1) / 1e6,
              ymin = 1.08, ymax = 1.42),
  major_events %>% filter(!is.na(Ppa_or_start_1), Ppa_or_length_bp > 0) %>%
    transmute(event_type, event_label, analysis_unit = "Ppa_or",
              xmin = (Ppa_or_start_1 - intervals$extracted_start_1[intervals$analysis_unit == "Ppa_or"]) / 1e6,
              xmax = (Ppa_or_end_1 - intervals$extracted_start_1[intervals$analysis_unit == "Ppa_or"] + 1) / 1e6,
              ymin = -0.02, ymax = 0.32)
) %>%
  mutate(event_display = if_else(event_type == "unaligned_chain_gap",
                                 "Complex replacement", "Assembly-specific insertion"))

ppa_labels <- event_rects %>%
  filter(analysis_unit == "Ppa", event_type == "Ppa_insertion", (xmax - xmin) >= 0.01) %>%
  mutate(x = (xmin + xmax) / 2, y = 1.49,
         label = paste0("+", sprintf("%.1f", (xmax - xmin) * 1000), " kb"))

palette <- c(
  "OR6C" = "#D55E00",
  "OR9K" = "#0072B2",
  "OR2AP" = "#7A5195",
  "Other functional OR" = "#3B8F7A",
  "Pseudogene" = "#8A8A8A",
  "Truncated" = "#2B2B2B",
  "Other catalog locus" = "#BDBDBD"
)

plot <- ggplot() +
  geom_polygon(data = ribbon_data,
               aes(x = x, y = y, group = block_id, fill = block_class),
               colour = NA, alpha = 0.24) +
  scale_fill_manual(values = c("Collinear alignment" = "#9EC9C0",
                               "Reverse repeat match" = "#E6A15A"),
                    guide = "none") +
  geom_rect(data = event_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                colour = event_display),
            fill = NA, linewidth = 0.45,
            linetype = ifelse(event_rects$event_type == "unaligned_chain_gap", "dashed", "solid")) +
  scale_colour_manual(
    name = NULL,
    values = c(palette, "Assembly-specific insertion" = "#C33C54",
               "Complex replacement" = "#555555"),
    breaks = c("OR6C", "OR9K", "OR2AP", "Other functional OR",
               "Pseudogene", "Truncated", "Assembly-specific insertion", "Complex replacement")
  ) +
  geom_segment(data = intervals,
               aes(x = 0, xend = extracted_bp / 1e6, y = track_y, yend = track_y),
               linewidth = 0.35, colour = "#242424") +
  geom_segment(data = catalog,
               aes(x = relative_mid_mb, xend = relative_mid_mb,
                   y = track_y - 0.075, yend = track_y + 0.075,
                   colour = display_class, linetype = locus_linetype),
               linewidth = 0.75) +
  scale_linetype_identity() +
  geom_text(data = ppa_labels, aes(x = x, y = y, label = label),
            size = 2.25, colour = "#A62C45") +
  annotate("text", x = 0.015, y = 1.58, label = "italic('P. p. japonensis')", parse = TRUE,
           hjust = 0, size = 2.7) +
  annotate("text", x = 0.015, y = -0.13, label = "italic('P. p. orientalis')", parse = TRUE,
           hjust = 0, size = 2.7) +
  annotate("text", x = 0.91, y = 0.70, label = "reverse repeat match",
           size = 2.1, colour = "#A45C18", angle = -15) +
  scale_x_continuous(name = "Position within extracted homologous region (Mb)",
                     limits = c(0, max(intervals$extracted_bp) / 1e6),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(NULL, breaks = NULL) +
  coord_cartesian(ylim = c(-0.2, 1.66), clip = "off") +
  labs(title = "Local remodeling of a homologous gamma-OR cluster") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(linetype = "solid", linewidth = 1.2))) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 6.5),
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 5.8),
    legend.key.width = grid::unit(4.2, "mm"),
    plot.margin = margin(5, 5, 2, 5)
  )

write.table(catalog %>%
              select(analysis_unit, species_label, catalog_id, classification,
                     final_subfamily, assignment_confidence, sseqid, start_1,
                     end_1, strand, relative_mid_mb, display_class),
            file.path(output_dir, "OR_cluster_loci.source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(paf, file.path(output_dir, "OR_cluster_alignment_blocks.source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(events, file.path(output_dir, "OR_cluster_structural_events.source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

prefix <- file.path(output_dir, "Ppa_Ppaor_OR_cluster_synteny")
width_in <- 183 / 25.4
height_in <- 90 / 25.4
svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in)
print(plot)
dev.off()
grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in,
                     family = "Arial")
print(plot)
dev.off()
ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, background = "white")
print(plot)
dev.off()
ragg::agg_png(paste0(prefix, ".png"), width = width_in, height = height_in,
              units = "in", res = 300, background = "white")
print(plot)
dev.off()

qa <- data.frame(
  metric = c("Ppa_catalog_loci", "Ppa_or_catalog_loci", "alignment_blocks",
             "forward_blocks", "reverse_blocks", "major_structural_intervals",
             "display_event_threshold_bp"),
  value = c(sum(catalog$analysis_unit == "Ppa"), sum(catalog$analysis_unit == "Ppa_or"),
            nrow(paf), sum(paf$strand == "+"), sum(paf$strand == "-"),
            nrow(major_events), 1000)
)
write.table(qa, file.path(output_dir, "OR_cluster_synteny.figure_QC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("FIGURE_COMPLETED\n", "output_prefix=", prefix, "\n", sep = "")
