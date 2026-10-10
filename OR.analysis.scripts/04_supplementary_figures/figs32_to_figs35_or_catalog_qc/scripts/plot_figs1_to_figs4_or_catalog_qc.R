#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript plot_figs1_to_figs4_or_catalog_qc.R INPUT_SNAPSHOT OUTPUT_DIR", call. = FALSE)

input_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

required <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork", "ragg", "grid")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(grid)
})

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing input: ", path, call. = FALSE)
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
}

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

units <- c("Fca", "Lge", "Pbe", "Ple", "Ppa", "Ppa_or", "Pti")
species_text <- c(
  Fca = "F. catus", Lge = "L. geoffroyi", Pbe = "P. bengalensis",
  Ple = "P. leo", Ppa = "P. p. japonensis", Ppa_or = "P. p. orientalis",
  Pti = "P. t. sumatrae"
)
species_expr <- parse(text = c(
  "italic(F.~catus)", "italic(L.~geoffroyi)", "italic(P.~bengalensis)",
  "italic(P.~leo)", "italic(P.~p.~japonensis)", "italic(P.~p.~orientalis)",
  "italic(P.~t.~sumatrae)"
))
names(species_expr) <- units

class_cols <- c(functional = "#168C86", pseudogene = "#E6A21A", truncated = "#A73778", ambiguous = "#B8B8B8")
species_cols <- c(Fca = "#1B9E77", Lge = "#D95F02", Pbe = "#7570B3", Ple = "#E7298A", Ppa = "#E6AB02", Ppa_or = "#66A61E", Pti = "#A6761D")
confidence_cols <- c(high = "#275D7A", medium = "#68A7B4", low = "#E2AD50", unresolved = "#B6B6B6")

theme_or <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.32, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(size = base_size + 0.8, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#4D4D4D"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid.major.y = element_line(linewidth = 0.22, colour = "#E5E5E5"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_or())

save_png <- function(plot, filename, width_mm, height_mm, dpi = 300) {
  ragg::agg_png(filename, width = width_mm / 25.4, height = height_mm / 25.4, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
}

summary_path <- file.path(input_root, "analysis", "OR_multispecies_20260825", "master_catalog_v2", "OR_count_matrix.tsv")
summary <- read_tsv(summary_path) |>
  mutate(analysis_unit = factor(analysis_unit, levels = units))
assert(nrow(summary) == 7L && sum(summary$total) == 7303L, "OR summary did not reconcile")

s1_long <- summary |>
  select(analysis_unit, functional, pseudogene, truncated) |>
  pivot_longer(-analysis_unit, names_to = "classification", values_to = "loci") |>
  mutate(classification = factor(classification, levels = c("truncated", "pseudogene", "functional")))

s1 <- ggplot(s1_long, aes(analysis_unit, loci, fill = classification)) +
  geom_col(width = 0.72) +
  geom_text(data = summary, aes(analysis_unit, total, label = total), inherit.aes = FALSE, vjust = -0.35, size = 2.35, family = "Arial") +
  scale_fill_manual(values = class_cols, breaks = c("functional", "pseudogene", "truncated"), labels = c("Functional", "Pseudogene", "Truncated"), name = "Locus class") +
  scale_x_discrete(labels = species_expr[units]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of OR loci") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
save_png(s1, file.path(output_dir, "Supplementary_Figure_S1.png"), 183, 105)

cluster_path <- file.path(input_root, "figure", "OR_cluster_architecture_v1", "OR_operational_clusters.source_data.tsv")
cluster_summary_path <- file.path(input_root, "figure", "OR_cluster_architecture_v1", "OR_cluster_architecture.summary.tsv")
threshold_path <- file.path(input_root, "figure", "OR_cluster_threshold_sensitivity_v1", "OR_cluster_threshold_sensitivity.source_data.tsv")
clusters <- read_tsv(cluster_path) |> mutate(analysis_unit = factor(analysis_unit, levels = units))
cluster_summary <- read_tsv(cluster_summary_path) |> mutate(analysis_unit = factor(analysis_unit, levels = units))
threshold <- read_tsv(threshold_path) |> mutate(analysis_unit = factor(analysis_unit, levels = units))
assert(nrow(clusters) == 270L && sum(clusters$member_count) == 7303L, "Cluster table did not reconcile")
assert(nrow(threshold) == 42L, "Threshold grid is incomplete")

p2a <- ggplot(clusters, aes(size_rank, member_count, colour = analysis_unit, group = analysis_unit)) +
  geom_line(linewidth = 0.55) + geom_point(size = 0.85) +
  scale_colour_manual(values = species_cols, labels = species_expr[units], name = "Analysis unit") +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200), labels = label_number()) +
  labs(x = "Within-genome cluster-size rank", y = "OR loci per operational cluster") +
  theme(legend.position = "none")

largest_long <- cluster_summary |>
  select(analysis_unit, largest_cluster_members, largest_cluster_functional, largest_cluster_pseudogene, largest_cluster_truncated, largest_cluster_ambiguous) |>
  pivot_longer(starts_with("largest_cluster_"), names_to = "classification", values_to = "loci") |>
  filter(classification != "largest_cluster_members") |>
  mutate(classification = sub("largest_cluster_", "", classification), classification = factor(classification, levels = c("truncated", "pseudogene", "functional", "ambiguous")))
p2b <- ggplot(largest_long, aes(analysis_unit, loci, fill = classification)) +
  geom_col(width = 0.72) +
  geom_text(data = cluster_summary, aes(analysis_unit, largest_cluster_members, label = largest_cluster_members), inherit.aes = FALSE, vjust = -0.35, size = 2.0, family = "Arial") +
  scale_fill_manual(values = class_cols, breaks = c("functional", "pseudogene", "truncated"), labels = c("Functional", "Pseudogene", "Truncated"), name = "Locus class", guide = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = species_expr[units]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.09))) +
  labs(x = NULL, y = "OR loci in largest cluster") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

threshold_scale <- scale_x_log10(breaks = c(50, 100, 250, 500, 1000, 2000), labels = c("50", "100", "250", "500", "1,000", "2,000"))
p2c <- ggplot(threshold, aes(threshold_kb, fraction_loci_in_multi_locus_clusters, colour = analysis_unit, group = analysis_unit)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.15) + threshold_scale +
  scale_colour_manual(values = species_cols, labels = units, name = "Analysis unit", guide = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.94, 1.005)) +
  labs(x = "Maximum adjacent-locus gap (kb)", y = "Loci in multi-locus clusters") +
  theme(legend.position = "none")
p2d <- ggplot(threshold, aes(threshold_kb, top_two_cluster_fraction, colour = analysis_unit, group = analysis_unit)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.15) + threshold_scale +
  scale_colour_manual(values = species_cols, labels = units, name = "Analysis unit", guide = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.18, 0.405)) +
  labs(x = "Maximum adjacent-locus gap (kb)", y = "Loci in two largest clusters") +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5), legend.key.width = unit(8, "pt"), legend.spacing.x = unit(2, "pt"))
s2 <- wrap_plots(p2a, p2b, p2c, p2d, ncol = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))
save_png(s2, file.path(output_dir, "Supplementary_Figure_S2.png"), 183, 205)

completeness_path <- file.path(input_root, "analysis", "OR_subfamily_copy_number_20260904", "results", "OR_subfamily_assignment_completeness.tsv")
complete <- read_tsv(completeness_path) |>
  mutate(analysis_unit = factor(analysis_unit, levels = units))
assert(nrow(complete) == 7L && sum(complete$functional_loci) == 6419L, "Subfamily completeness table did not reconcile")
assert(sum(complete$primary_high_medium_named) == 4082L && sum(complete$all_named) == 5935L && sum(complete$unresolved_or_unnamed) == 484L, "Subfamily assignment totals changed")

confidence_long <- complete |>
  select(analysis_unit, high_confidence_named, medium_confidence_named, low_confidence_named, unresolved_or_unnamed) |>
  rename(high = high_confidence_named, medium = medium_confidence_named, low = low_confidence_named, unresolved = unresolved_or_unnamed) |>
  pivot_longer(-analysis_unit, names_to = "assignment_layer", values_to = "loci") |>
  mutate(assignment_layer = factor(assignment_layer, levels = c("unresolved", "low", "medium", "high")))
p3a <- ggplot(confidence_long, aes(analysis_unit, loci, fill = assignment_layer)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = confidence_cols, breaks = c("high", "medium", "low", "unresolved"), labels = c("High", "Medium", "Low", "Unresolved/unnamed"), name = "Assignment") +
  scale_x_discrete(labels = species_expr[units]) +
  labs(x = NULL, y = "Functional OR loci") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

complete_long <- complete |>
  select(analysis_unit, primary_completeness_fraction, all_named_completeness_fraction) |>
  pivot_longer(-analysis_unit, names_to = "set", values_to = "fraction") |>
  mutate(set = factor(set, levels = c("primary_completeness_fraction", "all_named_completeness_fraction"), labels = c("Primary: high + medium", "Sensitivity: all named")))
p3b <- ggplot(complete_long, aes(fraction, analysis_unit, colour = set)) +
  geom_line(aes(group = analysis_unit), colour = "#BDBDBD", linewidth = 0.5) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c("Primary: high + medium" = "#275D7A", "Sensitivity: all named" = "#E2AD50"), name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.55, 0.96)) +
  scale_y_discrete(labels = species_expr[units]) +
  labs(x = "Fraction of functional loci with named subfamily", y = NULL) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
s3 <- (p3a | p3b) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))
save_png(s3, file.path(output_dir, "Supplementary_Figure_S3.png"), 183, 112)

block_path <- file.path(input_root, "figure", "Ppa_OR_cluster_synteny_v1", "OR_cluster_alignment_blocks.source_data.tsv")
event_path <- file.path(input_root, "figure", "Ppa_OR_cluster_synteny_v1", "OR_cluster_structural_events.source_data.tsv")
pair_path <- file.path(input_root, "analysis", "OR6C_synteny_audit_v1_20260902", "derived_v2", "Ppa_to_Ppa_or.OR_locus_pairs.tsv")
blocks <- read_tsv(block_path)
events <- read_tsv(event_path) |> mutate(event_id = seq_len(n()), assembly = ifelse(event_type == "Ppa_insertion", "Ppa", "Ppa_or"), event_length_bp = pmax(Ppa_length_bp, Ppa_or_length_bp))
pairs <- read_tsv(pair_path)
assert(nrow(blocks) == 6L && nrow(events) == 14L && nrow(pairs) == 65L, "OR6C evidence counts changed")

p4a <- ggplot(blocks) +
  geom_abline(slope = 1, intercept = 0, colour = "#C7C7C7", linewidth = 0.35, linetype = 2) +
  geom_segment(aes(x = tstart0 / 1e6, xend = tend0 / 1e6, y = qstart0 / 1e6, yend = qend0 / 1e6, colour = block_class), linewidth = 0.8, lineend = "round") +
  scale_colour_manual(values = c("Collinear alignment" = "#4E88A8", "Reverse repeat match" = "#D98A3A"), name = NULL) +
  coord_equal() +
  labs(x = expression(italic("P. p. orientalis") * " interval (Mb)"), y = expression(italic("P. p. japonensis") * " interval (Mb)")) +
  theme(legend.position = "none")

p4b <- ggplot(events, aes(event_id, event_length_bp, fill = assembly, alpha = displayed_in_figure)) +
  geom_col(width = 0.72) +
  scale_y_log10(labels = label_number(big.mark = ",")) +
  scale_fill_manual(values = c(Ppa = "#C73E59", Ppa_or = "#4E88A8"), labels = c(Ppa = "Ppa-relative", Ppa_or = "Ppa_or-relative"), name = "Assembly-relative segment") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.38), labels = c(`TRUE` = "Shown in main panel", `FALSE` = "Recorded in source data"), name = "Main-panel display") +
  scale_x_continuous(breaks = seq_len(nrow(events))) +
  labs(x = "Structural-event record", y = "Event length (bp, log scale)") +
  theme(legend.position = "none")

aligned_pairs <- pairs |> filter(mapping_status == "aligned", !is.na(mapped_Ppa_or_position_1)) |>
  mutate(display_group = case_when(Ppa_classification == "pseudogene" ~ "Pseudogene", Ppa_subfamily %in% c("OR6C", "OR9K", "OR2AP") ~ Ppa_subfamily, TRUE ~ "Other functional OR"))
assert(nrow(aligned_pairs) == 60L, "Expected 60 aligned positional pairs")
p4c <- ggplot(aligned_pairs, aes(mapped_Ppa_or_position_1 / 1e6, (Ppa_start_1 + Ppa_end_1) / 2e6, colour = display_group, shape = reciprocal_nearest_locus)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = c(OR6C = "#D55E00", OR9K = "#0072B2", OR2AP = "#7A5195", `Other functional OR` = "#2A9D8F", Pseudogene = "#8C8C8C"), name = "Ppa locus class") +
  scale_shape_manual(values = c(true = 16, false = 1), labels = c(true = "Reciprocal nearest", false = "Not reciprocal"), name = NULL) +
  labs(x = expression("Mapped position in " * italic("P. p. orientalis") * " (Mb)"), y = expression("Position in " * italic("P. p. japonensis") * " (Mb)")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15, label = "60 aligned pairs; 5 insertion-associated or unaligned", size = 2.1, family = "Arial") +
  theme(legend.position = "none")
s4 <- (p4a | p4b | p4c) + plot_layout(widths = c(0.82, 1, 1.18)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))
save_png(s4, file.path(output_dir, "Supplementary_Figure_S4.png"), 183, 96)

outputs <- file.path(output_dir, paste0("Supplementary_Figure_S", 1:4, ".png"))
assert(all(file.exists(outputs) & file.info(outputs)$size > 0), "One or more supplementary figures were not written")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_session_info.txt"))
cat(paste(outputs, collapse = "\n"), "\n")
