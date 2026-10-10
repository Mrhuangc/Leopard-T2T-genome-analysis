#!/usr/bin/env Rscript

# ============================================================
# Visualize high-divergence centromeric satellite array results
# Compatible with satellite_array_pipeline.optimized.py outputs.
# It can also read the earlier pipeline outputs; region-level and
# full-length-only panels are skipped when the corresponding files/columns
# do not exist.
#
# Usage:
#   Rscript plot_SAT113_array_results.R \
#     PPA.SAT113.array_pipeline/07_results \
#     PPA.SAT113 \
#     PPA.SAT113.array_pipeline/08_plots \
#     60
#
# Arguments:
#   1. result directory
#   2. prefix
#   3. output plot directory
#   4. maximum number of arrays in heatmap/PCoA (default 60)
# ============================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr", "patchwork", "scales"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install with:\n",
    "  install.packages(c(",
    paste(sprintf("'%s'", missing_pkgs), collapse = ", "), "))\n",
    "or use conda/mamba: r-ggplot2 r-dplyr r-tidyr r-readr r-stringr r-patchwork r-scales"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args) >= 1) args[[1]] else "PPA.SAT113.array_pipeline/07_results"
prefix <- if (length(args) >= 2) args[[2]] else "PPA.SAT113"
out_dir <- if (length(args) >= 3) args[[3]] else file.path(dirname(result_dir), "08_plots")
heatmap_top_n <- if (length(args) >= 4) as.integer(args[[4]]) else 60L
if (is.na(heatmap_top_n) || heatmap_top_n < 5) heatmap_top_n <- 60L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_of <- function(suffix) file.path(result_dir, paste0(prefix, suffix))
read_tsv_safe <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Required file not found: ", path)
    message("[skip] Optional file not found: ", path)
    return(NULL)
  }
  readr::read_tsv(path, na = c("NA", "NaN", "nan", ""), show_col_types = FALSE, progress = FALSE)
}

copy_file <- path_of(".copy_divergence.tsv")
array_file <- path_of(".array_divergence_summary.tsv")
array_bins_file <- path_of(".array_position_bins.tsv")
region_file <- path_of(".region_divergence_summary.tsv")
region_bins_file <- path_of(".region_position_bins.tsv")
cons_stats_file <- path_of(".final.consensus_stats.tsv")
matrix_file <- path_of(".inter_array_substitution_divergence.matrix.tsv")
run_summary_file <- path_of(".run_summary.tsv")

copy_df <- read_tsv_safe(copy_file)
array_df <- read_tsv_safe(array_file)
array_bins <- read_tsv_safe(array_bins_file)
region_df <- read_tsv_safe(region_file, required = FALSE)
region_bins <- read_tsv_safe(region_bins_file, required = FALSE)
cons_stats <- read_tsv_safe(cons_stats_file, required = FALSE)
run_summary <- read_tsv_safe(run_summary_file, required = FALSE)

# -------------------------------
# Helpers
# -------------------------------
first_existing <- function(df, candidates, required = TRUE) {
  found <- candidates[candidates %in% names(df)]
  if (length(found) > 0) return(found[[1]])
  if (required) stop("None of the expected columns exist: ", paste(candidates, collapse = ", "))
  return(NA_character_)
}

natural_chrom_levels <- function(x) {
  x <- unique(as.character(x))
  m <- stringr::str_match(x, "^Chr([A-Za-z]*)([0-9]+|X|Y)$")
  prefix_part <- ifelse(is.na(m[, 2]), x, m[, 2])
  number_part <- suppressWarnings(as.numeric(m[, 3]))
  number_part[m[, 3] == "X"] <- 1e6
  number_part[m[, 3] == "Y"] <- 1e6 + 1
  number_part[is.na(number_part)] <- 2e6
  x[order(prefix_part, number_part, x)]
}

set_chrom_factor <- function(df) {
  if (!"chrom" %in% names(df)) return(df)
  df$chrom <- factor(df$chrom, levels = natural_chrom_levels(df$chrom))
  df
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_q <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, names = FALSE, type = 7))
}

save_plot <- function(plot, stem, width = 10, height = 7) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", device = "pdf")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 300)
  message("[write] ", pdf_path)
}

theme_sat <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey60"),
      strip.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = rel(1.08)),
      plot.subtitle = element_text(color = "grey30"),
      legend.position = "right"
    )
}

copy_df <- set_chrom_factor(copy_df)
array_df <- set_chrom_factor(array_df)
if (!is.null(region_df)) region_df <- set_chrom_factor(region_df)
if (!is.null(array_bins)) array_bins <- set_chrom_factor(array_bins)
if (!is.null(region_bins)) region_bins <- set_chrom_factor(region_bins)

# Explicit metric columns; remain compatible with v1 results.
copy_sub_col <- first_existing(copy_df, c("substitution_divergence_pct"))
array_copy_all_col <- first_existing(array_df, c("copy_to_consensus_sub_div_median_all_pct"))
array_copy_full_col <- first_existing(array_df, c("copy_to_consensus_sub_div_median_full_length_pct"), required = FALSE)
pair_all_col <- first_existing(
  array_df,
  c("intra_array_pairwise_sub_div_median_all_pct", "intra_array_pairwise_sub_div_median_pct")
)
pair_full_col <- first_existing(
  array_df,
  c("intra_array_pairwise_sub_div_median_full_length_pct"),
  required = FALSE
)

# -------------------------------
# Summary tables
# -------------------------------
chrom_summary <- copy_df %>%
  group_by(chrom) %>%
  summarise(
    n_candidates = n(),
    n_full_length = sum(quality_tier == "full_length", na.rm = TRUE),
    full_length_fraction = mean(quality_tier == "full_length", na.rm = TRUE),
    median_copy_divergence_all_pct = safe_median(.data[[copy_sub_col]]),
    median_copy_divergence_full_length_pct = safe_median(
      .data[[copy_sub_col]][quality_tier == "full_length"]
    ),
    .groups = "drop"
  )
readr::write_tsv(chrom_summary, file.path(out_dir, paste0(prefix, ".chromosome_summary.tsv")))

edge_center_summary <- function(df, group_col, score_col, min_copies = 10) {
  if (!all(c(group_col, score_col, "quality_tier", copy_sub_col) %in% names(df))) return(NULL)
  df %>%
    filter(quality_tier == "full_length") %>%
    mutate(zone = case_when(
      .data[[score_col]] <= 0.25 ~ "edge",
      .data[[score_col]] >= 0.75 ~ "center",
      TRUE ~ "middle"
    )) %>%
    filter(zone %in% c("edge", "center")) %>%
    group_by(.data[[group_col]], chrom, zone) %>%
    summarise(n = n(), median_divergence_pct = safe_median(.data[[copy_sub_col]]), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = zone, values_from = c(n, median_divergence_pct)) %>%
    filter(coalesce(n_edge, 0) >= min_copies, coalesce(n_center, 0) >= min_copies) %>%
    mutate(edge_minus_center_pct = median_divergence_pct_edge - median_divergence_pct_center)
}

array_edge_center <- edge_center_summary(copy_df, "array_key", "edge_to_center_rank_score", min_copies = 3)
if (!is.null(array_edge_center)) {
  write_tsv(array_edge_center, file.path(out_dir, paste0(prefix, ".array_edge_center.tsv")))
}
region_edge_center <- NULL
if ("region_key" %in% names(copy_df) && "region_edge_to_center_rank_score" %in% names(copy_df)) {
  region_edge_center <- edge_center_summary(
    copy_df, "region_key", "region_edge_to_center_rank_score", min_copies = 10
  )
  if (!is.null(region_edge_center)) {
    write_tsv(region_edge_center, file.path(out_dir, paste0(prefix, ".region_edge_center.tsv")))
  }
}

# -------------------------------
# Figure 1: overview
# -------------------------------
p_chr_count <- ggplot(chrom_summary, aes(x = chrom, y = n_candidates)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = scales::comma(n_candidates)), vjust = -0.25, size = 3) +
  scale_y_continuous(trans = scales::log1p_trans(), labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "SAT113 candidate copy number",
    subtitle = "The y-axis uses log(1+x) transformation",
    x = NULL, y = "Candidate copies"
  ) +
  theme_sat() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_chr_full <- ggplot(chrom_summary, aes(x = chrom, y = full_length_fraction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Near-full-length monomer fraction", x = NULL, y = "Full-length fraction") +
  theme_sat() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

copy_plot_df <- copy_df %>% filter(is.finite(.data[[copy_sub_col]]))
p_chr_div <- ggplot(copy_plot_df, aes(x = chrom, y = .data[[copy_sub_col]])) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3) +
  coord_cartesian(ylim = c(0, min(60, safe_q(copy_plot_df[[copy_sub_col]], 0.995) + 2))) +
  labs(title = "Copy-to-local-consensus divergence", x = NULL, y = "Substitution divergence (%)") +
  theme_sat() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cons_len <- ggplot(array_df, aes(x = local_consensus_length)) +
  geom_histogram(binwidth = 1, boundary = 0.5, color = "white") +
  geom_vline(xintercept = 113, linetype = 2, linewidth = 0.5) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(title = "Array-specific consensus length", x = "Consensus length (bp)", y = "Array blocks") +
  theme_sat()

overview_plot <- (p_chr_count | p_chr_full) / (p_chr_div | p_cons_len) +
  plot_annotation(title = paste0(prefix, " array-specific SAT analysis"))
save_plot(overview_plot, paste0(prefix, ".01_overview"), width = 12, height = 9)

# -------------------------------
# Figure 2: array-level divergence
# -------------------------------
array_scatter <- array_df %>%
  mutate(label_flag = n_candidates >= safe_q(n_candidates, 0.95) |
           .data[[array_copy_all_col]] >= safe_q(.data[[array_copy_all_col]], 0.95) |
           .data[[array_copy_all_col]] <= safe_q(.data[[array_copy_all_col]], 0.05))

p_array_size <- ggplot(
  array_scatter,
  aes(x = n_candidates, y = .data[[array_copy_all_col]], size = full_length_fraction, color = chrom)
) +
  geom_point(alpha = 0.8) +
  geom_text(
    data = filter(array_scatter, label_flag), aes(label = array_key),
    size = 2.4, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_x_continuous(trans = "log10", labels = scales::comma) +
  scale_size_continuous(range = c(1.5, 6), labels = percent_format(accuracy = 1)) +
  labs(
    title = "Array size and local sequence divergence",
    x = "Candidate copies per TRASH2 block (log10)",
    y = "Median copy-to-consensus divergence (%)",
    size = "Full length", color = "Chromosome"
  ) + theme_sat()

if (!is.na(pair_full_col)) {
  pair_compare_df <- array_df %>%
    filter(is.finite(.data[[pair_all_col]]), is.finite(.data[[pair_full_col]])) %>%
    mutate(fragment_effect = .data[[pair_all_col]] - .data[[pair_full_col]])
  max_pair <- max(c(pair_compare_df[[pair_all_col]], pair_compare_df[[pair_full_col]]), na.rm = TRUE)
  p_pair_compare <- ggplot(
    pair_compare_df,
    aes(x = .data[[pair_full_col]], y = .data[[pair_all_col]], color = chrom, size = n_candidates)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) +
    geom_point(alpha = 0.8) +
    coord_equal(xlim = c(0, max_pair * 1.05), ylim = c(0, max_pair * 1.05)) +
    scale_size_continuous(range = c(1.5, 6), trans = "sqrt") +
    labs(
      title = "Effect of fragments on intra-array pairwise divergence",
      subtitle = "Points above the diagonal are inflated when all copies are included",
      x = "Full-length-only pairwise divergence (%)",
      y = "All-copy pairwise divergence (%)", color = "Chromosome", size = "Copies"
    ) + theme_sat()
} else {
  p_pair_compare <- ggplot(array_df, aes(x = n_candidates, y = .data[[pair_all_col]], color = chrom)) +
    geom_point(alpha = 0.8) +
    scale_x_continuous(trans = "log10", labels = comma) +
    labs(
      title = "Intra-array pairwise divergence",
      subtitle = "The current file contains the all-copy metric only",
      x = "Candidate copies per block (log10)", y = "Pairwise divergence (%)", color = "Chromosome"
    ) + theme_sat()
}

p_array_copy_pair <- ggplot(
  array_df,
  aes(x = .data[[array_copy_all_col]], y = .data[[pair_all_col]], color = chrom, size = n_candidates)
) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(1.5, 6), trans = "sqrt") +
  labs(
    title = "Consensus-based versus pairwise divergence",
    x = "Median copy-to-consensus divergence (%)",
    y = "Median intra-array pairwise divergence (%)",
    color = "Chromosome", size = "Copies"
  ) + theme_sat()

p_top <- array_df %>%
  slice_max(order_by = n_candidates, n = 20, with_ties = FALSE) %>%
  mutate(array_key = reorder(array_key, n_candidates)) %>%
  ggplot(aes(x = array_key, y = n_candidates, fill = chrom)) +
  geom_col() + coord_flip() +
  labs(title = "Largest SAT113 blocks", x = NULL, y = "Candidate copies", fill = "Chromosome") +
  theme_sat()

array_panel <- (p_array_size | p_pair_compare) / (p_array_copy_pair | p_top)
save_plot(array_panel, paste0(prefix, ".02_array_divergence"), width = 13, height = 10)

# -------------------------------
# Position-gradient plots
# -------------------------------
standardize_bin_columns <- function(df) {
  all_col <- first_existing(
    df,
    c("median_substitution_divergence_all_pct", "median_substitution_divergence_pct")
  )
  full_col <- first_existing(
    df,
    c("median_substitution_divergence_full_length_pct"),
    required = FALSE
  )
  df %>%
    mutate(
      bin_mid = (bin_start_frac + bin_end_frac) / 2,
      divergence_all = .data[[all_col]],
      divergence_full = if (!is.na(full_col)) .data[[full_col]] else NA_real_
    )
}

aggregate_bins <- function(df, value_col) {
  df %>%
    filter(is.finite(.data[[value_col]])) %>%
    group_by(bin_index, bin_mid) %>%
    summarise(
      n_groups = n(),
      median = safe_median(.data[[value_col]]),
      q25 = safe_q(.data[[value_col]], 0.25),
      q75 = safe_q(.data[[value_col]], 0.75),
      .groups = "drop"
    )
}

plot_position_gradient <- function(bin_df, title_prefix, stem) {
  std <- standardize_bin_columns(bin_df)
  all_agg <- aggregate_bins(std, "divergence_all") %>% mutate(dataset = "All candidates")
  full_agg <- aggregate_bins(std, "divergence_full") %>% mutate(dataset = "Full-length only")
  agg <- bind_rows(all_agg, full_agg) %>% filter(is.finite(median))

  p_global <- ggplot(agg, aes(x = bin_mid, y = median, group = dataset, linetype = dataset)) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = dataset), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = paste0(title_prefix, ": pooled positional gradient"),
      subtitle = "Line: median across groups; ribbon: inter-group IQR",
      x = "Relative position from region/block start", y = "Median substitution divergence (%)",
      linetype = NULL, fill = NULL
    ) + theme_sat()

  chrom_agg <- std %>%
    filter(is.finite(divergence_all)) %>%
    group_by(chrom, bin_index, bin_mid) %>%
    summarise(median = safe_median(divergence_all), n_groups = n(), .groups = "drop")
  p_chrom <- ggplot(chrom_agg, aes(x = bin_mid, y = median, group = chrom)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
    facet_wrap(~chrom, scales = "free_y", ncol = 4) +
    scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1)) +
    labs(
      title = paste0(title_prefix, ": chromosome-stratified gradient"),
      x = "Relative position", y = "Median divergence (%)"
    ) + theme_sat(base_size = 9) + theme(legend.position = "none")

  save_plot(p_global / p_chrom, stem, width = 12, height = 10)
}

plot_position_gradient(
  array_bins, "TRASH2 block", paste0(prefix, ".03_array_position_gradient")
)
if (!is.null(region_bins)) {
  plot_position_gradient(
    region_bins, "Merged SAT-rich region", paste0(prefix, ".04_region_position_gradient")
  )
}

# Edge-center plot
edge_plot_data <- bind_rows(
  if (!is.null(array_edge_center)) mutate(array_edge_center, level = "TRASH2 block") else NULL,
  if (!is.null(region_edge_center)) mutate(region_edge_center, level = "Merged region") else NULL
)
if (nrow(edge_plot_data) > 0) {
  p_edge <- ggplot(edge_plot_data, aes(x = chrom, y = edge_minus_center_pct, fill = level)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA, width = 0.65) +
    geom_point(
      aes(group = level), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      alpha = 0.55, size = 1.4
    ) +
    labs(
      title = "Edge-to-center divergence contrast",
      subtitle = "Positive values indicate more divergent edges",
      x = NULL, y = "Median divergence at edge − center (percentage points)", fill = NULL
    ) + theme_sat() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_edge, paste0(prefix, ".05_edge_center_contrast"), width = 11, height = 6.5)
}

# -------------------------------
# Corrected inter-array matrix: clustered heatmap and PCoA
# -------------------------------
if (file.exists(matrix_file)) {
  matrix_df <- read_tsv_safe(matrix_file)
  array_name_col <- names(matrix_df)[1]
  mat <- as.matrix(matrix_df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- matrix_df[[array_name_col]]
  colnames(mat) <- names(matrix_df)[-1]

  reliable <- array_df %>%
    mutate(consensus_N_fraction = NA_real_)
  if (!is.null(cons_stats) && all(c("array_key", "consensus_length", "consensus_N_columns") %in% names(cons_stats))) {
    reliable <- reliable %>%
      left_join(
        cons_stats %>% transmute(
          array_key,
          consensus_N_fraction = consensus_N_columns / pmax(consensus_length, 1)
        ),
        by = "array_key", suffix = c("", ".new")
      ) %>%
      mutate(consensus_N_fraction = coalesce(consensus_N_fraction.new, consensus_N_fraction)) %>%
      select(-consensus_N_fraction.new)
  }
  reliable <- reliable %>%
    filter(n_candidates >= 10, full_length_fraction >= 0.60) %>%
    filter(is.na(consensus_N_fraction) | consensus_N_fraction <= 0.20) %>%
    arrange(desc(n_candidates)) %>%
    slice_head(n = heatmap_top_n)

  ids <- intersect(reliable$array_key, intersect(rownames(mat), colnames(mat)))
  if (length(ids) >= 3) {
    mat2 <- mat[ids, ids, drop = FALSE]
    # Symmetrize and fill rare missing distances only for clustering.
    mat2 <- (mat2 + t(mat2)) / 2
    diag(mat2) <- 0
    finite_values <- mat2[is.finite(mat2)]
    fill_value <- if (length(finite_values) > 0) median(finite_values) else 0
    mat2[!is.finite(mat2)] <- fill_value

    hc <- hclust(as.dist(mat2), method = "average")
    ordered_ids <- hc$labels[hc$order]
    heat_long <- as.data.frame(as.table(mat2), stringsAsFactors = FALSE) %>%
      rename(array1 = Var1, array2 = Var2, divergence = Freq) %>%
      mutate(
        array1 = factor(array1, levels = ordered_ids),
        array2 = factor(array2, levels = rev(ordered_ids))
      )

    p_heat <- ggplot(heat_long, aes(x = array1, y = array2, fill = divergence)) +
      geom_tile() +
      scale_fill_viridis_c(option = "C", name = "Divergence\n(%)", na.value = "grey90") +
      coord_fixed() +
      labs(
        title = "Globally phase-normalized array-consensus divergence",
        subtitle = paste0("Average-linkage clustering; ", length(ids), " reliable/high-copy arrays"),
        x = NULL, y = NULL
      ) +
      theme_sat(base_size = 8) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
        axis.text.y = element_text(size = 5.5)
      )
    save_plot(p_heat, paste0(prefix, ".06_inter_array_heatmap"), width = 10, height = 9)

    pcoa <- cmdscale(as.dist(mat2), k = 2, eig = TRUE)
    pcoa_df <- tibble(
      array_key = rownames(pcoa$points),
      PCo1 = pcoa$points[, 1],
      PCo2 = pcoa$points[, 2]
    ) %>% left_join(array_df %>% select(array_key, chrom, n_candidates, full_length_fraction), by = "array_key")
    eig <- pcoa$eig[pcoa$eig > 0]
    explained <- if (length(eig) >= 2) 100 * eig[1:2] / sum(eig) else c(NA_real_, NA_real_)

    p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = chrom, size = n_candidates)) +
      geom_point(alpha = 0.85) +
      geom_text(
        data = pcoa_df %>% slice_max(abs(PCo1) + abs(PCo2), n = 10, with_ties = FALSE),
        aes(label = array_key), size = 2.5, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE
      ) +
      scale_size_continuous(range = c(1.8, 7), trans = "sqrt") +
      labs(
        title = "PCoA of array-consensus divergence",
        x = ifelse(is.finite(explained[1]), sprintf("PCo1 (%.1f%%)", explained[1]), "PCo1"),
        y = ifelse(is.finite(explained[2]), sprintf("PCo2 (%.1f%%)", explained[2]), "PCo2"),
        color = "Chromosome", size = "Copies"
      ) + theme_sat()
    save_plot(p_pcoa, paste0(prefix, ".07_inter_array_PCoA"), width = 9, height = 7)
  } else {
    message("[skip] Fewer than three reliable arrays remain for heatmap/PCoA.")
  }
}

# -------------------------------
# Compact statistics text
# -------------------------------
stats_lines <- c(
  paste0("prefix\t", prefix),
  paste0("candidate_copies\t", nrow(copy_df)),
  paste0("array_blocks\t", nrow(array_df)),
  paste0("full_length_fraction\t", sprintf("%.6f", mean(copy_df$quality_tier == "full_length", na.rm = TRUE))),
  paste0("median_copy_divergence_all_pct\t", sprintf("%.6f", safe_median(copy_df[[copy_sub_col]]))),
  paste0("median_array_copy_divergence_pct\t", sprintf("%.6f", safe_median(array_df[[array_copy_all_col]]))),
  paste0("median_array_pairwise_divergence_all_pct\t", sprintf("%.6f", safe_median(array_df[[pair_all_col]])))
)
if (!is.na(pair_full_col)) {
  stats_lines <- c(
    stats_lines,
    paste0(
      "median_array_pairwise_divergence_full_length_pct\t",
      sprintf("%.6f", safe_median(array_df[[pair_full_col]]))
    )
  )
}
if (!is.null(region_df)) stats_lines <- c(stats_lines, paste0("merged_regions\t", nrow(region_df)))
writeLines(stats_lines, file.path(out_dir, paste0(prefix, ".plot_summary.txt")))

message("[done] Plots and summary tables written to: ", normalizePath(out_dir, mustWork = FALSE))
