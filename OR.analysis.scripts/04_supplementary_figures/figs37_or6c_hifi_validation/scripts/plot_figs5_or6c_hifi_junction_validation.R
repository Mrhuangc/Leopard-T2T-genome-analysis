#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript plot_figs5_or6c_hifi_junction_validation.R <analysis_directory>")
}

analysis_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
input_dir <- file.path(analysis_dir, "input_snapshot")
figure_dir <- file.path(analysis_dir, "figure")
table_dir <- file.path(analysis_dir, "table")
qa_dir <- file.path(analysis_dir, "qa")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

events <- data.frame(
  event_id = c("OR6C_INS27kb", "OR6C_INS17kb"),
  event_label = c("27.3-kb segment", "17.3-kb segment"),
  panel_label = c("a", "b"),
  breakpoint = c(80631604L, 81020332L),
  depth_file = c(
    "OR6C_INS27kb.MAPQ20.depth.tsv.gz",
    "OR6C_INS17kb.MAPQ20.depth.tsv.gz"
  ),
  read_file = c(
    "OR6C_INS27kb.direct_primary_MAPQ20.tsv",
    "OR6C_INS17kb.direct_primary_MAPQ20.tsv"
  ),
  stringsAsFactors = FALSE
)

required_files <- c(
  events$depth_file,
  events$read_file,
  "event_evidence_summary.tsv",
  "coverage_summary.tsv"
)
missing_files <- required_files[!file.exists(file.path(input_dir, required_files))]
if (length(missing_files) > 0L) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "))
}

evidence <- read.delim(
  file.path(input_dir, "event_evidence_summary.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
coverage_windows <- read.delim(
  file.path(input_dir, "coverage_summary.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

depth_list <- vector("list", nrow(events))
read_list <- vector("list", nrow(events))

for (i in seq_len(nrow(events))) {
  depth_i <- read.delim(
    gzfile(file.path(input_dir, events$depth_file[[i]])),
    header = FALSE,
    col.names = c("contig", "position_1based", "depth"),
    stringsAsFactors = FALSE
  )
  stopifnot(nrow(depth_i) == 100001L)
  stopifnot(length(unique(depth_i$contig)) == 1L)
  stopifnot(all(diff(depth_i$position_1based) == 1L))
  depth_i$event_id <- events$event_id[[i]]
  depth_i$position_relative_bp <- depth_i$position_1based - events$breakpoint[[i]]
  depth_i$position_relative_kb <- depth_i$position_relative_bp / 1000
  depth_i <- depth_i[
    depth_i$position_relative_bp >= -12000L &
      depth_i$position_relative_bp <= 12000L,
  ]
  stopifnot(nrow(depth_i) == 24001L)
  stopifnot(sum(depth_i$depth == 0L) == 0L)
  depth_list[[i]] <- depth_i

  reads_i <- read.delim(
    file.path(input_dir, events$read_file[[i]]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(nrow(reads_i) > 0L)
  stopifnot(length(unique(reads_i$read_name)) == nrow(reads_i))
  stopifnot(all(reads_i$MAPQ >= 20L))
  reads_i$event_id <- events$event_id[[i]]
  reads_i$strand <- ifelse(bitwAnd(as.integer(reads_i$flag), 16L) > 0L, "Reverse", "Forward")
  reads_i$start_relative_kb <- (reads_i$alignment_start_1based - events$breakpoint[[i]]) / 1000
  reads_i$end_relative_kb <- (reads_i$alignment_end_1based - events$breakpoint[[i]]) / 1000
  reads_i <- reads_i[order(reads_i$alignment_start_1based, reads_i$alignment_end_1based), ]
  reads_i$display_row <- seq_len(nrow(reads_i))
  reads_i$display_start_kb <- pmax(reads_i$start_relative_kb, -12)
  reads_i$display_end_kb <- pmin(reads_i$end_relative_kb, 12)
  read_list[[i]] <- reads_i
}

depth_plot_data <- bind_rows(depth_list)
read_plot_data <- bind_rows(read_list)

depth_binned <- depth_plot_data %>%
  mutate(bin_center_bp = floor((position_relative_bp + 50) / 100) * 100) %>%
  group_by(event_id, bin_center_bp) %>%
  summarise(
    position_relative_kb = first(bin_center_bp) / 1000,
    mean_depth = mean(depth),
    .groups = "drop"
  )

write.table(
  depth_binned,
  file.path(table_dir, "Source_data_depth_100bp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  read_plot_data[, c(
    "event_id", "read_name", "strand", "MAPQ", "CIGAR",
    "alignment_start_1based", "alignment_end_1based",
    "start_relative_kb", "end_relative_kb", "read_length_bp",
    "soft_clipped_bp", "NM"
  )],
  file.path(table_dir, "Source_data_spanning_reads.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

palette <- c(
  neutral_light = "#D8D8D8",
  neutral_mid = "#737373",
  neutral_dark = "#2F2F2F",
  teal = "#168C88",
  teal_light = "#73BDB8",
  blue = "#4C6F91",
  breakpoint = "#C33D5A",
  center_fill = "#F5DDE3"
)

theme_pub <- function(base_size = 7.2, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(5, "mm"),
      legend.key.height = grid::unit(2.4, "mm"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.margin = margin(3, 4, 2, 4),
      panel.grid = element_blank()
    )
}

coverage_max <- max(depth_plot_data$depth)
y_limit <- ceiling(coverage_max / 10) * 10

make_event_column <- function(event_row, depth_raw, depth_mean, reads, evidence_row) {
  n_reads <- nrow(reads)
  unique_starts <- length(unique(reads$alignment_start_1based))
  stopifnot(n_reads == evidence_row$MAPQ20_direct_spanning_reads)
  stopifnot(unique_starts == evidence_row$MAPQ20_unique_alignment_starts)

  annotation_text <- sprintf(
    "n = %d; MAPQ = %d; center/flank = %.2f",
    n_reads,
    min(reads$MAPQ),
    evidence_row$center_to_flank_depth_ratio_MAPQ20
  )

  p_depth <- ggplot() +
    annotate(
      "rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = Inf,
      fill = palette[["center_fill"]], alpha = 0.55
    ) +
    geom_line(
      data = depth_raw,
      aes(x = position_relative_kb, y = depth),
      colour = palette[["neutral_light"]],
      linewidth = 0.25
    ) +
    geom_line(
      data = depth_mean,
      aes(x = position_relative_kb, y = mean_depth),
      colour = palette[["teal"]],
      linewidth = 0.65
    ) +
    geom_vline(
      xintercept = 0,
      colour = palette[["breakpoint"]],
      linetype = "22",
      linewidth = 0.55
    ) +
    annotate(
      "text", x = -11.6, y = y_limit * 0.95,
      label = annotation_text,
      hjust = 0, vjust = 1,
      family = "Arial", size = 2.45,
      colour = palette[["neutral_dark"]]
    ) +
    scale_x_continuous(
      limits = c(-12, 12),
      breaks = c(-12, -6, 0, 6, 12),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, y_limit),
      breaks = seq(0, y_limit, by = 20),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = paste0(event_row$panel_label, "  ", event_row$event_label),
      x = NULL,
      y = "HiFi depth"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 8, face = "bold")
    )

  p_reads <- ggplot(reads) +
    geom_vline(
      xintercept = 0,
      colour = palette[["breakpoint"]],
      linetype = "22",
      linewidth = 0.55
    ) +
    geom_segment(
      aes(
        x = display_start_kb,
        xend = display_end_kb,
        y = display_row,
        yend = display_row,
        colour = strand
      ),
      linewidth = 0.62,
      lineend = "round",
      alpha = 0.94
    ) +
    scale_colour_manual(
      values = c("Forward" = palette[["teal"]], "Reverse" = palette[["blue"]]),
      breaks = c("Forward", "Reverse")
    ) +
    scale_x_continuous(
      limits = c(-12, 12),
      breaks = c(-12, -6, 0, 6, 12),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0.5, n_reads + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Position relative to junction (kb)",
      y = "Spanning HiFi reads",
      colour = NULL
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom"
    )

  p_depth / p_reads + plot_layout(heights = c(0.43, 1))
}

columns <- vector("list", nrow(events))
for (i in seq_len(nrow(events))) {
  event_id <- events$event_id[[i]]
  columns[[i]] <- make_event_column(
    events[i, ],
    depth_plot_data[depth_plot_data$event_id == event_id, ],
    depth_binned[depth_binned$event_id == event_id, ],
    read_plot_data[read_plot_data$event_id == event_id, ],
    evidence[evidence$event_id == event_id, ]
  )
}

figure <- (columns[[1]] | columns[[2]]) +
  plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "bottom")

output_stem <- file.path(figure_dir, "Supplementary_Figure_OR6C_HiFi_junction_validation")
width_mm <- 183
height_mm <- 122
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(
  paste0(output_stem, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(figure)
dev.off()

grDevices::cairo_pdf(
  paste0(output_stem, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Arial",
  bg = "white"
)
print(figure)
dev.off()

ragg::agg_png(
  paste0(output_stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white"
)
print(figure)
dev.off()

ragg::agg_tiff(
  paste0(output_stem, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white",
  compression = "lzw"
)
print(figure)
dev.off()

qa <- data.frame(
  metric = c(
    "figure_width_mm", "figure_height_mm", "coverage_plot_window_bp",
    "coverage_bin_bp", "depth_rows_27kb", "depth_rows_17kb",
    "spanning_reads_27kb", "spanning_reads_17kb",
    "zero_coverage_bases_27kb", "zero_coverage_bases_17kb"
  ),
  value = c(
    width_mm, height_mm, 24001, 100,
    nrow(depth_list[[1]]), nrow(depth_list[[2]]),
    nrow(read_list[[1]]), nrow(read_list[[2]]),
    sum(depth_list[[1]]$depth == 0), sum(depth_list[[2]]$depth == 0)
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa,
  file.path(qa_dir, "figure_data_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Created: ", output_stem, ".{svg,pdf,png,tiff}")
