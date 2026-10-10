#!/usr/bin/env Rscript

# Rebuild panel a from the archived, filtered primary alignment blocks.
# The source table was derived from the six primary minimap2 PAF files using
# block length >= 10 kb, identity >= 90%, and MAPQ >= 20.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)][1])
cli_args <- commandArgs(trailingOnly = TRUE)
if (!length(cli_args) %in% c(0L, 2L)) stop("Usage: Rscript SCRIPT.R [INPUT_DIR OUTPUT_DIR]")
input_dir <- normalizePath(if (length(cli_args)) cli_args[[1]] else dirname(script_arg), winslash = "/", mustWork = TRUE)
out_dir <- if (length(cli_args)) cli_args[[2]] else dirname(script_arg)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

comparison_units <- c("Fca", "Lge", "Pbe", "Ple", "Ppa_or", "Pti")
species <- c("F. catus", "L. geoffroyi", "P. bengalensis", "P. leo",
             "P. p. orientalis", "P. tigris")
blocks_path <- file.path(input_dir, "panel_a_filtered_primary_blocks.tsv")
if (!file.exists(blocks_path)) stop("Missing filtered primary-block table: ", blocks_path)
blocks <- read.delim(blocks_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("unit", "target_start", "target_end", "query_oriented_start",
                      "query_oriented_end", "identity", "block_length", "mapq")
if (!all(required_columns %in% names(blocks))) stop("Filtered block table is missing required columns")
if (nrow(blocks) != 3800L || !setequal(unique(blocks$unit), comparison_units)) {
  stop("Filtered block table does not match the six-comparison, 3,800-block snapshot")
}
if (any(blocks$block_length < 10000 | blocks$identity < 0.90 | blocks$mapq < 20)) {
  stop("Filtered block table contains an alignment below the declared thresholds")
}

blue <- colorRampPalette(c("#BFD3E6", "#123B5D"))(101)
block_color <- blue[pmax(1L, pmin(101L, round((blocks$identity - 0.90) / 0.10 * 100) + 1L))]

draw_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 1), yscale = c(0, 1)))

  grid::grid.text("a", x = 0.010, y = 0.985, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold"))
  grid::grid.text("Macrosynteny across the two largest felid OR arrays",
                  x = 0.045, y = 0.985, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold"))
  grid::grid.text("Primary alignment blocks are shown relative to the P. p. japonensis region.",
                  x = 0.045, y = 0.935, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 7.5, col = "#4D4D4D"))

  col_x <- c(0.095, 0.345, 0.595)
  plot_w <- 0.22
  row_y <- c(0.515, 0.15)
  plot_h <- 0.29
  breaks <- c(0, 15, 30, 45)

  for (i in seq_along(comparison_units)) {
    row <- if (i <= 3L) 1L else 2L
    col <- ((i - 1L) %% 3L) + 1L
    x_left <- col_x[col]
    y_bottom <- row_y[row]
    d <- blocks[blocks$unit == comparison_units[i], , drop = FALSE]
    ix <- match(rownames(d), rownames(blocks))

    grid::grid.rect(x = x_left + plot_w / 2, y = y_bottom + plot_h + 0.027,
                    width = grid::unit(plot_w, "npc"), height = grid::unit(0.028, "npc"),
                    gp = grid::gpar(fill = NA, col = "#222222", lwd = 0.7))
    grid::grid.text(species[i], x = x_left + plot_w / 2, y = y_bottom + plot_h + 0.027,
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 8.2, fontface = "italic"))

    grid::pushViewport(grid::viewport(
      x = x_left + plot_w / 2, y = y_bottom + plot_h / 2,
      width = grid::unit(plot_w, "npc"), height = grid::unit(plot_h, "npc"),
      xscale = c(0, 47), yscale = c(0, 47), clip = "on"
    ))
    grid::grid.segments(x0 = 0, y0 = 0, x1 = 47, y1 = 47,
                        default.units = "native",
                        gp = grid::gpar(col = "#BDBDBD", lwd = 0.45, lty = 2))
    grid::grid.segments(
      x0 = d$target_start / 1e6, y0 = d$query_oriented_start / 1e6,
      x1 = d$target_end / 1e6, y1 = d$query_oriented_end / 1e6,
      default.units = "native",
      gp = grid::gpar(col = block_color[ix], lwd = 0.42, alpha = 0.68, lineend = "round")
    )
    grid::grid.lines(x = c(0, 47), y = c(0, 0), default.units = "native",
                     gp = grid::gpar(col = "#222222", lwd = 0.7))
    grid::grid.lines(x = c(0, 0), y = c(0, 47), default.units = "native",
                     gp = grid::gpar(col = "#222222", lwd = 0.7))
    for (v in breaks) {
      grid::grid.lines(x = c(v, v), y = c(0, -0.65), default.units = "native",
                       gp = grid::gpar(col = "#222222", lwd = 0.55))
      grid::grid.lines(x = c(-0.65, 0), y = c(v, v), default.units = "native",
                       gp = grid::gpar(col = "#222222", lwd = 0.55))
    }
    grid::popViewport()

    if (col == 1L) {
      for (v in breaks) {
        yy <- y_bottom + (v / 47) * plot_h
        grid::grid.text(as.character(v), x = x_left - 0.012, y = yy,
                        just = "right", gp = grid::gpar(fontfamily = "Arial", fontsize = 7.4))
      }
    }
    if (row == 2L) {
      for (v in breaks) {
        xx <- x_left + (v / 47) * plot_w
        grid::grid.text(as.character(v), x = xx, y = y_bottom - 0.020,
                        gp = grid::gpar(fontfamily = "Arial", fontsize = 7.4))
      }
    }
  }

  grid::grid.text("Position in the P. p. japonensis region (Mb)",
                  x = 0.455, y = 0.055,
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 8.2))
  grid::grid.text("Position in the comparison region (Mb)",
                  x = 0.025, y = 0.50, rot = 90,
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 8.2))

  # Identity colour scale (90-100%), matching panel b's blue palette.
  cb_x <- 0.895
  cb_y0 <- 0.235
  cb_y1 <- 0.765
  for (j in seq_len(100)) {
    ya <- cb_y0 + (j - 1) / 100 * (cb_y1 - cb_y0)
    yb <- cb_y0 + j / 100 * (cb_y1 - cb_y0)
    grid::grid.rect(x = cb_x, y = (ya + yb) / 2,
                    width = grid::unit(5.0, "mm"), height = grid::unit(yb - ya, "npc"),
                    gp = grid::gpar(fill = blue[j], col = NA))
  }
  for (v in c(90, 92, 95, 98, 100)) {
    yy <- cb_y0 + (v - 90) / 10 * (cb_y1 - cb_y0)
    grid::grid.lines(x = grid::unit(c(cb_x - 0.014, cb_x - 0.008), "npc"),
                     y = grid::unit(c(yy, yy), "npc"),
                     gp = grid::gpar(col = "white", lwd = 0.7))
    grid::grid.text(paste0(v, "%"), x = cb_x + 0.022, y = yy, just = "left",
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 6.8))
  }
  grid::grid.text("Nucleotide\nidentity", x = cb_x, y = 0.815,
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 7.5), just = "center")
  grid::popViewport()
}

png(file.path(out_dir, "panel_a_macro_alignment.png"),
    width = 2161, height = 1511, res = 300, bg = "white", type = "cairo")
draw_panel()
dev.off()

cairo_pdf(file.path(out_dir, "panel_a_macro_alignment.pdf"),
          width = 183 / 25.4, height = 128 / 25.4, family = "Arial")
draw_panel()
dev.off()

cat("Filtered primary blocks:", nrow(blocks), "\n")
cat("Output width: 183 mm (2161 px at 300 dpi); output height: 128 mm\n")
