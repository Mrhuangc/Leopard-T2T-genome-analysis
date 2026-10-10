#!/usr/bin/env Rscript

# Rebuild panel b at the same 183 mm width as panel a, using the original
# seven-genome macro-synteny summary table copied into this output directory.

script_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
cli_args <- commandArgs(trailingOnly = TRUE)
if (!length(cli_args) %in% c(0L, 2L)) stop("Usage: Rscript SCRIPT.R [INPUT_DIR OUTPUT_DIR]")
input_dir <- normalizePath(if (length(cli_args)) cli_args[[1]] else dirname(script_arg), winslash = "/", mustWork = TRUE)
out_dir <- if (length(cli_args)) cli_args[[2]] else dirname(script_arg)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

metrics <- read.delim(file.path(input_dir, "Ppa_vs_felids.macro_synteny.summary.tsv"),
                      sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
comparison_units <- c("Fca", "Lge", "Pbe", "Ple", "Ppa_or", "Pti")
species <- c("F. catus", "L. geoffroyi", "P. bengalensis", "P. leo",
             "P. p. orientalis", "P. tigris")
metrics <- metrics[match(comparison_units, metrics$comparison_unit), ]
if (anyNA(metrics$comparison_unit)) stop("One or more comparison units are missing from the source summary")
coverage <- pmin(metrics$query_coverage, metrics$target_coverage) * 100
identity <- metrics$weighted_identity * 100

source_values <- data.frame(
  comparison_unit = comparison_units,
  species = species,
  query_coverage_percent = metrics$query_coverage * 100,
  target_coverage_percent = metrics$target_coverage * 100,
  minimum_reciprocal_regional_coverage_percent = coverage,
  weighted_nucleotide_identity_percent = identity,
  formula = "minimum(query_coverage, target_coverage); weighted_identity",
  source = "Ppa_vs_felids.macro_synteny.summary.tsv"
)
write.table(source_values, file.path(out_dir, "panel_b_source_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

identity_col <- colorRampPalette(c("#BDD1E3", "#163C60"))(101)
point_col <- function(x) {
  identity_col[pmax(1L, pmin(101L, round((x - 94) / 6 * 100) + 1L))]
}

draw_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    xscale = c(0, 1), yscale = c(0, 1),
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc")
  ))

  # Header shares the typographic hierarchy of panel a.
  grid::grid.text("b", x = 0.010, y = 0.985, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold"))
  grid::grid.text("Regional sequence homology to P. p. japonensis",
                  x = 0.045, y = 0.985, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold"))
  grid::grid.text("Points show minimum reciprocal coverage; labels report weighted nucleotide identity.",
                  x = 0.045, y = 0.935, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 7.5, col = "#4D4D4D"))

  x0 <- 0.245
  x1 <- 0.755
  x_at <- function(v) x0 + (v - 70) / 30 * (x1 - x0)
  ys <- seq(0.785, 0.285, length.out = 6)

  # Species labels and coverage tracks.
  for (i in seq_along(species)) {
    grid::grid.text(species[i], x = 0.222, y = ys[i], just = "right",
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 8.2, fontface = "italic"))
    grid::grid.lines(x = grid::unit(c(x_at(70), x_at(coverage[i])), "npc"),
                     y = grid::unit(c(ys[i], ys[i]), "npc"),
                     gp = grid::gpar(col = "#BDBDBD", lwd = 1.2))
    grid::grid.circle(x = x_at(coverage[i]), y = ys[i], r = grid::unit(2.0, "mm"),
                      gp = grid::gpar(fill = point_col(identity[i]), col = "#303840", lwd = 0.7))
    grid::grid.text(sprintf("%.1f%% identity", identity[i]),
                    x = x_at(coverage[i]) + 0.016, y = ys[i], just = "left",
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 7.6, col = "#1A1A1A"))
  }

  # Shared x-axis for minimum reciprocal coverage.
  grid::grid.lines(x = grid::unit(c(x_at(70), x_at(100)), "npc"),
                   y = grid::unit(c(0.185, 0.185), "npc"),
                   gp = grid::gpar(col = "#202020", lwd = 0.8))
  ticks <- seq(70, 100, 10)
  for (v in ticks) {
    grid::grid.lines(x = grid::unit(c(x_at(v), x_at(v)), "npc"),
                     y = grid::unit(c(0.185, 0.173), "npc"),
                     gp = grid::gpar(col = "#202020", lwd = 0.7))
    grid::grid.text(paste0(v, "%"), x = x_at(v), y = 0.145,
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 7.4))
  }
  grid::grid.text("Minimum reciprocal regional coverage",
                  x = (x0 + x1) / 2, y = 0.075,
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 8.2))

  # Identity color scale, matching the blue sequential scale in panel a.
  cb_x <- 0.875
  cb_y0 <- 0.285
  cb_y1 <- 0.785
  for (j in seq_len(100)) {
    ya <- cb_y0 + (j - 1) / 100 * (cb_y1 - cb_y0)
    yb <- cb_y0 + j / 100 * (cb_y1 - cb_y0)
    grid::grid.rect(x = cb_x, y = (ya + yb) / 2,
                    width = grid::unit(5.0, "mm"), height = grid::unit(yb - ya, "npc"),
                    gp = grid::gpar(fill = identity_col[j], col = NA))
  }
  for (v in 94:100) {
    yy <- cb_y0 + (v - 94) / 6 * (cb_y1 - cb_y0)
    grid::grid.lines(x = grid::unit(c(cb_x - 0.014, cb_x - 0.008), "npc"),
                     y = grid::unit(c(yy, yy), "npc"),
                     gp = grid::gpar(col = "white", lwd = 0.7))
    grid::grid.text(paste0(v, "%"), x = cb_x + 0.022, y = yy, just = "left",
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 7.2))
  }
  grid::grid.text("Nucleotide\nidentity", x = cb_x, y = 0.835,
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 7.5), just = "center")
  grid::popViewport()
}

png(file.path(out_dir, "panel_b_regional_homology.png"),
    width = 2161, height = 1087, res = 300, bg = "white", type = "cairo")
draw_panel()
dev.off()

cairo_pdf(file.path(out_dir, "panel_b_regional_homology.pdf"),
          width = 183 / 25.4, height = 92 / 25.4, family = "Arial")
draw_panel()
dev.off()

cat("Source summary: Ppa_vs_felids.macro_synteny.summary.tsv\n")
cat("Output width: 183 mm (2161 px at 300 dpi); output height: 92 mm\n")
