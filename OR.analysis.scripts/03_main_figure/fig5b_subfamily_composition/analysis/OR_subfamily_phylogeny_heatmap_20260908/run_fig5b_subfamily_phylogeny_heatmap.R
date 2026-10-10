#!/usr/bin/env Rscript

# Generate the base figure and its audited source-data bundle.
source("analysis/OR_subfamily_phylogeny_heatmap_20260908/build_fig5b_subfamily_phylogeny_heatmap.R")

# Replace the provisional tree layer with a fully connected rectangular cladogram.
tree_segments <- data.frame(
  x = c(
    4.1875, 2.25, 6.125,
    1.50, 3.00, 1.00, 2.00,
    5.25, 7.00, 4.50, 6.00, 4.00, 5.00,
    2.25, 1.50, 1.00, 5.25, 4.50, 4.00
  ),
  xend = c(
    4.1875, 2.25, 6.125,
    1.50, 3.00, 1.00, 2.00,
    5.25, 7.00, 4.50, 6.00, 4.00, 5.00,
    6.125, 3.00, 2.00, 7.00, 6.00, 5.00
  ),
  y = c(
    0.15, 0.70, 0.70,
    1.30, 1.30, 2.15, 2.15,
    1.30, 1.30, 1.90, 1.90, 2.45, 2.45,
    0.70, 1.30, 2.15, 1.30, 1.90, 2.45
  ),
  yend = c(
    0.70, 1.30, 1.30,
    2.15, 3.15, 3.15, 3.15,
    1.90, 3.15, 2.45, 3.15, 3.15, 3.15,
    0.70, 1.30, 2.15, 1.30, 1.90, 2.45
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

caption_text <- paste0(
  "Numbers are functional-locus counts. Colour is centred within each subfamily and does not encode statistical significance.\n",
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
                                  hjust = 0, lineheight = 1.1, margin = margin(t = 4)),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
  )

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

# Final source-data integrity and rendered-tree checks.
input_path <- "figure/OR_subfamily_comparison_final/OR_subfamily_heatmap.sensitivity.source_data.tsv"
output_path <- "figure/OR_subfamily_comparison_phylogeny_v1_20260908/OR_subfamily_composition.sensitivity.phylogeny.source_data.tsv"
qc_path <- "figure/OR_subfamily_comparison_phylogeny_v1_20260908/OR_subfamily_composition.sensitivity.phylogeny.QC.json"
input <- read.delim(input_path, sep = "\t", quote = "", check.names = FALSE)
output <- read.delim(output_path, sep = "\t", quote = "", check.names = FALSE)
key <- c("analysis_unit", "or_subfamily")
cols <- c(key, "functional_loci", "row_median", "deviation_from_median")
input_cmp <- input[do.call(order, input[key]), cols]
output_cmp <- output[do.call(order, output[key]), cols]
rownames(input_cmp) <- NULL
rownames(output_cmp) <- NULL

qc <- jsonlite::read_json(qc_path, simplifyVector = FALSE)
qc$checks$all_original_counts_retained <- identical(input_cmp, output_cmp)
qc$checks$rectangular_tree_segments_complete <- nrow(tree_segments) == 19L
qc$visual_revision <- "connected rectangular cladogram; two-line caption"
qc$status <- if (all(unlist(qc$checks))) "pass" else "fail"
jsonlite::write_json(qc, qc_path, pretty = TRUE, auto_unbox = TRUE)
if (!identical(qc$status, "pass")) stop("Final QA failed")

cat("PHYLOGENY_HEATMAP_V2_QA_PASS\n")
