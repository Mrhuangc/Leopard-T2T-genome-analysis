args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
    stop(
        "usage: nc_tree_parser.R TREE TYPE1_IDS TYPE2_IDS PYTHON_IDS OUTPUT_PREFIX"
    )
}

suppressPackageStartupMessages(library(ape))

tree_file <- args[[1L]]
type1_file <- args[[2L]]
type2_file <- args[[3L]]
python_ids_file <- args[[4L]]
output_prefix <- args[[5L]]

tree <- read.tree(tree_file)
tree <- makeNodeLabel(tree, method = "number", prefix = "Node")
type1 <- scan(type1_file, what = "character", quiet = TRUE)
type2 <- scan(type2_file, what = "character", quiet = TRUE)
known <- c(type1, type2)

if (anyDuplicated(known)) {
    stop("known Type I and Type II OR identifiers are not unique")
}
missing_known <- setdiff(known, tree$tip.label)
if (length(missing_known) > 0L) {
    stop(sprintf("%d known OR identifiers are absent from the tree", length(missing_known)))
}

# This is the operative logic of the published Tree_parser_2.R. ape::getMRCA
# replaces phytools::findMRCA(..., type='node') without changing the node used
# by ape::root and ape::extract.clade.
mrca_or <- getMRCA(tree, known)
if (is.null(mrca_or)) {
    stop("could not identify the initial known-OR MRCA")
}
rooted <- root(tree, node = mrca_or, resolve.root = TRUE)
mrca_or_rooted <- getMRCA(rooted, known)
if (is.null(mrca_or_rooted)) {
    stop("could not identify the known-OR MRCA after rerooting")
}
or_genes <- extract.clade(rooted, node = mrca_or_rooted)$tip.label
nc_r_ids <- setdiff(or_genes, known)

python_ids <- scan(python_ids_file, what = "character", quiet = TRUE)
r_only <- sort(setdiff(nc_r_ids, python_ids))
python_only <- sort(setdiff(python_ids, nc_r_ids))
parity_status <- if (
    length(r_only) == 0L && length(python_only) == 0L
) "pass" else "disagree"

writeLines(nc_r_ids, sprintf("%s.nc_r.ids.txt", output_prefix))
writeLines(r_only, sprintf("%s.nc_r_only.ids.txt", output_prefix))
writeLines(python_only, sprintf("%s.python_only.ids.txt", output_prefix))
write.tree(rooted, file = sprintf("%s.nc_r.rooted.treefile", output_prefix))

summary_lines <- c(
    "{",
    sprintf('  "status": "%s",', parity_status),
    '  "method": "ape root(node=MRCA, resolve.root=TRUE)",',
    sprintf('  "r_version": "%s",', paste(R.version$major, R.version$minor, sep = ".")),
    sprintf('  "ape_version": "%s",', as.character(packageVersion("ape"))),
    sprintf('  "known_type1_or": %d,', length(type1)),
    sprintf('  "known_type2_or": %d,', length(type2)),
    sprintf('  "nc_r_retained_candidates": %d,', length(nc_r_ids)),
    sprintf('  "python_retained_candidates": %d,', length(python_ids)),
    sprintf('  "nc_r_only_candidates": %d,', length(r_only)),
    sprintf('  "python_only_candidates": %d', length(python_only)),
    "}"
)
writeLines(summary_lines, sprintf("%s.nc_r_parity.summary.json", output_prefix))
