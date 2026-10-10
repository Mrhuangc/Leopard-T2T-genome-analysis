#!/usr/bin/env Rscript

# SAT113_CDR_enrichment.v1.R
# Chromosome-stratified randomization test for association between final CDRs
# and SAT113 within curated centromere intervals.
#
# Null model: each observed CDR is relocated uniformly within the SAME curated
# centromere, preserving its length and chromosome. SAT113 coordinates are fixed.
# The supplied SAT113 file is treated as SAT113-only; its GFF Name attribute is
# not used as a filter (important for PPA TRASH output where Name is an array ID).

options(stringsAsFactors = FALSE, scipen = 999)

usage <- function(status=0L) {
  cat(paste0(
    "Usage:\n",
    "  Rscript SAT113_CDR_enrichment.v1.R \\\n",
    "    --cdr PPA.primary_CDR.bed \\\n",
    "    --sat113 PPA.SAT113.gff \\\n",
    "    --centromere PPA.centromere.region.bed \\\n",
    "    --prefix PPA --outdir PPA.SAT113_test [--nperm 10000] [--seed 20260818]\n\n",
    "Options:\n",
    "  --cdr FILE            Final Primary CDR BED (>=3 columns).\n",
    "  --sat113 FILE         SAT113-only BED or GFF/GFF3. Entire file is SAT113.\n",
    "  --centromere FILE     Curated centromere BED, one interval per chromosome.\n",
    "  --sat-format FORMAT   auto, bed, or gff [auto].\n",
    "  --prefix STR          Output prefix.\n",
    "  --outdir DIR          Output directory.\n",
    "  --nperm INT           Number of within-centromere permutations [10000].\n",
    "  --seed INT            RNG seed [20260818].\n",
    "  --no-plots            Skip PDF/PNG null-distribution plots.\n"
  ))
  quit(save="no", status=status)
}

parse_args <- function(x) {
  out <- list(`sat-format`="auto", nperm=10000L, seed=20260818L, `no-plots`=FALSE)
  i <- 1L
  value_opts <- c("cdr","sat113","centromere","sat-format","prefix","outdir","nperm","seed")
  while (i <= length(x)) {
    a <- x[[i]]
    if (a %in% c("-h","--help")) usage(0L)
    if (a == "--no-plots") { out[["no-plots"]] <- TRUE; i <- i + 1L; next }
    if (!startsWith(a, "--")) stop("Unexpected positional argument: ", a, call.=FALSE)
    key <- substring(a, 3L)
    if (!(key %in% value_opts)) stop("Unknown option: ", a, call.=FALSE)
    if (i == length(x)) stop("Missing value for ", a, call.=FALSE)
    out[[key]] <- x[[i+1L]]
    i <- i + 2L
  }
  req <- c("cdr","sat113","centromere","prefix","outdir")
  miss <- req[!vapply(req, function(z) !is.null(out[[z]]) && nzchar(out[[z]]), logical(1))]
  if (length(miss)) stop("Missing required option(s): --", paste(miss, collapse=", --"), call.=FALSE)
  out$nperm <- as.integer(out$nperm)
  out$seed <- as.integer(out$seed)
  if (!is.finite(out$nperm) || out$nperm < 100L) stop("--nperm must be >=100", call.=FALSE)
  if (!(out[["sat-format"]] %in% c("auto","bed","gff"))) stop("--sat-format must be auto, bed, or gff", call.=FALSE)
  out
}

args <- parse_args(commandArgs(trailingOnly=TRUE))
for (f in c(args$cdr, args$sat113, args$centromere)) if (!file.exists(f)) stop("Input file not found: ", f, call.=FALSE)
dir.create(args$outdir, recursive=TRUE, showWarnings=FALSE)

read_tab <- function(path) {
  read.delim(path, header=FALSE, sep="\t", quote="", comment.char="#", fill=TRUE,
             check.names=FALSE, stringsAsFactors=FALSE)
}

read_bed3 <- function(path, label="BED") {
  x <- read_tab(path)
  if (ncol(x) < 3L) stop(label, " requires >=3 columns: ", path, call.=FALSE)
  y <- data.frame(chrom=as.character(x[[1]]), start=suppressWarnings(as.numeric(x[[2]])),
                  end=suppressWarnings(as.numeric(x[[3]])), stringsAsFactors=FALSE)
  y <- y[is.finite(y$start) & is.finite(y$end) & y$start >= 0 & y$end > y$start, , drop=FALSE]
  if (!nrow(y)) stop("No valid intervals in ", path, call.=FALSE)
  y$start <- as.integer(y$start); y$end <- as.integer(y$end)
  y
}

read_cdr <- function(path) {
  x <- read_tab(path)
  if (ncol(x) < 3L) stop("CDR BED requires >=3 columns", call.=FALSE)
  out <- data.frame(chrom=as.character(x[[1]]), start=suppressWarnings(as.numeric(x[[2]])),
                    end=suppressWarnings(as.numeric(x[[3]])), stringsAsFactors=FALSE)
  out$id <- if (ncol(x) >= 4L) as.character(x[[4]]) else paste0(out$chrom, ".CDR", seq_len(nrow(out)))
  out$role <- if (ncol(x) >= 5L) as.character(x[[5]]) else "CDR"
  out$confidence <- if (ncol(x) >= 6L) as.character(x[[6]]) else "NA"
  out <- out[is.finite(out$start) & is.finite(out$end) & out$start >= 0 & out$end > out$start, , drop=FALSE]
  out$start <- as.integer(out$start); out$end <- as.integer(out$end)
  if (!nrow(out)) stop("No valid CDR intervals", call.=FALSE)
  out
}

infer_sat_format <- function(path, requested) {
  if (requested != "auto") return(requested)
  z <- tolower(path)
  if (grepl("\\.gff3?(\\.gz)?$", z)) return("gff")
  # Peek at the first data line: >=9 columns with numeric GFF start/end suggests GFF.
  con <- if (grepl("\\.gz$", z)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add=TRUE)
  repeat {
    line <- readLines(con, n=1L, warn=FALSE)
    if (!length(line)) break
    if (!nzchar(line) || startsWith(line, "#")) next
    f <- strsplit(line, "\t", fixed=TRUE)[[1]]
    if (length(f) >= 9L && suppressWarnings(!is.na(as.numeric(f[[4]]))) && suppressWarnings(!is.na(as.numeric(f[[5]])))) return("gff")
    break
  }
  "bed"
}

read_sat <- function(path, fmt) {
  fmt <- infer_sat_format(path, fmt)
  x <- read_tab(path)
  if (fmt == "gff") {
    if (ncol(x) < 5L) stop("GFF requires >=5 columns", call.=FALSE)
    # GFF: 1-based inclusive -> BED: 0-based half-open.
    out <- data.frame(chrom=as.character(x[[1]]), start=suppressWarnings(as.numeric(x[[4]]))-1,
                      end=suppressWarnings(as.numeric(x[[5]])), stringsAsFactors=FALSE)
  } else {
    if (ncol(x) < 3L) stop("SAT113 BED requires >=3 columns", call.=FALSE)
    out <- data.frame(chrom=as.character(x[[1]]), start=suppressWarnings(as.numeric(x[[2]])),
                      end=suppressWarnings(as.numeric(x[[3]])), stringsAsFactors=FALSE)
  }
  out <- out[is.finite(out$start) & is.finite(out$end) & out$start >= 0 & out$end > out$start, , drop=FALSE]
  out$start <- as.integer(out$start); out$end <- as.integer(out$end)
  if (!nrow(out)) stop("No valid SAT113 intervals", call.=FALSE)
  attr(out, "format") <- fmt
  out
}

merge_intervals <- function(x) {
  if (!nrow(x)) return(x)
  ans <- list(); k <- 0L
  for (chr in unique(x$chrom)) {
    z <- x[x$chrom == chr, , drop=FALSE]
    z <- z[order(z$start, z$end), , drop=FALSE]
    cs <- z$start[[1]]; ce <- z$end[[1]]
    if (nrow(z) > 1L) for (i in 2:nrow(z)) {
      if (z$start[[i]] <= ce) ce <- max(ce, z$end[[i]]) else {
        k <- k+1L; ans[[k]] <- data.frame(chrom=chr,start=cs,end=ce)
        cs <- z$start[[i]]; ce <- z$end[[i]]
      }
    }
    k <- k+1L; ans[[k]] <- data.frame(chrom=chr,start=cs,end=ce)
  }
  do.call(rbind, ans)
}

clip_to_centromere <- function(x, cen, label) {
  out <- list(); k <- 0L
  for (i in seq_len(nrow(x))) {
    j <- which(cen$chrom == x$chrom[[i]])
    if (!length(j)) next
    # one curated interval per chromosome is enforced below
    s <- max(x$start[[i]], cen$start[[j[[1]]]])
    e <- min(x$end[[i]], cen$end[[j[[1]]]])
    if (e > s) {
      k <- k+1L
      row <- x[i,,drop=FALSE]; row$start <- as.integer(s); row$end <- as.integer(e)
      out[[k]] <- row
    }
  }
  if (!length(out)) stop("No ", label, " intervals overlap curated centromeres", call.=FALSE)
  do.call(rbind, out)
}

ov_bp_one <- function(s, e, sat_chr) {
  if (!nrow(sat_chr)) return(0)
  sum(pmax(0, pmin(e, sat_chr$end) - pmax(s, sat_chr$start)))
}

distance_one <- function(s, e, sat_chr) {
  if (!nrow(sat_chr)) return(NA_real_)
  ov <- pmax(0, pmin(e, sat_chr$end) - pmax(s, sat_chr$start))
  if (any(ov > 0)) return(0)
  d <- ifelse(sat_chr$end <= s, s - sat_chr$end,
              ifelse(sat_chr$start >= e, sat_chr$start - e, 0))
  min(d)
}

cen <- read_bed3(args$centromere, "centromere BED")
if (anyDuplicated(cen$chrom)) stop("Centromere BED must contain exactly one interval per chromosome", call.=FALSE)
cdr <- read_cdr(args$cdr)
sat <- read_sat(args$sat113, args[["sat-format"]])
sat_format_used <- attr(sat, "format")

cdr <- clip_to_centromere(cdr, cen, "CDR")
sat <- clip_to_centromere(sat, cen, "SAT113")
sat <- merge_intervals(sat)

# Per-CDR observed statistics.
obs_rows <- vector("list", nrow(cdr))
for (i in seq_len(nrow(cdr))) {
  sc <- sat[sat$chrom == cdr$chrom[[i]],,drop=FALSE]
  ov <- ov_bp_one(cdr$start[[i]], cdr$end[[i]], sc)
  dist <- distance_one(cdr$start[[i]], cdr$end[[i]], sc)
  len <- cdr$end[[i]] - cdr$start[[i]]
  obs_rows[[i]] <- data.frame(
    chrom=cdr$chrom[[i]], start=cdr$start[[i]], end=cdr$end[[i]], id=cdr$id[[i]],
    role=cdr$role[[i]], confidence=cdr$confidence[[i]], length_bp=len,
    SAT113_overlap_bp=ov, SAT113_overlap_fraction=ov/len,
    SAT113_distance_bp=dist, SAT113_present_on_chrom=nrow(sc)>0,
    stringsAsFactors=FALSE
  )
}
obs <- do.call(rbind, obs_rows)
write.table(obs, file=file.path(args$outdir, paste0(args$prefix, ".SAT113_observed_per_CDR.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)

# Coverage of SAT113 in the actual eligible search universe: only centromeres
# carrying an input CDR are relevant because permutations preserve chromosome.
test_chroms <- unique(cdr$chrom)
cen_test <- cen[cen$chrom %in% test_chroms, , drop=FALSE]
sat_test <- sat[sat$chrom %in% test_chroms, , drop=FALSE]
cen_total_bp <- sum(cen_test$end - cen_test$start)
sat_total_bp <- sum(sat_test$end - sat_test$start)
sat_cov_fraction <- sat_total_bp / cen_total_bp

calc_metrics <- function(pos_start, pos_end) {
  ovs <- numeric(length(pos_start)); ds <- rep(NA_real_, length(pos_start))
  for (i in seq_along(pos_start)) {
    sc <- sat[sat$chrom == cdr$chrom[[i]],,drop=FALSE]
    ovs[[i]] <- ov_bp_one(pos_start[[i]], pos_end[[i]], sc)
    ds[[i]] <- distance_one(pos_start[[i]], pos_end[[i]], sc)
  }
  lens <- pos_end - pos_start
  frac_overlap <- mean(ovs > 0)
  bp_frac <- sum(ovs) / sum(lens)
  good_d <- ds[is.finite(ds)]
  c(n_overlap=sum(ovs > 0), fraction_CDR_overlapping=frac_overlap,
    total_overlap_bp=sum(ovs), CDR_bp_overlap_fraction=bp_frac,
    mean_distance_bp=if(length(good_d)) mean(good_d) else NA_real_,
    median_distance_bp=if(length(good_d)) median(good_d) else NA_real_,
    bp_fold_enrichment=if(sat_cov_fraction > 0) bp_frac/sat_cov_fraction else NA_real_,
    n_distance_evaluable=length(good_d))
}

observed_metrics <- calc_metrics(cdr$start, cdr$end)

# Precompute permissible start ranges.
cen_idx <- match(cdr$chrom, cen$chrom)
if (anyNA(cen_idx)) stop("Internal error: CDR chromosome absent from centromere BED", call.=FALSE)
cdr_len <- cdr$end - cdr$start
max_start <- cen$end[cen_idx] - cdr_len
min_start <- cen$start[cen_idx]
if (any(max_start < min_start)) stop("A CDR is longer than its curated centromere interval", call.=FALSE)

set.seed(args$seed)
metric_names <- names(observed_metrics)
null <- matrix(NA_real_, nrow=args$nperm, ncol=length(metric_names), dimnames=list(NULL, metric_names))

# Per-CDR null data needed for chromosome-specific empirical tests.
null_overlap_bp <- matrix(0, nrow=args$nperm, ncol=nrow(cdr))
null_distance <- matrix(NA_real_, nrow=args$nperm, ncol=nrow(cdr))

message("Running ", format(args$nperm, big.mark=","), " chromosome-stratified permutations...")
for (b in seq_len(args$nperm)) {
  starts <- integer(nrow(cdr))
  for (i in seq_len(nrow(cdr))) {
    if (max_start[[i]] == min_start[[i]]) starts[[i]] <- min_start[[i]] else
      starts[[i]] <- sample.int(max_start[[i]] - min_start[[i]] + 1L, 1L) - 1L + min_start[[i]]
  }
  ends <- starts + cdr_len
  null[b,] <- calc_metrics(starts, ends)
  for (i in seq_len(nrow(cdr))) {
    sc <- sat[sat$chrom == cdr$chrom[[i]],,drop=FALSE]
    null_overlap_bp[b,i] <- ov_bp_one(starts[[i]], ends[[i]], sc)
    null_distance[b,i] <- distance_one(starts[[i]], ends[[i]], sc)
  }
}

emp_upper <- function(nullv, obs) {
  z <- nullv[is.finite(nullv)]
  if (!is.finite(obs) || !length(z)) return(NA_real_)
  (1 + sum(z >= obs)) / (length(z) + 1)
}
emp_lower <- function(nullv, obs) {
  z <- nullv[is.finite(nullv)]
  if (!is.finite(obs) || !length(z)) return(NA_real_)
  (1 + sum(z <= obs)) / (length(z) + 1)
}

# Global tests.
global <- data.frame(
  metric=metric_names,
  observed=as.numeric(observed_metrics),
  null_mean=vapply(seq_along(metric_names), function(j) mean(null[,j], na.rm=TRUE), numeric(1)),
  null_median=vapply(seq_along(metric_names), function(j) median(null[,j], na.rm=TRUE), numeric(1)),
  null_q025=vapply(seq_along(metric_names), function(j) unname(quantile(null[,j], .025, na.rm=TRUE)), numeric(1)),
  null_q975=vapply(seq_along(metric_names), function(j) unname(quantile(null[,j], .975, na.rm=TRUE)), numeric(1)),
  alternative=c("greater","greater","greater","greater","less","less","greater","two.sided")
)
global$empirical_p <- NA_real_
for (j in seq_len(nrow(global))) {
  if (global$alternative[[j]] == "greater") global$empirical_p[[j]] <- emp_upper(null[,j], global$observed[[j]])
  else if (global$alternative[[j]] == "less") global$empirical_p[[j]] <- emp_lower(null[,j], global$observed[[j]])
}
global$centromere_total_bp <- cen_total_bp
global$SAT113_total_bp <- sat_total_bp
global$centromere_SAT113_coverage_fraction <- sat_cov_fraction
global$n_permutations <- args$nperm
global$seed <- args$seed
write.table(global, file=file.path(args$outdir, paste0(args$prefix, ".SAT113_global_test.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)

# Per-CDR/chromosome tests.
per <- obs
per$null_overlap_bp_mean <- colMeans(null_overlap_bp, na.rm=TRUE)
per$overlap_bp_fold_vs_null <- ifelse(per$null_overlap_bp_mean > 0, per$SAT113_overlap_bp/per$null_overlap_bp_mean, NA_real_)
per$p_overlap_bp_greater <- vapply(seq_len(nrow(per)), function(i) emp_upper(null_overlap_bp[,i], per$SAT113_overlap_bp[[i]]), numeric(1))
per$null_distance_median_bp <- vapply(seq_len(nrow(per)), function(i) median(null_distance[,i], na.rm=TRUE), numeric(1))
per$p_distance_less <- vapply(seq_len(nrow(per)), function(i) emp_lower(null_distance[,i], per$SAT113_distance_bp[[i]]), numeric(1))
per$p_overlap_bp_BH <- p.adjust(per$p_overlap_bp_greater, method="BH")
per$p_distance_BH <- p.adjust(per$p_distance_less, method="BH")
write.table(per, file=file.path(args$outdir, paste0(args$prefix, ".SAT113_per_CDR_test.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)

# Null summary and optionally useful complete null matrix (gzipped).
null_df <- as.data.frame(null)
write.table(data.frame(metric=metric_names,
                       mean=vapply(null_df, mean, numeric(1), na.rm=TRUE),
                       median=vapply(null_df, median, numeric(1), na.rm=TRUE),
                       sd=vapply(null_df, sd, numeric(1), na.rm=TRUE),
                       q025=vapply(null_df, quantile, numeric(1), probs=.025, na.rm=TRUE),
                       q975=vapply(null_df, quantile, numeric(1), probs=.975, na.rm=TRUE)),
            file=file.path(args$outdir, paste0(args$prefix, ".SAT113_permutation_summary.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)
gz <- gzfile(file.path(args$outdir, paste0(args$prefix, ".SAT113_global_null.tsv.gz")), "wt")
write.table(null_df, file=gz, sep="\t", quote=FALSE, row.names=FALSE)
close(gz)

plot_null <- function(v, obsval, main, xlab) {
  v <- v[is.finite(v)]
  if (!length(v)) { plot.new(); title(main=paste0(main, " (not evaluable)")); return() }
  if (length(unique(v)) < 2L) {
    plot(1, v[[1]], pch=16, xlim=c(.5,1.5), xaxt="n", xlab="", ylab=xlab, main=main)
    if (is.finite(obsval)) abline(h=obsval, lty=2, lwd=2)
  } else {
    hist(v, breaks="FD", main=main, xlab=xlab, border="white")
    if (is.finite(obsval)) abline(v=obsval, lty=2, lwd=2)
  }
}

if (!args[["no-plots"]]) {
  draw <- function() {
    old <- par(no.readonly=TRUE); on.exit(par(old))
    par(mfrow=c(2,2), mar=c(4.2,4.2,2.8,1.0))
    plot_null(null[,"CDR_bp_overlap_fraction"], observed_metrics[["CDR_bp_overlap_fraction"]],
              "CDR bp overlapping SAT113", "Overlap fraction")
    plot_null(null[,"fraction_CDR_overlapping"], observed_metrics[["fraction_CDR_overlapping"]],
              "CDRs overlapping SAT113", "Fraction of CDRs")
    plot_null(null[,"mean_distance_bp"], observed_metrics[["mean_distance_bp"]],
              "Mean CDR-SAT113 distance", "Distance (bp)")
    plot_null(null[,"bp_fold_enrichment"], observed_metrics[["bp_fold_enrichment"]],
              "SAT113 bp enrichment", "Observed/background fold")
    mtext(paste0(args$prefix, ": within-centromere SAT113 randomization (n=", args$nperm, ")"),
          outer=TRUE, line=-1.2, cex=1.0)
  }
  pdf(file.path(args$outdir, paste0(args$prefix, ".SAT113_randomization.pdf")), width=8.5, height=7.2)
  draw(); dev.off()
  png(file.path(args$outdir, paste0(args$prefix, ".SAT113_randomization.png")), width=1800, height=1500, res=200)
  draw(); dev.off()
}

methods <- paste0(
  "SAT113-CDR spatial association test\n\n",
  "The supplied SAT113 annotation was treated as a SAT113-only interval set; GFF feature names were not used as an additional filter. GFF coordinates were converted from 1-based inclusive to 0-based half-open coordinates, and overlapping SAT113 intervals were merged within each chromosome. Final CDRs and SAT113 intervals were restricted to the curated centromere intervals.\n\n",
  "Association was evaluated using a chromosome-stratified randomization test. In each of ", args$nperm,
  " permutations, every CDR was relocated uniformly within the same curated centromere while preserving its length and chromosome, with the observed SAT113 annotation held fixed. Global statistics included the fraction of CDRs overlapping SAT113, the fraction of total CDR bp overlapping SAT113, total overlap bp, mean and median nearest distance to SAT113, and bp-level fold enrichment relative to the fraction of the curated centromeric search space occupied by SAT113. Empirical one-sided P values were calculated as (1 + number of null statistics at least as extreme as observed)/(1 + number of permutations). Per-CDR overlap and distance tests were also calculated and adjusted across CDRs using the Benjamini-Hochberg procedure. SAT113 was not used to define or rank Primary CDRs.\n"
)
writeLines(methods, con=file.path(args$outdir, paste0(args$prefix, ".SAT113_TEST_METHODS.txt")))

message("Done: ", normalizePath(args$outdir, mustWork=FALSE))
message("SAT113 format used: ", sat_format_used)
