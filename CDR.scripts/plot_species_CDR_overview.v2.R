#!/usr/bin/env Rscript

# plot_species_CDR_overview.v2.R
# Portrait supplementary overview inspired by the giant-panda centromere CDR layout.
# Each chromosome uses only two visual layers:
#   (1) satellite strip (SAT113/SAT483)
#   (2) enlarged methylation profile with Primary CDR indicated by a black bar
# Raw fragmented CDR-Finder calls are NOT plotted.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c('-h','--help'))) {
cat("Usage:
  Rscript plot_species_CDR_overview.v2.R \\
    --species PTI \\
    --methyl PTI_methyl.bed \\
    --chrom-summary PTI.chromosome_summary.tsv \\
    --graded PTI.CDR.graded.tsv \\
    --out-prefix PTI.CDR_overview.v2 \\
    [--sat113 PTI.SAT113.gff] [--sat483 PTI.SAT483.gff] \\
    [--bin 5000] [--smooth-bins 3] [--pad 0] \\
    [--show-secondary FALSE] [--show-uncertain TRUE] \\
    [--width 7.2] [--row-height 0.62]\n\n",
"Design notes:\n",
"  * Final Primary CDRs only; raw CDR-Finder fragments are not shown.\n",
"  * Methylation is coverage-weighted within bins.\n",
"  * --smooth-bins is display-only rolling median smoothing (default 3 bins).\n",
"  * Primary CDR = black bar under methylation; uncertain = open/hatched-like bar.\n",
"  * PDF and PNG are written explicitly to <out-prefix>.* (no Rplots.pdf).\n", sep='')
quit(save='no', status=0)
}

parse_args <- function(x) {
  z <- list(); i <- 1
  while (i <= length(x)) {
    k <- sub('^--','',x[i])
    if (i == length(x) || grepl('^--', x[i+1])) {z[[k]] <- TRUE; i <- i+1}
    else {z[[k]] <- x[i+1]; i <- i+2}
  }
  z
}
o <- parse_args(args)
req <- c('species','methyl','chrom-summary','graded','out-prefix')
miss <- req[!req %in% names(o)]
if (length(miss)) stop('Missing required arguments: ', paste(miss, collapse=', '))

species <- o[['species']]
methyl_file <- o[['methyl']]
summary_file <- o[['chrom-summary']]
graded_file <- o[['graded']]
out_prefix <- o[['out-prefix']]
sat113_file <- if ('sat113' %in% names(o)) o[['sat113']] else NA_character_
sat483_file <- if ('sat483' %in% names(o)) o[['sat483']] else NA_character_
bin_size <- if ('bin' %in% names(o)) as.integer(o[['bin']]) else 5000L
smooth_bins <- if ('smooth-bins' %in% names(o)) as.integer(o[['smooth-bins']]) else 3L
pad_bp <- if ('pad' %in% names(o)) as.integer(o[['pad']]) else 0L
show_secondary <- if ('show-secondary' %in% names(o)) toupper(o[['show-secondary']]) == 'TRUE' else FALSE
show_uncertain <- if ('show-uncertain' %in% names(o)) toupper(o[['show-uncertain']]) == 'TRUE' else TRUE
fig_width <- if ('width' %in% names(o)) as.numeric(o[['width']]) else 7.2
row_height <- if ('row-height' %in% names(o)) as.numeric(o[['row-height']]) else 0.62

read_tsv <- function(f) read.table(f, sep='\t', header=TRUE, quote='', comment.char='', stringsAsFactors=FALSE, check.names=FALSE)

read_sat <- function(path, label) {
  if (is.na(path) || !file.exists(path)) return(data.frame(chrom=character(), start=integer(), end=integer(), family=character()))
  x <- read.table(path, sep='\t', header=FALSE, quote='', comment.char='#', stringsAsFactors=FALSE)
  if (ncol(x) >= 9) {
    d <- data.frame(chrom=as.character(x[[1]]), start=as.integer(x[[4]])-1L, end=as.integer(x[[5]]), family=label)
  } else if (ncol(x) >= 3) {
    d <- data.frame(chrom=as.character(x[[1]]), start=as.integer(x[[2]]), end=as.integer(x[[3]]), family=label)
  } else stop('Cannot parse satellite file: ', path)
  d[d$end > d$start,,drop=FALSE]
}

natural_chr_key <- function(x) {
  # Works for ChrA1...ChrF2, ChrX, ChrY and chr1...chr22.
  y <- sub('^(Chr|chr)','',x)
  y[y=='X'] <- 'ZZX'
  y[y=='Y'] <- 'ZZY'
  # Pad trailing numeric part so A2 precedes A10.
  m <- regexec('^([^0-9]*)([0-9]+)$', y)
  rr <- regmatches(y,m)
  out <- y
  for (i in seq_along(rr)) if (length(rr[[i]])==3) out[i] <- paste0(rr[[i]][2], sprintf('%04d', as.integer(rr[[i]][3])))
  out
}

weighted_bins <- function(df, st, en, bin) {
  nbin <- max(1L, ceiling((en-st)/bin))
  bs <- st + (0:(nbin-1))*bin
  be <- pmin(bs+bin,en)
  out <- data.frame(start=bs,end=be,mid=(bs+be)/2,methyl=NA_real_,cov=0)
  if (!nrow(df)) return(out)
  idx <- floor((df$pos-st)/bin)+1L
  keep <- idx>=1L & idx<=nbin & is.finite(df$cov) & is.finite(df$meth)
  if (!any(keep)) return(out)
  df <- df[keep,,drop=FALSE]; idx <- idx[keep]
  cs <- tapply(df$cov, idx, sum)
  mw <- tapply(df$meth*df$cov, idx, sum)
  ii <- as.integer(names(cs))
  out$cov[ii] <- as.numeric(cs)
  out$methyl[ii] <- as.numeric(mw/cs)
  out
}

roll_median <- function(x, k=3L) {
  if (k <= 1L || length(x) < 2L) return(x)
  if (k %% 2L == 0L) k <- k + 1L
  h <- floor(k/2L)
  y <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    lo <- max(1L, i-h); hi <- min(length(x), i+h)
    v <- x[lo:hi]
    if (any(is.finite(v))) y[i] <- median(v, na.rm=TRUE)
  }
  y
}

# Inputs
sm <- read_tsv(summary_file)
gr <- read_tsv(graded_file)
sm$plot_start <- pmax(0, as.integer(sm$centromere_start)-pad_bp)
sm$plot_end <- as.integer(sm$centromere_end)+pad_bp
sm <- sm[order(natural_chr_key(sm$chrom)),,drop=FALSE]

sat113 <- read_sat(sat113_file,'SAT113')
sat483 <- read_sat(sat483_file,'SAT483')

m <- read.table(methyl_file, sep='\t', header=FALSE, quote='', comment.char='', stringsAsFactors=FALSE)
if (ncol(m) < 12) stop('Methyl BED requires at least 12 columns: ', methyl_file)
colnames(m)[1:12] <- c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb','cov','meth','meth_reads')
m$chrom <- as.character(m$chrom); m$pos <- as.integer(m$start); m$cov <- as.numeric(m$cov); m$meth <- as.numeric(m$meth)

# Pre-bin by chromosome.
binned <- vector('list',nrow(sm)); names(binned) <- sm$chrom
for (i in seq_len(nrow(sm))) {
  chr <- sm$chrom[i]; st <- sm$plot_start[i]; en <- sm$plot_end[i]
  z <- m[m$chrom==chr & m$pos>=st & m$pos<en,c('pos','cov','meth'),drop=FALSE]
  b <- weighted_bins(z,st,en,bin_size)
  b$display_methyl <- roll_median(b$methyl,smooth_bins)
  binned[[chr]] <- b
}

write.table(sm[,c('chrom','plot_start','plot_end','primary_id','primary_role','primary_confidence')], paste0(out_prefix,'.plot_ranges.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

# Colors chosen to remain close to the existing centromere composition figures.
COL_SAT113 <- '#E76F9A'
COL_SAT483 <- '#36A9A6'
COL_PRIMARY <- '#111111'
COL_SECONDARY <- '#7A7A7A'
COL_BG <- '#BDBDBD'
COL_LINE <- '#111111'
COL_UNCERTAIN <- '#FFFFFF'

plot_page <- function(filename, type=c('pdf','png')) {
  type <- match.arg(type)
  n <- nrow(sm)
  fig_height <- max(7.5, 1.25 + n*row_height)
  if (type=='pdf') pdf(filename, width=fig_width, height=fig_height, useDingbats=FALSE)
  else png(filename, width=fig_width, height=fig_height, units='in', res=300)
  on.exit(dev.off(), add=TRUE)

  par(mar=c(3.0,5.7,2.0,2.0), xpd=NA, lend='butt')
  total_h <- n*1.0 + 1.1
  plot(0,0,type='n',xlim=c(0,1),ylim=c(0,total_h),axes=FALSE,xlab='',ylab='')
  title(main=paste0(species,' centromere methylation and CDRs'), cex.main=0.95, font.main=2)

  # compact legend at top-right
  yleg <- total_h-0.25
  x0 <- 0.58
  rect(x0,yleg-0.045,x0+0.025,yleg+0.045,col=COL_SAT483,border=NA); text(x0+0.032,yleg,'SAT483',adj=c(0,0.5),cex=0.62)
  rect(x0+0.15,yleg-0.045,x0+0.175,yleg+0.045,col=COL_SAT113,border=NA); text(x0+0.182,yleg,'SAT113',adj=c(0,0.5),cex=0.62)
  segments(x0+0.30,yleg,x0+0.35,yleg,lwd=3.2,col=COL_PRIMARY); text(x0+0.36,yleg,'Primary CDR',adj=c(0,0.5),cex=0.62)

  for (ii in seq_len(n)) {
    i <- ii
    ybase <- (n-i)*1.0 + 0.45
    ysat <- ybase + 0.50
    ymet0 <- ybase - 0.26
    ymet1 <- ybase + 0.34
    ycdr <- ymet0 - 0.075

    chr <- sm$chrom[i]; st <- sm$plot_start[i]; en <- sm$plot_end[i]
    span <- en-st
    tx <- function(x) (x-st)/span

    # labels
    text(-0.035,ysat,chr,adj=c(1,0.5),cex=0.72)
    text(-0.035,ymet0+0.30,paste0(round(span/1e6,2),' Mb'),adj=c(1,0.5),cex=0.46,col='grey40')

    # satellite strip, visually similar to panda supplementary layout
    rect(0,ysat-0.07,1,ysat+0.07,col=COL_BG,border=NA)
    s <- sat483[sat483$chrom==chr & sat483$end>st & sat483$start<en,,drop=FALSE]
    if (nrow(s)) for (k in seq_len(nrow(s))) rect(tx(max(st,s$start[k])),ysat-0.07,tx(min(en,s$end[k])),ysat+0.07,col=COL_SAT483,border=NA)
    s <- sat113[sat113$chrom==chr & sat113$end>st & sat113$start<en,,drop=FALSE]
    if (nrow(s)) for (k in seq_len(nrow(s))) rect(tx(max(st,s$start[k])),ysat-0.07,tx(min(en,s$end[k])),ysat+0.07,col=COL_SAT113,border=NA)

    # methylation panel: taller than v1, with subtle 50% guide only.
    segments(0,ymet0,1,ymet0,col='grey70',lwd=0.45)
    segments(0,ymet0+0.5*(ymet1-ymet0),1,ymet0+0.5*(ymet1-ymet0),col='grey88',lwd=0.45)
    segments(0,ymet1,1,ymet1,col='grey85',lwd=0.35)

    b <- binned[[chr]]
    ok <- is.finite(b$display_methyl)
    if (sum(ok)>=2L) {
      xx <- tx(b$mid[ok]); yy <- ymet0 + (b$display_methyl[ok]/100)*(ymet1-ymet0)
      lines(xx,yy,col=COL_LINE,lwd=0.72)
    }

    # final CDR marker directly beneath methylation; no separate CDR track.
    ss <- gr[gr$chrom==chr,,drop=FALSE]
    pri <- ss[ss$role=='PRIMARY',,drop=FALSE]
    unc <- ss[ss$role=='PRIMARY_UNCERTAIN',,drop=FALSE]
    sec <- if (show_secondary) ss[ss$role=='SECONDARY' & ss$confidence %in% c('HIGH','MODERATE'),,drop=FALSE] else ss[FALSE,,drop=FALSE]

    if (nrow(sec)) for (k in seq_len(nrow(sec))) segments(tx(sec$start[k]),ycdr,tx(sec$end[k]),ycdr,col=COL_SECONDARY,lwd=2.2)
    if (nrow(pri)) for (k in seq_len(nrow(pri))) segments(tx(pri$start[k]),ycdr,tx(pri$end[k]),ycdr,col=COL_PRIMARY,lwd=3.4)
    if (show_uncertain && nrow(unc)) for (k in seq_len(nrow(unc))) {
      rect(tx(unc$start[k]),ycdr-0.018,tx(unc$end[k]),ycdr+0.018,col=COL_UNCERTAIN,border='#D95F0E',lwd=0.7)
    }

    # very compact y labels, like the panda figure.
    text(-0.006,ymet0,'0',adj=c(1,0.5),cex=0.38,col='grey45')
    text(-0.006,ymet1,'100',adj=c(1,0.5),cex=0.38,col='grey45')

    conf <- sm$primary_confidence[i]
    role <- sm$primary_role[i]
    lab <- ifelse(is.na(role) || role=='' || role=='NO_PRIMARY','No primary',paste(role,conf))
    text(1.012,ymet0+0.30*(ymet1-ymet0),lab,adj=c(0,0.5),cex=0.46,
         col=ifelse(conf=='HIGH','grey20',ifelse(conf=='MODERATE','#D95F0E','grey48')))
  }

  mtext(paste0('Relative position within each plotted centromere interval; methylation shown as ',bin_size/1000,'-kb weighted bins',
               if (smooth_bins>1) paste0(' (',smooth_bins,'-bin rolling median for display)') else ''),
        side=1,line=1.9,cex=0.58,col='grey30')
}

plot_page(paste0(out_prefix,'.pdf'),'pdf')
plot_page(paste0(out_prefix,'.png'),'png')
