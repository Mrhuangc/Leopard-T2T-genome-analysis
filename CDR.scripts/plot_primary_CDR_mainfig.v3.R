#!/usr/bin/env Rscript

# plot_primary_CDR_mainfig.v3.R
# CDR-Finder-like BAR methylation panels for final post-processed Primary CDRs.
# No satellite track and no separate CDR track: a transparent red rectangle
# marks the final Primary CDR directly behind black methylation bars.
# Ideal for manual placement below centromere_merge_plot.

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0 || any(args %in% c('-h','--help'))) {
cat("Usage:
  Rscript plot_primary_CDR_mainfig.v3.R \\
    --species PPA \\
    --methyl PPA_methyl.bed \\
    --graded PPA.CDR.graded.tsv \\
    --chrom-summary PPA.chromosome_summary.tsv \\
    --out-prefix PPA.primary_CDR_main \\
    [--chroms ChrB2,ChrB3,ChrX] \\
    [--window centromere|primary] [--pad 250000] \\
    [--bin 5000] [--smooth-bins 1] \\
    [--show-secondary FALSE] [--per-chrom TRUE] [--minimal TRUE]\n\n",
"The red translucent rectangle represents the FINAL post-processed CDR, not raw CDR-Finder calls.\n",
"--window centromere is recommended when stacking below a centromere_merge_plot that uses the same interval.\n",
"--minimal TRUE removes most titles/legends so per-chromosome panels can be assembled externally.\n", sep='')
quit(save='no',status=0)
}

parse_args <- function(x){z<-list();i<-1;while(i<=length(x)){k<-sub('^--','',x[i]);if(i==length(x)||grepl('^--',x[i+1])){z[[k]]<-TRUE;i<-i+1}else{z[[k]]<-x[i+1];i<-i+2}};z}
o<-parse_args(args)
req<-c('species','methyl','graded','chrom-summary','out-prefix')
miss<-req[!req%in%names(o)]; if(length(miss)) stop('Missing: ',paste(miss,collapse=', '))

species<-o[['species']]; methyl_file<-o[['methyl']]; graded_file<-o[['graded']]; summary_file<-o[['chrom-summary']]; out_prefix<-o[['out-prefix']]
chroms_user<-if('chroms'%in%names(o)) strsplit(o[['chroms']],',')[[1]] else NULL
window_mode<-if('window'%in%names(o)) o[['window']] else 'centromere'
pad_bp<-if('pad'%in%names(o)) as.integer(o[['pad']]) else 250000L
bin_size<-if('bin'%in%names(o)) as.integer(o[['bin']]) else 5000L
smooth_bins<-if('smooth-bins'%in%names(o)) as.integer(o[['smooth-bins']]) else 1L
show_secondary<-if('show-secondary'%in%names(o)) toupper(o[['show-secondary']])=='TRUE' else FALSE
per_chrom<-if('per-chrom'%in%names(o)) toupper(o[['per-chrom']])=='TRUE' else TRUE
minimal<-if('minimal'%in%names(o)) toupper(o[['minimal']])=='TRUE' else TRUE

read_tsv<-function(f) read.table(f,sep='\t',header=TRUE,quote='',comment.char='',stringsAsFactors=FALSE,check.names=FALSE)
weighted_bins<-function(df,st,en,bin){nb<-max(1L,ceiling((en-st)/bin));bs<-st+(0:(nb-1))*bin;be<-pmin(bs+bin,en);out<-data.frame(start=bs,end=be,mid=(bs+be)/2,methyl=NA_real_,cov=0);if(!nrow(df))return(out);idx<-floor((df$pos-st)/bin)+1L;keep<-idx>=1&idx<=nb&is.finite(df$cov)&is.finite(df$meth);if(!any(keep))return(out);df<-df[keep,,drop=FALSE];idx<-idx[keep];cs<-tapply(df$cov,idx,sum);mw<-tapply(df$meth*df$cov,idx,sum);ii<-as.integer(names(cs));out$cov[ii]<-as.numeric(cs);out$methyl[ii]<-as.numeric(mw/cs);out}
roll_median<-function(x,k=1L){if(k<=1L)return(x);if(k%%2L==0L)k<-k+1L;h<-floor(k/2L);y<-rep(NA_real_,length(x));for(i in seq_along(x)){v<-x[max(1,i-h):min(length(x),i+h)];if(any(is.finite(v)))y[i]<-median(v,na.rm=TRUE)};y}

sm<-read_tsv(summary_file); gr<-read_tsv(graded_file)
m<-read.table(methyl_file,sep='\t',header=FALSE,quote='',comment.char='',stringsAsFactors=FALSE)
if(ncol(m)<12)stop('Methyl BED requires >=12 columns')
colnames(m)[1:12]<-c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb','cov','meth','meth_reads')
m$chrom<-as.character(m$chrom);m$pos<-as.integer(m$start);m$cov<-as.numeric(m$cov);m$meth<-as.numeric(m$meth)

if(is.null(chroms_user)) chroms<-sm$chrom[sm$primary_role%in%c('PRIMARY','PRIMARY_UNCERTAIN')] else chroms<-chroms_user
regions<-data.frame()
for(chr in chroms){s<-sm[sm$chrom==chr,,drop=FALSE];if(!nrow(s))next;if(window_mode=='primary'){g<-gr[gr$chrom==chr & gr$role%in%c('PRIMARY','PRIMARY_UNCERTAIN'),,drop=FALSE];if(nrow(g)){st<-max(0,min(g$start)-pad_bp);en<-max(g$end)+pad_bp}else{st<-s$centromere_start[1];en<-s$centromere_end[1]}}else{st<-s$centromere_start[1];en<-s$centromere_end[1]};regions<-rbind(regions,data.frame(chrom=chr,start=as.integer(st),end=as.integer(en)))}
write.table(regions,paste0(out_prefix,'.regions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

RED_FILL<-rgb(1,0,0,alpha=0.24)
RED_BORDER<-rgb(0.85,0,0,alpha=0.65)
ORANGE_FILL<-rgb(1,0.5,0,alpha=0.16)

plot_one<-function(chr,st,en,show_title=TRUE){
  span<-en-st; tx<-function(x)(x-st)/span
  par(mar=if(minimal)c(2.7,3.8,0.4,0.4) else c(3.2,4.4,2.1,0.8),xpd=NA)
  plot(0,0,type='n',xlim=c(0,1),ylim=c(0,100),axes=FALSE,xlab='',ylab='')

  # final post-processed CDR shading FIRST, so methylation line remains visible above it.
  gg<-gr[gr$chrom==chr,,drop=FALSE]
  sec<-if(show_secondary) gg[gg$role=='SECONDARY' & gg$confidence%in%c('HIGH','MODERATE'),,drop=FALSE] else gg[FALSE,,drop=FALSE]
  if(nrow(sec))for(i in seq_len(nrow(sec)))rect(tx(max(st,sec$start[i])),0,tx(min(en,sec$end[i])),100,col=ORANGE_FILL,border=NA)
  pri<-gg[gg$role=='PRIMARY',,drop=FALSE]
  if(nrow(pri))for(i in seq_len(nrow(pri)))rect(tx(max(st,pri$start[i])),0,tx(min(en,pri$end[i])),100,col=RED_FILL,border=RED_BORDER,lwd=0.5)
  unc<-gg[gg$role=='PRIMARY_UNCERTAIN',,drop=FALSE]
  if(nrow(unc))for(i in seq_len(nrow(unc)))rect(tx(max(st,unc$start[i])),0,tx(min(en,unc$end[i])),100,col=rgb(1,0.65,0,alpha=0.12),border='#D95F0E',lwd=0.7,lty=2)

  segments(0,0,1,0,col='grey35',lwd=0.5)
  segments(0,50,1,50,col='grey88',lwd=0.45)
  segments(0,100,1,100,col='grey88',lwd=0.35)

  z<-m[m$chrom==chr & m$pos>=st & m$pos<en,c('pos','cov','meth'),drop=FALSE]
  b<-weighted_bins(z,st,en,bin_size); b$methyl<-roll_median(b$methyl,smooth_bins)
  ok<-is.finite(b$methyl)
  if(any(ok)){
    ii<-which(ok)
    for(j in ii){
      x1<-tx(b$start[j]); x2<-tx(b$end[j])
      if(x2>x1 && b$methyl[j]>0) rect(x1,0,x2,b$methyl[j],col='black',border='black',lwd=0.15)
    }
  }

  axis(2,at=c(0,50,100),labels=c('0','50','100'),las=1,cex.axis=0.65,tck=-0.015)
  pretty_pos<-pretty(c(st,en),n=4);pretty_pos<-pretty_pos[pretty_pos>=st&pretty_pos<=en]
  axis(1,at=tx(pretty_pos),labels=round(pretty_pos/1e6,2),cex.axis=0.63,tck=-0.015)
  if(!minimal){mtext('Mean CpG methylation (%)',side=2,line=2.5,cex=0.78);mtext('Position (Mb)',side=1,line=2.0,cex=0.76)}
  if(show_title) title(main=paste0(species,' ',chr),cex.main=0.82)
}

# combined panel
n<-nrow(regions); nc<-if(n<=3)1 else if(n<=6)2 else 3; nr<-ceiling(n/nc)
pdf(paste0(out_prefix,'.pdf'),width=4.2*nc,height=2.05*nr,useDingbats=FALSE)
op<-par(mfrow=c(nr,nc));for(i in seq_len(n))plot_one(regions$chrom[i],regions$start[i],regions$end[i],TRUE);par(op);dev.off()
png(paste0(out_prefix,'.png'),width=4.2*nc,height=2.05*nr,units='in',res=300)
op<-par(mfrow=c(nr,nc));for(i in seq_len(n))plot_one(regions$chrom[i],regions$start[i],regions$end[i],TRUE);par(op);dev.off()

# per-chromosome, minimal track for centromere_merge_plot assembly
if(per_chrom){for(i in seq_len(n)){chr<-regions$chrom[i];pdf(paste0(out_prefix,'.',chr,'.pdf'),width=6.2,height=1.55,useDingbats=FALSE);plot_one(chr,regions$start[i],regions$end[i],FALSE);dev.off();png(paste0(out_prefix,'.',chr,'.png'),width=6.2,height=1.55,units='in',res=300,bg='transparent');plot_one(chr,regions$start[i],regions$end[i],FALSE);dev.off()}}
