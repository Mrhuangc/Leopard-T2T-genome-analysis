#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CDR="${REPO_ROOT}/CDR.scripts"

#Generate methyl.bed and take a preliminary look at the methylation status.
bash "${CDR}/run_PPA_HiFi_CDR.primary_only.v2.sh" -c PPA.centromere.region.bed
bash "${CDR}/run_sumatran_tiger_ONT_CDR.v4.sh" PTI.genome.fna PTI.centromere.region.bed

#Based on the XXX_CDR.bed and XXX_methyl.bed generated in the previous step, merge and filter the CDRs to detect reliable CDRs.
bash "${CDR}/run_carnivora_primary_CDR_pipeline.v2.sh" single \
    PPA \
    PPA_CDR.bed \
    PPA_methyl.bed \
    PPA.centromere.region.bed \
    PPA.final.genome.fna_repeats1.SAT113.gff \
    PPA.primary_CDR_results

#plot
Rscript "${CDR}/plot_species_CDR_overview.v2.R" \
  --species PPA \
  --methyl PPA_methyl.bed \
  --chrom-summary PPA.primary_CDR_results/01.primary_selection/PPA.chromosome_summary.tsv \
  --graded PPA.primary_CDR_results/01.primary_selection/PPA.CDR.graded.tsv \
  --sat113 PPA.final.genome.fna_repeats1.SAT113.gff \
  --sat483 PPA.final.genome.fna_repeats1.SAT483.gff \
  --out-prefix PPA.CDR_overview.v2 \
  --bin 5000 \
  --smooth-bins 3 \
  --show-secondary FALSE

Rscript "${CDR}/plot_primary_CDR_mainfig.v3.R" \
  --species PPA \
  --methyl PPA_methyl.bed \
  --chrom-summary PPA.primary_CDR_results/01.primary_selection/PPA.chromosome_summary.tsv \
  --graded PPA.primary_CDR_results/01.primary_selection/PPA.CDR.graded.tsv \
  --sat113 PPA.final.genome.fna_repeats1.SAT113.gff \
  --sat483 PPA.final.genome.fna_repeats1.SAT483.gff \
  --chroms ChrB2,ChrB3,ChrX \
  --window centromere \
  --bin 5000 \
  --smooth-bins 1 \
  --minimal TRUE \
  --per-chrom TRUE \
  --out-prefix PPA.primary_CDR_main.v2
