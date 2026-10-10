#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUX="${REPO_ROOT}/scripts"

meryl count k=15 output my_merylDB PPA.final.genome.fna
meryl print greater-than distinct=0.9998 my_merylDB > repetitive_k15.txt
winnowmap -t 30 -k 15 -W repetitive_k15.txt -x asm5 PPA.final.genome.fna Ppar.fa > PPA.Ppar.winnowmap.paf
python "${AUX}/paf_to_PUR_bed.modified.py" -i PPA.Ppar.winnowmap.paf -o PPA.Ppar.PUR.bed --min_mapq 20 --min_len 1000 -s PPA.genome.sizes --primary_only --long_primary_min_len 1000000 --kept_paf PPA.Ppar.PUR.used_alignments.tsv
