#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUX="${REPO_ROOT}/scripts"

biser -o PPA.biser -t 20 --gc-heap 128G path/PPA_carnivora/huabao_carnivora.softmasked.fasta

python "${AUX}/filter_biser_sd_sedeflike_v2.py" \
  -i PPA.biser \
  -o PPA.BISER.SD85 \
  --satellite-bed path/PPA_carnivora/huabao_carnivora.filteredRepeats.Satellite.bed \
  --satellite-mode either \
  --min-aln-len 1000 \
  --min-frac-match 0.85 \
  --min-frac-match-indel 0.50 \
  --prefix PPA

python "${AUX}/biser_sd_stats_sedeflike_v2.py" \
  -i PPA.BISER.SD85.final.bedpe.tsv \
  -o PPA.BISER.SD85.stats \
  --fai path/PPA_carnivora/huabao_carnivora.softmasked.fasta.fai
