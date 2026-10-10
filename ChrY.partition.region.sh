#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YTOOLS="${REPO_ROOT}/ChrY.analysis.scripts"

# Candidate sequence-class windows
PARTITION_CONFIG="${PARTITION_CONFIG:-${YTOOLS}/partition.scripts/PPA.example.yaml}"
REVIEWED="${YTOOLS}/results/04_comparative"
GENE_CONTEXT="${YTOOLS}/par_xdeg_4tracks"

samtools faidx PPA.ChrY.fa
minimap2 \
  -t 32 \
  -DP \
  -k19 \
  -w19 \
  -m200 \
  PPA.ChrY.fa \
  PPA.ChrY.fa \
  > PPA.Y.self.paf

minimap2 -x asm5 -c --eqx -N 200 PPA.ChrX.fa PPA.ChrY.fa > PPA.XY.paf
python3 "${YTOOLS}/PAR.scripts/call_par.py" \
  --paf PPA.XY.paf --out PPA.PAR.candidates.tsv \
  --min-aln 10000 --min-identity 0.95 --end-window 10000000 --merge-gap 100000

python3 "${YTOOLS}/partition.scripts/repeat_window_summary.py" \
  --config "${PARTITION_CONFIG}" \
  --fai PPA.ChrY.fa.fai \
  --out PPA.repeat_windows.tsv

python3 "${YTOOLS}/partition.scripts/classify_y_windows.py" \
  --config "${PARTITION_CONFIG}" \
  --fai PPA.ChrY.fa.fai \
  --xy-paf PPA.XY.paf \
  --self-paf PPA.Y.self.paf \
  --par PPA.PAR.candidates.tsv \
  --repeats PPA.repeat_windows.tsv \
  --out PPA.sequence_class_windows.tsv

#reclassify
python3 "${YTOOLS}/partition.scripts/reclassify_xdeg_ampliconic.py" --partitions-only
python3 "${YTOOLS}/partition.scripts/validate_partition.py" \
  --bed "${REVIEWED}/PPA.ChrY.structural_partition.revised.bed" --length 44117358
python3 "${YTOOLS}/partition.scripts/validate_partition.py" \
  --bed "${REVIEWED}/FCA.ChrY.structural_partition.revised.bed" --length 45999067

python3 "${YTOOLS}/partition.scripts/build_par_xdeg_from_annotations.py" \
  --fca-y-high "${YTOOLS}/annotation.scripts/results/FCA/FCA.ChrY.protein_coding.high_confidence.gff3" \
  --fca-y-provisional "${YTOOLS}/annotation.scripts/results/FCA/FCA.ChrY.protein_coding.provisional.gff3" \
  --fca-partition "${REVIEWED}/FCA.ChrY.structural_partition.revised.bed" \
  --ppa-y-high "${YTOOLS}/annotation.scripts/results/PPA/PPA.ChrY.protein_coding.high_confidence.gff3" \
  --ppa-y-provisional "${YTOOLS}/annotation.scripts/results/PPA/PPA.ChrY.protein_coding.provisional.gff3" \
  --ppa-partition "${REVIEWED}/PPA.ChrY.structural_partition.revised.bed" \
  --x-catalog path/ChrX.gff3 \
  --y-order-output "${GENE_CONTEXT}/PAR_XDG_gene_order.tsv" \
  --output "${GENE_CONTEXT}/base_input.tsv"

python3 "${YTOOLS}/partition.scripts/reclassify_xdeg_ampliconic.py" \
  --fca-fasta FCA.ChrY.fa \
  --ppa-fasta PPA.ChrY.fa \
  --base-input "${GENE_CONTEXT}/base_input.tsv"
