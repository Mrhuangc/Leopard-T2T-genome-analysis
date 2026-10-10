#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YTOOLS="${REPO_ROOT}/ChrY.analysis.scripts"

#leopard AMP07 region
MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/analyze_chry_ampliconic_units.0907.py" \
  --fasta PPA.AMP07.flank50k.new.fa \
  --repeat-gff PPA.AMP07.repeat.gff \
  --coding-gff3 PPA.ChrY.protein_coding.high_confidence.gff3 \
  --pseudogene-gff3 PPA.ChrY.reliable_pseudogenes.gff3 \
  --stainedglass-bed PPA_ChrY.5000.10000.bed.gz \
  --family RND-5_FAMILY-1989 \
  --target-region ChrY:11927310-12914358 \
  --prefix PPA_AMP07 \
  --outdir results/PPA_AMP07_new \
  --plot-all-units \
  --threads 16

#leopard CPLX01 region
MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/analyze_chry_ampliconic_units.0907.py" \
  --fasta PPA.ChrY.repeat.rich.fa \
  --repeat-gff PPA.repeat.ChrY.gff \
  --coding-gff3 PPA.ChrY.protein_coding.high_confidence.gff3 \
  --pseudogene-gff3 PPA.ChrY.reliable_pseudogenes.gff3 \
  --stainedglass-bed PPA_ChrY.5000.10000.bed.gz \
  --target-region ChrY:12914359-44117358 \
  --prefix PPA_CPLX01 \
  --outdir results_940000_all \
  --period-bp 940000 \
  --plot-all-units \
  --threads 16

#domestic cat AMP03C region
MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/analyze_chry_ampliconic_units.0907.py" \
  --fasta FCA.AMP03C.u.fa \
  --repeat-gff FCA.repeat.ChrY.gff \
  --coding-gff3 FCA.ChrY.protein_coding.high_confidence.gff3 \
  --pseudogene-gff3 FCA.ChrY.reliable_pseudogenes.gff3 \
  --stainedglass-bed FCA_ChrY.5000.10000.bed.gz \
  --target-region ChrY:12800001-19293000 \
  --prefix FCA.AMP03C \
  --plot-all-units \
  --threads 16 \
  --outdir results/FCA.AMP03C.allunit

#domestic cat CPLX01 region
MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/analyze_chry_ampliconic_units.0907.py" \
  --fasta FCA.ChrY.repeat.rich.fa \
  --repeat-gff FCA.repeat.ChrY.gff \
  --coding-gff3 FCA.ChrY.protein_coding.high_confidence.gff3 \
  --pseudogene-gff3 FCA.ChrY.reliable_pseudogenes.gff3 \
  --stainedglass-bed FCA_ChrY.5000.10000.bed.gz \
  --target-region ChrY:19293002-45999067 \
  --prefix FCA_CPLX01 \
  --outdir results.1430000.new \
  --period-bp 1430000 \
  --plot-all-units \
  --threads 16

#units identity
MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/calculate_unit_identity.py" \
  --fasta results/PPA_AMP07_new/PPA_AMP07.units.oriented.fa \
  --prefix PPA_AMP07 \
  --outdir identity_results \
  --threads 32

MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/calculate_unit_identity.py" \
  --fasta results_940000_all/PPA_CPLX01.units.oriented.fa \
  --prefix PPA_CPLX01 \
  --outdir identity_results \
  --threads 32

MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/calculate_unit_identity.py" \
  --fasta results/FCA.AMP03C.allunit/FCA.AMP03C.units.oriented.fa \
  --prefix FCA_AMP03C \
  --outdir identity_results \
  --threads 32

MPLBACKEND=Agg python3 "${YTOOLS}/repeat.units.scripts/calculate_unit_identity.py" \
  --fasta results.1430000.new/FCA_CPLX01.units.oriented.fa \
  --prefix FCA_CPLX01 \
  --outdir identity_results \
  --threads 32
