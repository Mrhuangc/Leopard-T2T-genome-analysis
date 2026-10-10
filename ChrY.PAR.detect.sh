#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YTOOLS="${REPO_ROOT}/ChrY.analysis.scripts"

minimap2 -x asm5 -c --eqx -N 200 PPA.ChrX.fa PPA.ChrY.fa > PPA.XY.paf
python3 "${YTOOLS}/PAR.scripts/call_par.py" --paf PPA.XY.paf --out PPA.PAR.candidates.tsv \
  --min-aln 10000 --min-identity 0.95 --end-window 10000000 --merge-gap 100000

python3 "${YTOOLS}/PAR.scripts/validate_reviewed_par.py" \
  --par par_region/PPA.ChrY.PAR.final.bed \
  --cap par_region/PPA.ChrY.satellite_cap.final.bed --length 44117358
