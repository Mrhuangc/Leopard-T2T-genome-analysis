#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUX="${REPO_ROOT}/scripts"

set -euo pipefail

# ============================================================
# PPA SAT113 array-specific divergence pipeline v2
# - reuses existing per-array MAFFT alignments with --resume
# - always rebuilds globally phase-normalized inter-array alignment
# - produces all-copy and full-length-only pairwise metrics
# - merges nearby TRASH2 blocks into SAT-rich regions
# - optionally runs the R visualization script
# ============================================================

INPUT="${INPUT:-PPA.final.genome.fna_repeats_with_seq.SAT113.csv}"   #produce from TRASH2
OUTDIR="${OUTDIR:-PPA.SAT113.array_pipeline}"
PREFIX="${PREFIX:-PPA.SAT113}"
THREADS="${THREADS:-20}"
MAFFT="${MAFFT:-mafft}"
PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"
RUN_PLOTS="${RUN_PLOTS:-1}"
HEATMAP_TOP_N="${HEATMAP_TOP_N:-60}"

PIPELINE="${PIPELINE:-${AUX}/satellite_array_pipeline_v2.py}"
PLOT_SCRIPT="${PLOT_SCRIPT:-${AUX}/plot_SAT113_array_results.R}"

"$PYTHON" "$PIPELINE" \
  -i "$INPUT" \
  -o "$OUTDIR" \
  --prefix "$PREFIX" \
  --min-array-copies 5 \
  --threads "$THREADS" \
  --mafft "$MAFFT" \
  --sequence-orientation as_is \
  --phase-k 9 \
  --phase-min-identity 0.35 \
  --inter-array-phase-min-identity 0.30 \
  --min-column-occupancy 0.50 \
  --min-consensus-support 0.50 \
  --full-length-min-ratio 0.80 \
  --full-length-max-ratio 1.20 \
  --full-length-min-coverage 0.80 \
  --max-pairwise-pairs 10000 \
  --position-bins 20 \
  --merge-block-gap 10000 \
  --resume

if [[ "$RUN_PLOTS" == "1" ]]; then
  "$RSCRIPT" "$PLOT_SCRIPT" \
    "$OUTDIR/07_results" \
    "$PREFIX" \
    "$OUTDIR/08_plots" \
    "$HEATMAP_TOP_N"
fi
