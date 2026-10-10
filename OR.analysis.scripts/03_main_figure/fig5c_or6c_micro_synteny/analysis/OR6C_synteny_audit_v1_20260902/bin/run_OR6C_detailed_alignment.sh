#!/usr/bin/env bash
set -Eeuo pipefail

# Add breakpoint-resolved PAF output to a completed OR6C synteny audit.
OUT_ROOT="${OUT_ROOT:?Set OUT_ROOT to the completed OR6C audit directory}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
MINIMAP2="${MINIMAP2:-$(command -v minimap2 || true)}"

QUERY="$OUT_ROOT/regions/Ppa.OR6C_region.fa"
TARGET="$OUT_ROOT/regions/Ppa_or.OR6C_region.fa"
PAF="$OUT_ROOT/alignments/Ppa_vs_Ppa_or.OR6C.asm5.detailed.paf"
LOG="$OUT_ROOT/logs/minimap2.detailed.log"
STATE="$OUT_ROOT/state/DETAILED_ALIGNMENT_COMPLETE.tsv"

mkdir -p "$OUT_ROOT"/{alignments,logs,state}

fail() {
    printf 'status\tfail\nmessage\t%s\n' "$1" > "$OUT_ROOT/state/DETAILED_ALIGNMENT_FAILED.tsv"
    printf 'FATAL: %s\n' "$1" >&2
    exit 1
}

[[ -n "$MINIMAP2" && -x "$MINIMAP2" ]] || fail "minimap2 was not found"
[[ -s "$QUERY" ]] || fail "missing query region: $QUERY"
[[ -s "$TARGET" ]] || fail "missing target region: $TARGET"

"$MINIMAP2" -x asm5 -c --cs=short --eqx -t "$THREADS" \
    "$TARGET" "$QUERY" > "$PAF" 2> "$LOG"

[[ -s "$PAF" ]] || fail "detailed minimap2 PAF is empty"
[[ $(awk -F '\t' '$5=="+" && $12==60 && $0~/\tcg:Z:/ && $0~/\tcs:Z:/ {n++} END{print n+0}' "$PAF") -ge 1 ]] || \
    fail "no plus-strand MAPQ60 alignment with cg and cs tags"

{
    printf 'status\tpass\n'
    printf 'method\tminimap2 asm5 with CIGAR and short cs tag\n'
    printf 'coordinate_system\tPAF 0-based half-open; FASTA headers retain source 1-based inclusive interval\n'
    printf 'paf\t%s\n' "$PAF"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$STATE"

rm -f "$OUT_ROOT/state/DETAILED_ALIGNMENT_FAILED.tsv"
printf 'OR6C DETAILED ALIGNMENT COMPLETED\noutput=%s\n' "$PAF"
