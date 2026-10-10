#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_ROOT="${PACKAGE_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd -P)}"
REGIONS_TSV="${REGIONS_TSV:-$PACKAGE_ROOT/input/regions.tsv}"
OUT_ROOT="${OUT_ROOT:-$PACKAGE_ROOT/results}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-32}}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools || true)}"
MINIMAP2="${MINIMAP2:-$(command -v minimap2 || true)}"
PYTHON="${PYTHON:-$(command -v python || true)}"

# The regions manifest uses relative genome paths under PACKAGE_ROOT.
cd "$PACKAGE_ROOT"

mkdir -p "$OUT_ROOT"/{regions,alignments,tables,qc,logs,state}
rm -f "$OUT_ROOT/state/FAILED.tsv" "$OUT_ROOT/state/COMPLETE.tsv"

fail() {
    printf 'status\tfail\nmessage\t%s\n' "$1" > "$OUT_ROOT/state/FAILED.tsv"
    printf 'FATAL: %s\n' "$1" >&2
    exit 1
}
trap 'fail "command failed at line $LINENO"' ERR

[[ -s "$REGIONS_TSV" ]] || fail "regions manifest is missing or empty: $REGIONS_TSV"
[[ -n "$SAMTOOLS" && -x "$SAMTOOLS" ]] || fail "samtools was not found"
[[ -n "$MINIMAP2" && -x "$MINIMAP2" ]] || fail "minimap2 was not found"
[[ -n "$PYTHON" && -x "$PYTHON" ]] || fail "python was not found"

{
    "$SAMTOOLS" --version | head -1
    "$MINIMAP2" --version
    "$PYTHON" --version
} > "$OUT_ROOT/qc/tool_versions.txt" 2>&1

printf 'analysis_unit\tscientific_name\tgenome_fasta\tcontig\trequested_start_1\trequested_end_1\textracted_start_1\textracted_end_1\textracted_length_bp\tn_fraction\n' > "$OUT_ROOT/qc/region_sequence_qc.tsv"

tail -n +2 "$REGIONS_TSV" | while IFS=$'\t' read -r unit species genome contig start end left_id right_id; do
    [[ -s "$genome" ]] || fail "genome FASTA is missing or empty for $unit: $genome"
    [[ -s "$genome.fai" ]] || "$SAMTOOLS" faidx "$genome"
    contig_len="$(awk -F '\t' -v id="$contig" '$1==id {print $2; exit}' "$genome.fai")"
    [[ "$contig_len" =~ ^[0-9]+$ ]] || fail "contig $contig is absent from $genome.fai"
    (( start < 1 )) && start=1
    (( end > contig_len )) && end="$contig_len"
    (( start <= end )) || fail "invalid clipped interval for $unit"
    fasta="$OUT_ROOT/regions/$unit.ChrD1_macro_region.fa"
    "$SAMTOOLS" faidx "$genome" "$contig:$start-$end" \
      | sed "1s/.*/>$unit|$contig:$start-$end/" > "$fasta"
    [[ -s "$fasta" ]] || fail "empty extracted FASTA for $unit"
    read -r seq_len n_count < <(awk '!/^>/{gsub(/[[:space:]]/,""); n+=length($0); s=toupper($0); gsub(/[^N]/,"",s); a+=length(s)} END{print n+0,a+0}' "$fasta")
    expected=$((end - start + 1))
    (( seq_len == expected )) || fail "length mismatch for $unit: expected $expected, observed $seq_len"
    n_fraction="$(awk -v n="$n_count" -v total="$seq_len" 'BEGIN{printf "%.8f", n/total}')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$unit" "$species" "$genome" "$contig" "$start" "$end" "$start" "$end" "$seq_len" "$n_fraction" \
      >> "$OUT_ROOT/qc/region_sequence_qc.tsv"
done

reference="$OUT_ROOT/regions/Ppa.ChrD1_macro_region.fa"
[[ -s "$reference" ]] || fail "Ppa reference region was not extracted"
for unit in Fca Lge Pbe Ple Ppa_or Pti; do
    query="$OUT_ROOT/regions/$unit.ChrD1_macro_region.fa"
    paf="$OUT_ROOT/alignments/Ppa_vs_${unit}.asm5.paf"
    "$MINIMAP2" -x asm5 -c --cs=short --eqx -t "$THREADS" \
      "$reference" "$query" > "$paf" 2> "$OUT_ROOT/logs/minimap2.${unit}.log"
    [[ -s "$paf" ]] || fail "minimap2 produced an empty PAF for Ppa vs $unit"
done

"$PYTHON" "$PACKAGE_ROOT/bin/summarize_paf.py" \
  --paf-dir "$OUT_ROOT/alignments" \
  --output-tsv "$OUT_ROOT/tables/Ppa_vs_felids.macro_synteny.summary.tsv" \
  --output-json "$OUT_ROOT/qc/OR_ChrD1_macro_synteny.summary.json"

{
    printf 'status\tpass\n'
    printf 'comparisons\t6\n'
    printf 'reference\tPpa\n'
    printf 'coordinate_system\tmanifest 1-based inclusive; PAF 0-based half-open\n'
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$OUT_ROOT/state/COMPLETE.tsv"
rm -f "$OUT_ROOT/state/FAILED.tsv"
printf 'CHR D1 MACRO-SYNTENY ANALYSIS COMPLETED\noutput=%s\n' "$OUT_ROOT"
