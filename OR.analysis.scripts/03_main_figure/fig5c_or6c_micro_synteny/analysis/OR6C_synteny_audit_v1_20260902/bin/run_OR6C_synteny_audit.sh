#!/usr/bin/env bash
set -Eeuo pipefail

# Ppa vs Ppa_or OR6C locus audit. Coordinates are 1-based inclusive.
PPA_FA="${PPA_FA:?Set PPA_FA to the P. p. japonensis genome FASTA}"
PPA_OR_FA="${PPA_OR_FA:?Set PPA_OR_FA to the P. p. orientalis genome FASTA}"
OUT_ROOT="${OUT_ROOT:?Set OUT_ROOT to a writable run directory}"
FLANK_BP="${FLANK_BP:-500000}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools || true)}"
MINIMAP2="${MINIMAP2:-$(command -v minimap2 || true)}"

PPA_CONTIG="ChrB4"
PPA_CORE_START=81450618
PPA_CORE_END=82626544
PPA_OR_CONTIG="HiC_scaffold_8"
PPA_OR_CORE_START=80447609
PPA_OR_CORE_END=81587262

mkdir -p "$OUT_ROOT"/{alignments,logs,qc,regions,state,tables}
exec > >(tee "$OUT_ROOT/logs/run_OR6C_synteny_audit.log") 2>&1

fail() {
    printf 'status\tfail\nmessage\t%s\n' "$1" > "$OUT_ROOT/state/FAILED.tsv"
    printf 'FATAL: %s\n' "$1" >&2
    exit 1
}
trap 'rc=$?; if (( rc != 0 )) && [[ ! -s "$OUT_ROOT/state/FAILED.tsv" ]]; then printf "status\tfail\nline\t%s\nexit_code\t%s\n" "$LINENO" "$rc" > "$OUT_ROOT/state/FAILED.tsv"; fi' EXIT

[[ -n "$SAMTOOLS" && -x "$SAMTOOLS" ]] || fail "samtools was not found"
[[ -n "$MINIMAP2" && -x "$MINIMAP2" ]] || fail "minimap2 was not found"
[[ -s "$PPA_FA" ]] || fail "missing Ppa FASTA: $PPA_FA"
[[ -s "$PPA_OR_FA" ]] || fail "missing Ppa_or FASTA: $PPA_OR_FA"
[[ "$FLANK_BP" =~ ^[0-9]+$ ]] || fail "FLANK_BP must be a non-negative integer"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || fail "THREADS must be a positive integer"

for fa in "$PPA_FA" "$PPA_OR_FA"; do
    if [[ ! -s "${fa}.fai" ]]; then
        "$SAMTOOLS" faidx "$fa"
    fi
    [[ -s "${fa}.fai" ]] || fail "FASTA index was not created: ${fa}.fai"
done

extract_region() {
    local label="$1" fa="$2" contig="$3" core_start="$4" core_end="$5"
    local contig_len start end region out
    contig_len=$(awk -F '\t' -v c="$contig" '$1==c {print $2; exit}' "${fa}.fai")
    [[ "$contig_len" =~ ^[1-9][0-9]*$ ]] || fail "contig $contig was not found in ${fa}.fai"
    (( core_start >= 1 && core_end >= core_start && core_end <= contig_len )) || \
        fail "core interval $contig:$core_start-$core_end is outside contig length $contig_len"
    start=$(( core_start - FLANK_BP ))
    (( start < 1 )) && start=1
    end=$(( core_end + FLANK_BP ))
    (( end > contig_len )) && end=$contig_len
    region="${contig}:${start}-${end}"
    out="$OUT_ROOT/regions/${label}.OR6C_region.fa"
    "$SAMTOOLS" faidx "$fa" "$region" |
        awk -v h=">${label}|${region}|1-based-inclusive" 'NR==1 {$0=h} {print}' > "$out"
    [[ $(grep -c '^>' "$out") -eq 1 ]] || fail "region extraction did not yield one sequence: $label"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$fa" "$contig" "$contig_len" "$core_start" "$core_end" "$start" "$end" \
        >> "$OUT_ROOT/tables/extracted_intervals.tsv"
}

printf 'analysis_unit\tgenome_fasta\tcontig\tcontig_length\tcore_start_1\tcore_end_1\textracted_start_1\textracted_end_1\n' \
    > "$OUT_ROOT/tables/extracted_intervals.tsv"
extract_region Ppa "$PPA_FA" "$PPA_CONTIG" "$PPA_CORE_START" "$PPA_CORE_END"
extract_region Ppa_or "$PPA_OR_FA" "$PPA_OR_CONTIG" "$PPA_OR_CORE_START" "$PPA_OR_CORE_END"

printf 'analysis_unit\tsequence_bp\tN_bp\tN_fraction\n' > "$OUT_ROOT/qc/region_sequence_qc.tsv"
for label in Ppa Ppa_or; do
    awk -v label="$label" '
        BEGIN {bp=0; n=0}
        !/^>/ {
            line=toupper($0); bp+=length(line); n+=gsub(/N/, "", line)
        }
        END {
            if (bp==0) exit 2
            printf "%s\t%d\t%d\t%.8f\n", label, bp, n, n/bp
        }
    ' "$OUT_ROOT/regions/${label}.OR6C_region.fa" >> "$OUT_ROOT/qc/region_sequence_qc.tsv" || \
        fail "empty extracted sequence: $label"
done

PAF="$OUT_ROOT/alignments/Ppa_vs_Ppa_or.OR6C.asm5.paf"
FILTERED="$OUT_ROOT/alignments/Ppa_vs_Ppa_or.OR6C.asm5.identity90.min1kb.paf"
"$MINIMAP2" -x asm5 --eqx -t "$THREADS" \
    "$OUT_ROOT/regions/Ppa_or.OR6C_region.fa" \
    "$OUT_ROOT/regions/Ppa.OR6C_region.fa" > "$PAF" \
    2> "$OUT_ROOT/logs/minimap2.log"
[[ -s "$PAF" ]] || fail "minimap2 returned an empty PAF"

awk -F '\t' 'BEGIN{OFS="\t"} $11>0 && $10>=1000 && ($10/$11)>=0.90 {print}' "$PAF" > "$FILTERED"
[[ -s "$FILTERED" ]] || fail "no alignment passed >=1 kb and >=90% identity"

printf 'metric\tvalue\n' > "$OUT_ROOT/qc/alignment.summary.tsv"
awk -F '\t' '
    BEGIN {n=0; matches=0; block=0; plus=0; minus=0}
    {
        n++; matches+=$10; block+=$11
        if ($5=="+") plus++; else if ($5=="-") minus++
    }
    END {
        print "filtered_alignments\t" n
        print "summed_matching_bases\t" matches
        print "summed_alignment_block_bases\t" block
        print "weighted_identity\t" (block ? matches/block : 0)
        print "plus_strand_alignments\t" plus
        print "minus_strand_alignments\t" minus
    }
' "$FILTERED" >> "$OUT_ROOT/qc/alignment.summary.tsv"

{
    printf 'status\tpass\n'
    printf 'coordinate_system\t1-based-inclusive\n'
    printf 'flank_bp\t%s\n' "$FLANK_BP"
    printf 'threads\t%s\n' "$THREADS"
    printf 'samtools\t%s\n' "$($SAMTOOLS --version | head -1)"
    printf 'minimap2\t%s\n' "$($MINIMAP2 --version 2>&1 | head -1)"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$OUT_ROOT/state/COMPLETE.tsv"

rm -f "$OUT_ROOT/state/FAILED.tsv"
printf 'OR6C SYNTENY EXTRACTION AND ALIGNMENT COMPLETED\noutput=%s\n' "$OUT_ROOT"
