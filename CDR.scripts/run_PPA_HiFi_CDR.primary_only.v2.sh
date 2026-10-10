#!/usr/bin/env bash
#SBATCH --job-name=PPA_HiFi_CDR
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=PPA_HiFi_CDR.%j.out
#SBATCH --error=PPA_HiFi_CDR.%j.err

# ==============================================================================
# Panthera pardus HiFi methylation + CDR workflow
#
# Main steps:
#   1. Discover all PacBio HiFi BAMs and verify MM/ML tags.
#   2. Use a BAM FOFN as pbmm2 input to jointly align all runs to the T2T genome.
#      The output is one coordinate-sorted, merged modBAM.
#   3. Filter the merged BAM with SAMtools -F 2308, retaining mapped primary
#      alignments only (remove unmapped, secondary and supplementary records).
#   4. Verify read counts, mapping statistics, and preservation of MM/ML tags.
#   5. Prepare centromere regions with configurable flanks.
#   6. Generate a CpG-combined bedMethyl track from the filtered BAM.
#   7. Run CDR-Finder with restrict_alr=false for non-alpha-satellite centromeres.
#
# Required:
#   pbmm2, samtools, bedtools, modkit, tabix, snakemake, conda, git
#
# Usage:
#   sbatch run_PPA_HiFi_CDR.sh \
#       -c /path/to/PPA.centromere.region.bed
#
# Optional:
#   sbatch run_PPA_HiFi_CDR.sh \
#       -c /path/to/PPA.centromere.region.bed \
#       -o /path/to/PPA.HiFi.CDR \
#       -t 48 \
#       -s 500000
#
# Notes:
#   * Input HiFi BAMs must contain MM/ML or Mm/Ml tags.
#   * Input PacBio HiFi BAMs are expected to be unaligned; input integrity
#     is checked with `samtools quickcheck -u`. Aligned output BAMs are
#     checked with strict `samtools quickcheck`.
#   * The centromere BED must have one interval per chromosome/contig.
#   * Coordinates are assumed to be standard 0-based, half-open BED coordinates.
#   * Set FORCE=1 in the environment, or use -f, to rebuild existing results.
# ==============================================================================

set -Eeuo pipefail
umask 002

trap 'echo "[ERROR] line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# ----------------------------- fixed input paths -------------------------------

HIFI_DIR="${HIFI_DIR:-/path/to/carnivora_project/huabao/raw_data/PPA_hifi}"
REF="${REF:-/path/to/carnivora_project/huabao/CDR/PPA/PPA.final.genome.fna}"

# Default centromere BED. Override with: -c /path/file.bed
CENTROMERE_BED="${CENTROMERE_BED:-/path/to/carnivora_project/huabao/CDR/PPA/PPA.centromere.region.bed}"

# Default output directory. Override with: -o /path/output
OUTDIR="${OUTDIR:-/path/to/carnivora_project/huabao/CDR/PPA}"

SAMPLE="PPA"
THREADS="${SLURM_CPUS_PER_TASK:-48}"
REGION_SLOP=500000
FORCE="${FORCE:-0}"

# CDR-Finder settings.
CDR_FINDER_DIR="${CDR_FINDER_DIR:-${OUTDIR}/software/CDR-Finder}"
CDR_FINDER_GIT="${CDR_FINDER_GIT:-https://github.com/EichlerLab/CDR-Finder.git}"
AUTO_CLONE_CDR_FINDER="${AUTO_CLONE_CDR_FINDER:-1}"
REPEATMASKER_SPECIES="${REPEATMASKER_SPECIES:-Carnivora}"

WINDOW_SIZE="${WINDOW_SIZE:-5000}"
HEIGHT_THRESHOLD="${HEIGHT_THRESHOLD:-0.34}"
PROM_THRESHOLD="${PROM_THRESHOLD:-0.30}"
BASELINE_AVG_METHYL="${BASELINE_AVG_METHYL:-0.40}"
BP_MERGE="${BP_MERGE:-1}"
EXTEND_EDGES_STD="${EXTEND_EDGES_STD:--1}"

# pbmm2 sorting settings.
SORT_THREADS="${SORT_THREADS:-8}"
SORT_MEMORY_PER_THREAD="${SORT_MEMORY_PER_THREAD:-8G}"

usage() {
    cat <<'EOF'
Usage:
  sbatch run_PPA_HiFi_CDR.sh -c CENTROMERE.bed [options]

Options:
  -c FILE   Centromere BED; one interval per chromosome/contig.
  -o DIR    Output directory.
  -t INT    Threads.
  -s INT    Flanking bases added to each centromere interval [500000].
  -f        Force rebuilding existing major outputs.
  -h        Show this help.

Environment overrides:
  CDR_FINDER_DIR, AUTO_CLONE_CDR_FINDER, REPEATMASKER_SPECIES,
  WINDOW_SIZE, HEIGHT_THRESHOLD, PROM_THRESHOLD, BASELINE_AVG_METHYL,
  BP_MERGE, EXTEND_EDGES_STD, SORT_THREADS, SORT_MEMORY_PER_THREAD.
EOF
}

while getopts ":c:o:t:s:fh" opt; do
    case "${opt}" in
        c) CENTROMERE_BED="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        s) REGION_SLOP="${OPTARG}" ;;
        f) FORCE=1 ;;
        h) usage; exit 0 ;;
        :) echo "[ERROR] -${OPTARG} requires an argument." >&2; usage; exit 2 ;;
        \?) echo "[ERROR] Unknown option: -${OPTARG}" >&2; usage; exit 2 ;;
    esac
done

# Recalculate default software path if OUTDIR was changed and CDR_FINDER_DIR was
# not explicitly exported before launching the script.
if [[ -z "${CDR_FINDER_DIR_EXPLICIT:-}" && "${CDR_FINDER_DIR}" == */PPA.HiFi.CDR/software/CDR-Finder ]]; then
    CDR_FINDER_DIR="${OUTDIR}/software/CDR-Finder"
fi

# ------------------------------- helper functions ------------------------------

timestamp() {
    date '+%F %T'
}

log() {
    echo "[$(timestamp)] $*"
}

die() {
    echo "[$(timestamp)] [ERROR] $*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found in PATH: $1"
}

run_step() {
    local marker="$1"
    shift

    if [[ "${FORCE}" -eq 0 && -s "${marker}" ]]; then
        log "Skip completed step: ${marker}"
        return 0
    fi

    "$@"
    date '+%F %T' > "${marker}"
}

sample_bam_records() {
    local bam="$1"
    local out="$2"
    local n="${3:-1000}"

    # samtools normally receives SIGPIPE when head exits. Temporarily disable
    # pipefail so that this expected condition does not terminate the workflow.
    set +o pipefail
    samtools view "${bam}" | head -n "${n}" > "${out}"
    local rc_samtools=${PIPESTATUS[0]}
    set -o pipefail

    if [[ "${rc_samtools}" -ne 0 && "${rc_samtools}" -ne 141 ]]; then
        die "samtools view failed while sampling ${bam}; exit=${rc_samtools}"
    fi
}

count_mod_tags_in_sam() {
    local sam="$1"
    awk -F'\t' '
        {
            mm=0; ml=0; mn=0
            for (i=12; i<=NF; i++) {
                if ($i ~ /^(MM|Mm):Z:/) mm=1
                if ($i ~ /^(ML|Ml):B:/) ml=1
                if ($i ~ /^MN:i:/) mn=1
            }
            n++
            n_mm += mm
            n_ml += ml
            n_mn += mn
        }
        END {
            printf "%d\t%d\t%d\t%d\n", n+0, n_mm+0, n_ml+0, n_mn+0
        }
    ' "${sam}"
}

# ---------------------------------- checks -------------------------------------

for cmd in pbmm2 samtools bedtools modkit tabix snakemake conda git awk sed sort uniq realpath; do
    require_cmd "${cmd}"
done

[[ -d "${HIFI_DIR}" ]] || die "HiFi directory not found: ${HIFI_DIR}"
[[ -s "${REF}" ]] || die "Reference FASTA not found or empty: ${REF}"
[[ -s "${CENTROMERE_BED}" ]] || die \
    "Centromere BED not found: ${CENTROMERE_BED}. Supply it with -c /path/PPA.centromere.region.bed"

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer."
[[ "${REGION_SLOP}" =~ ^[0-9]+$ ]] || die "REGION_SLOP must be a non-negative integer."
[[ "${SORT_THREADS}" =~ ^[1-9][0-9]*$ ]] || die "SORT_THREADS must be a positive integer."

mkdir -p \
    "${OUTDIR}"/{00.input,01.reference,02.alignment,03.regions,04.methylation,05.cdr_finder,06.final,logs,tmp,software,conda_envs,markers}

export TMPDIR="${OUTDIR}/tmp"
export LC_ALL=C

# Convert all relevant paths to absolute paths.
HIFI_DIR="$(realpath "${HIFI_DIR}")"
REF="$(realpath "${REF}")"
CENTROMERE_BED="$(realpath "${CENTROMERE_BED}")"
OUTDIR="$(realpath "${OUTDIR}")"
CDR_FINDER_DIR="${CDR_FINDER_DIR/#\~/${HOME}}"

FOFN="${OUTDIR}/00.input/${SAMPLE}.hifi.bam.fofn"
INPUT_TAG_QC="${OUTDIR}/00.input/${SAMPLE}.input.MM_ML.QC.tsv"
INPUT_READ_COUNTS="${OUTDIR}/00.input/${SAMPLE}.input.read_counts.tsv"

REF_MMI="${OUTDIR}/01.reference/${SAMPLE}.hifi.mmi"

MERGED_BAM="${OUTDIR}/02.alignment/${SAMPLE}.hifi.all_runs.merged.aligned.bam"
PRIMARY_BAM="${OUTDIR}/02.alignment/${SAMPLE}.hifi.all_runs.primary.mapped.bam"
ALIGN_FLAGSTAT="${OUTDIR}/02.alignment/${SAMPLE}.hifi.all_runs.flagstat.txt"
ALIGN_IDXSTAT="${OUTDIR}/02.alignment/${SAMPLE}.hifi.all_runs.idxstats.tsv"
ALIGN_FILTER_STATS="${OUTDIR}/02.alignment/${SAMPLE}.alignment_filter.stats.tsv"
ALIGN_TAG_QC="${OUTDIR}/02.alignment/${SAMPLE}.primary.mapped.MM_ML.QC.tsv"

CENTROMERE_3COL="${OUTDIR}/03.regions/${SAMPLE}.centromere.input.3col.bed"
CENTROMERE_SLOP_BED="${OUTDIR}/03.regions/${SAMPLE}.centromere.plus${REGION_SLOP}.bed"

BEDMETHYL="${OUTDIR}/04.methylation/${SAMPLE}.centromere.CpG.combined.bed.gz"
MODKIT_SUMMARY="${OUTDIR}/04.methylation/${SAMPLE}.centromere.modkit.summary.tsv"

CDR_CONFIG="${OUTDIR}/05.cdr_finder/config.${SAMPLE}.yaml"
CDR_RESULTS="${OUTDIR}/05.cdr_finder/results"
CDR_LOGS="${OUTDIR}/05.cdr_finder/logs"
CDR_BENCHMARKS="${OUTDIR}/05.cdr_finder/benchmarks"
CDR_BED="${CDR_RESULTS}/bed/${SAMPLE}_CDR.bed"

log "Output directory: ${OUTDIR}"
log "Reference: ${REF}"
log "Centromere BED: ${CENTROMERE_BED}"
log "Threads: ${THREADS}"

# Record tool versions.
{
    echo -e "tool\tversion"
    echo -e "pbmm2\t$(pbmm2 --version 2>&1 | head -n1)"
    echo -e "samtools\t$(samtools --version | head -n1)"
    echo -e "bedtools\t$(bedtools --version | head -n1)"
    echo -e "modkit\t$(modkit --version 2>&1 | head -n1)"
    echo -e "snakemake\t$(snakemake --version 2>&1 | head -n1)"
    echo -e "conda\t$(conda --version 2>&1 | head -n1)"
} > "${OUTDIR}/00.input/software_versions.tsv"

# --------------------------- discover and check BAMs ---------------------------

log "Discovering HiFi BAM files."

mapfile -t HIFI_BAMS < <(
    find "${HIFI_DIR}" -maxdepth 1 -type f \
        -name 'm84104_*.hifi_reads.bc2012.bam' \
        -print | sort
)

[[ "${#HIFI_BAMS[@]}" -gt 0 ]] || die "No matching HiFi BAM files found in ${HIFI_DIR}"

printf '%s\n' "${HIFI_BAMS[@]}" > "${FOFN}"

log "Found ${#HIFI_BAMS[@]} HiFi BAM files:"
printf '  %s\n' "${HIFI_BAMS[@]}"

{
    echo -e "bam\ttotal_sampled_records\tMM_or_Mm_records\tML_or_Ml_records\tMN_records"
} > "${INPUT_TAG_QC}"

{
    echo -e "bam\tread_records"
} > "${INPUT_READ_COUNTS}"

INPUT_TOTAL_READS=0

for bam in "${HIFI_BAMS[@]}"; do
    # PacBio hifi_reads BAM is an unaligned BAM and may legitimately have no
    # @SQ target sequences in its header.  The -u option tells quickcheck to
    # expect unmapped input while still checking the BAM header and EOF block.
    samtools quickcheck -u -v "${bam}" || \
        die "samtools quickcheck -u failed for unaligned HiFi BAM: ${bam}"

    base="$(basename "${bam}" .bam)"
    tmp_sam="${OUTDIR}/tmp/${base}.first1000.sam"

    sample_bam_records "${bam}" "${tmp_sam}" 1000
    read -r n sampled_mm sampled_ml sampled_mn < <(count_mod_tags_in_sam "${tmp_sam}")
    rm -f "${tmp_sam}"

    echo -e "${bam}\t${n}\t${sampled_mm}\t${sampled_ml}\t${sampled_mn}" >> "${INPUT_TAG_QC}"

    [[ "${n}" -gt 0 ]] || die "No records sampled from ${bam}"
    [[ "${sampled_mm}" -gt 0 ]] || die "No MM/Mm tags detected in the first ${n} records: ${bam}"
    [[ "${sampled_ml}" -gt 0 ]] || die "No ML/Ml tags detected in the first ${n} records: ${bam}"

    bam_count="$(samtools view -c "${bam}")"
    echo -e "${bam}\t${bam_count}" >> "${INPUT_READ_COUNTS}"
    INPUT_TOTAL_READS=$((INPUT_TOTAL_READS + bam_count))
done

echo -e "TOTAL\t${INPUT_TOTAL_READS}" >> "${INPUT_READ_COUNTS}"
log "Total input HiFi read records: ${INPUT_TOTAL_READS}"

# ----------------------------- reference indexing ------------------------------

if [[ "${FORCE}" -eq 1 || ! -s "${REF}.fai" ]]; then
    log "Indexing reference FASTA with samtools faidx."
    samtools faidx "${REF}"
fi

if [[ "${FORCE}" -eq 1 || ! -s "${REF_MMI}" ]]; then
    log "Building pbmm2 HIFI index."
    rm -f "${REF_MMI}"
    pbmm2 index "${REF}" "${REF_MMI}" --preset HIFI \
        2>&1 | tee "${OUTDIR}/logs/pbmm2.index.log"
fi

# -------------------- joint alignment and merged modBAM ------------------------

align_all_hifi() {
    log "Jointly aligning all HiFi BAMs and producing one merged sorted modBAM."
    rm -f "${MERGED_BAM}" "${MERGED_BAM}.bai" "${MERGED_BAM}.csi"

    pbmm2 align \
        "${REF_MMI}" \
        "${FOFN}" \
        "${MERGED_BAM}" \
        --preset HIFI \
        --sort \
        -j "${THREADS}" \
        -J "${SORT_THREADS}" \
        -m "${SORT_MEMORY_PER_THREAD}" \
        --bam-index BAI \
        --log-level INFO \
        2>&1 | tee "${OUTDIR}/logs/pbmm2.align.log"

    samtools quickcheck -v "${MERGED_BAM}"
    [[ -s "${MERGED_BAM}.bai" ]] || samtools index -@ "${THREADS}" "${MERGED_BAM}"
}

run_step "${OUTDIR}/markers/02.pbmm2_align.done" align_all_hifi

# --------------------- retain mapped primary alignments ------------------------

filter_primary_mapped_bam() {
    log "Filtering merged BAM with SAMtools -F 2308."
    log "Retaining mapped primary alignments; removing unmapped, secondary and supplementary records."

    rm -f "${PRIMARY_BAM}" "${PRIMARY_BAM}.bai" "${PRIMARY_BAM}.csi"

    # 2308 = 0x4 (unmapped) + 0x100 (secondary) + 0x800 (supplementary).
    # samtools view preserves the coordinate order of the already sorted input,
    # so an additional samtools sort step is unnecessary. MM/ML tags are retained.
    samtools view \
        -@ "${THREADS}" \
        -bh \
        -F 2308 \
        -o "${PRIMARY_BAM}" \
        "${MERGED_BAM}"

    samtools quickcheck -v "${PRIMARY_BAM}"
    samtools index -@ "${THREADS}" "${PRIMARY_BAM}"
}

run_step "${OUTDIR}/markers/02.filter_primary_mapped.done" filter_primary_mapped_bam

# ------------------------------ alignment QC -----------------------------------

log "Generating alignment and filtering statistics."

samtools flagstat -@ "${THREADS}" "${MERGED_BAM}" > "${ALIGN_FLAGSTAT}"
samtools idxstats "${MERGED_BAM}" > "${ALIGN_IDXSTAT}"

ALL_ALIGNMENT_RECORDS="$(samtools view -c "${MERGED_BAM}")"
UNMAPPED_RECORDS="$(samtools view -c -f 4 "${MERGED_BAM}")"
SECONDARY_RECORDS="$(samtools view -c -f 256 "${MERGED_BAM}")"
SUPPLEMENTARY_RECORDS="$(samtools view -c -f 2048 "${MERGED_BAM}")"
PRIMARY_MAPPED_READS="$(samtools view -c -F 2308 "${MERGED_BAM}")"
FILTERED_BAM_RECORDS="$(samtools view -c "${PRIMARY_BAM}")"

{
    echo -e "category\tcount"
    echo -e "input_hifi_read_records\t${INPUT_TOTAL_READS}"
    echo -e "all_alignment_records\t${ALL_ALIGNMENT_RECORDS}"
    echo -e "unmapped_records_flag_0x4\t${UNMAPPED_RECORDS}"
    echo -e "secondary_records_flag_0x100\t${SECONDARY_RECORDS}"
    echo -e "supplementary_records_flag_0x800\t${SUPPLEMENTARY_RECORDS}"
    echo -e "primary_mapped_records_F2308\t${PRIMARY_MAPPED_READS}"
    echo -e "filtered_bam_records\t${FILTERED_BAM_RECORDS}"
} > "${ALIGN_FILTER_STATS}"

[[ "${FILTERED_BAM_RECORDS}" -eq "${PRIMARY_MAPPED_READS}" ]] || \
    die "Filtered BAM record count does not match samtools view -F 2308 count."

[[ "$(samtools view -c -f 4 "${PRIMARY_BAM}")" -eq 0 ]] || \
    die "Filtered BAM still contains unmapped records."
[[ "$(samtools view -c -f 256 "${PRIMARY_BAM}")" -eq 0 ]] || \
    die "Filtered BAM still contains secondary alignments."
[[ "$(samtools view -c -f 2048 "${PRIMARY_BAM}")" -eq 0 ]] || \
    die "Filtered BAM still contains supplementary alignments."

{
    echo -e "metric\tvalue"
    echo -e "input_read_records\t${INPUT_TOTAL_READS}"
    echo -e "output_primary_mapped_records\t${PRIMARY_MAPPED_READS}"
    echo -e "difference\t$((PRIMARY_MAPPED_READS - INPUT_TOTAL_READS))"
} > "${OUTDIR}/02.alignment/${SAMPLE}.input_vs_output_primary_mapped.tsv"

if [[ "${PRIMARY_MAPPED_READS}" -ne "${INPUT_TOTAL_READS}" ]]; then
    log "WARNING: mapped primary record count differs from input read count. This can occur when reads do not map."
    log "Inspect: ${ALIGN_FILTER_STATS}"
fi

tmp_aligned_sam="${OUTDIR}/tmp/${SAMPLE}.primary_mapped.first1000.sam"
sample_bam_records "${PRIMARY_BAM}" "${tmp_aligned_sam}" 1000
read -r n aligned_mm aligned_ml aligned_mn < <(count_mod_tags_in_sam "${tmp_aligned_sam}")
rm -f "${tmp_aligned_sam}"

{
    echo -e "bam\ttotal_sampled_records\tMM_or_Mm_records\tML_or_Ml_records\tMN_records"
    echo -e "${PRIMARY_BAM}\t${n}\t${aligned_mm}\t${aligned_ml}\t${aligned_mn}"
} > "${ALIGN_TAG_QC}"

[[ "${n}" -gt 0 ]] || die "No records were sampled from the filtered primary BAM."
[[ "${aligned_mm}" -gt 0 ]] || die "MM/Mm tags were not detected in the filtered primary BAM."
[[ "${aligned_ml}" -gt 0 ]] || die "ML/Ml tags were not detected in the filtered primary BAM."

log "Mapped primary BAM passed flag and MM/ML checks."

# --------------------------- prepare target regions ----------------------------

log "Validating and preparing centromere target regions."

awk -v OFS='\t' '
    NR==FNR {
        chrlen[$1]=$2
        next
    }
    /^#/ || NF==0 { next }
    NF < 3 {
        printf("[ERROR] BED line has fewer than 3 columns: %s\n", $0) > "/dev/stderr"
        err=1
        next
    }
    {
        chr=$1
        st=$2
        en=$3

        if (!(chr in chrlen)) {
            printf("[ERROR] Contig not found in FASTA index: %s\n", chr) > "/dev/stderr"
            err=1
            next
        }
        if (st !~ /^[0-9]+$/ || en !~ /^[0-9]+$/ || st < 0 || en <= st) {
            printf("[ERROR] Invalid BED coordinates: %s\n", $0) > "/dev/stderr"
            err=1
            next
        }
        if (en > chrlen[chr]) {
            printf("[ERROR] BED end exceeds contig length: %s; length=%d\n", $0, chrlen[chr]) > "/dev/stderr"
            err=1
            next
        }

        print chr, st, en
    }
    END {
        if (err) exit 1
    }
' "${REF}.fai" "${CENTROMERE_BED}" > "${CENTROMERE_3COL}"

[[ -s "${CENTROMERE_3COL}" ]] || die "No valid centromere intervals were produced."

duplicate_contigs="$(
    cut -f1 "${CENTROMERE_3COL}" | sort | uniq -d | paste -sd, -
)"

if [[ -n "${duplicate_contigs}" ]]; then
    die "CDR-Finder input must contain one target interval per contig. Duplicates: ${duplicate_contigs}"
fi

bedtools slop \
    -i "${CENTROMERE_3COL}" \
    -g "${REF}.fai" \
    -b "${REGION_SLOP}" |
bedtools sort \
    -faidx "${REF}.fai" \
    -i - \
    > "${CENTROMERE_SLOP_BED}"

[[ -s "${CENTROMERE_SLOP_BED}" ]] || die "Failed to create expanded centromere BED."

# ---------------------------- modkit bedMethyl ---------------------------------

make_bedmethyl() {
    log "Generating CpG-combined bedMethyl within expanded centromere regions."
    rm -f "${BEDMETHYL}" "${BEDMETHYL}.tbi"

    modkit pileup \
        "${PRIMARY_BAM}" \
        "${BEDMETHYL}" \
        --modified-bases 5mC \
        --cpg \
        --combine-strands \
        --include-bed "${CENTROMERE_SLOP_BED}" \
        --reference "${REF}" \
        --threads "${THREADS}" \
        --bgzf \
        --log-filepath "${OUTDIR}/logs/modkit.pileup.log"

    tabix -f -p bed "${BEDMETHYL}"
}

run_step "${OUTDIR}/markers/04.modkit_pileup.primary_mapped.done" make_bedmethyl

log "Generating modkit summary for the target regions."

modkit summary \
    "${PRIMARY_BAM}" \
    --reference "${REF}" \
    --include-bed "${CENTROMERE_SLOP_BED}" \
    --mapped-only \
    --tsv \
    > "${MODKIT_SUMMARY}" \
    2> "${OUTDIR}/logs/modkit.summary.log"

# ------------------------------- CDR-Finder ------------------------------------

if [[ ! -s "${CDR_FINDER_DIR}/workflow/Snakefile" ]]; then
    if [[ "${AUTO_CLONE_CDR_FINDER}" -eq 1 ]]; then
        log "CDR-Finder was not found; cloning the official repository."
        mkdir -p "$(dirname "${CDR_FINDER_DIR}")"
        git clone --depth 1 "${CDR_FINDER_GIT}" "${CDR_FINDER_DIR}" \
            2>&1 | tee "${OUTDIR}/logs/git.clone.CDR-Finder.log"
    else
        die "CDR-Finder not found: ${CDR_FINDER_DIR}"
    fi
fi

CDR_FINDER_DIR="$(realpath "${CDR_FINDER_DIR}")"
mkdir -p "${CDR_RESULTS}" "${CDR_LOGS}" "${CDR_BENCHMARKS}"

cat > "${CDR_CONFIG}" <<EOF
samples:
  ${SAMPLE}:
    fasta: "${REF}"
    regions: "${CENTROMERE_SLOP_BED}"
    bam: "${PRIMARY_BAM}"

output_dir: "${CDR_RESULTS}"
log_dir: "${CDR_LOGS}"
benchmark_dir: "${CDR_BENCHMARKS}"

species: "${REPEATMASKER_SPECIES}"

# Panthera pardus centromeres are not human ALR/Alpha arrays.
restrict_alr: false

window_size: ${WINDOW_SIZE}
alr_threshold: 100000
bp_merge: ${BP_MERGE}

height_perc_valley_threshold: ${HEIGHT_THRESHOLD}
prom_perc_valley_threshold: ${PROM_THRESHOLD}
baseline_avg_methyl: ${BASELINE_AVG_METHYL}

bp_alr_merge: 1000

# The supplied region was extended by REGION_SLOP on each side. Excluding the
# same number of bases from the search protects against flank-derived valleys.
bp_edge: ${REGION_SLOP}

extend_edges_std: ${EXTEND_EDGES_STD}
add_hbar: true
overlap_opacity: 0.34
EOF

run_cdr_finder() {
    log "Running CDR-Finder."

    local -a deploy_args
    if snakemake --help 2>&1 | grep -q -- '--software-deployment-method'; then
        deploy_args=(--software-deployment-method conda)
    else
        deploy_args=(--use-conda)
    fi

    snakemake \
        --snakefile "${CDR_FINDER_DIR}/workflow/Snakefile" \
        --configfile "${CDR_CONFIG}" \
        --cores "${THREADS}" \
        "${deploy_args[@]}" \
        --conda-prefix "${OUTDIR}/conda_envs" \
        --rerun-incomplete \
        --latency-wait 120 \
        --printshellcmds \
        2>&1 | tee "${OUTDIR}/logs/CDR-Finder.snakemake.log"
}

run_step "${OUTDIR}/markers/05.CDR_Finder.primary_mapped.done" run_cdr_finder

[[ -e "${CDR_BED}" ]] || die "CDR-Finder finished but the expected CDR BED is missing: ${CDR_BED}"

# ------------------------------- final outputs ---------------------------------

log "Collecting key outputs."

ln -sfn "${PRIMARY_BAM}" "${OUTDIR}/06.final/${SAMPLE}.primary.mapped.modBAM.bam"
ln -sfn "${PRIMARY_BAM}.bai" "${OUTDIR}/06.final/${SAMPLE}.primary.mapped.modBAM.bam.bai"
ln -sfn "${ALIGN_FILTER_STATS}" "${OUTDIR}/06.final/${SAMPLE}.alignment_filter.stats.tsv"
ln -sfn "${CENTROMERE_SLOP_BED}" "${OUTDIR}/06.final/${SAMPLE}.CDR.search_regions.bed"
ln -sfn "${BEDMETHYL}" "${OUTDIR}/06.final/${SAMPLE}.centromere.CpG.bed.gz"
ln -sfn "${BEDMETHYL}.tbi" "${OUTDIR}/06.final/${SAMPLE}.centromere.CpG.bed.gz.tbi"
ln -sfn "${CDR_BED}" "${OUTDIR}/06.final/${SAMPLE}.CDR.bed"

{
    echo -e "metric\tvalue"
    echo -e "hifi_bam_count\t${#HIFI_BAMS[@]}"
    echo -e "input_read_records\t${INPUT_TOTAL_READS}"
    echo -e "output_primary_mapped_records\t${PRIMARY_MAPPED_READS}"

    if [[ -s "${CDR_BED}" ]]; then
        awk '
            BEGIN { n=0; bp=0 }
            NF>=3 && $0 !~ /^#/ {
                n++
                bp += $3-$2
            }
            END {
                printf "cdr_count\t%d\ncdr_total_bp\t%d\n", n, bp
            }
        ' "${CDR_BED}"
    else
        echo -e "cdr_count\t0"
        echo -e "cdr_total_bp\t0"
    fi
} > "${OUTDIR}/06.final/${SAMPLE}.workflow.summary.tsv"

cat > "${OUTDIR}/06.final/README.txt" <<EOF
Key results
===========

Filtered mapped-primary HiFi modBAM (-F 2308):
  ${PRIMARY_BAM}

Unfiltered merged aligned HiFi modBAM:
  ${MERGED_BAM}

Alignment filtering statistics:
  ${ALIGN_FILTER_STATS}

Expanded centromere search regions:
  ${CENTROMERE_SLOP_BED}

CpG-combined bedMethyl:
  ${BEDMETHYL}

CDR-Finder candidate regions:
  ${CDR_BED}

CDR-Finder plots:
  ${CDR_RESULTS}/plot/${SAMPLE}/

Workflow summary:
  ${OUTDIR}/06.final/${SAMPLE}.workflow.summary.tsv

Important interpretation:
  restrict_alr=false was used because Panthera pardus centromeres are not
  human ALR/Alpha arrays. Candidate CDRs should be checked against coverage,
  satellite/HOR structure, and CENP-A enrichment.
EOF

log "Workflow completed successfully."
log "Unfiltered merged modBAM: ${MERGED_BAM}"
log "Filtered mapped-primary modBAM: ${PRIMARY_BAM}"
log "Alignment filter statistics: ${ALIGN_FILTER_STATS}"
log "CpG bedMethyl: ${BEDMETHYL}"
log "CDR BED: ${CDR_BED}"
log "Final directory: ${OUTDIR}/06.final"
