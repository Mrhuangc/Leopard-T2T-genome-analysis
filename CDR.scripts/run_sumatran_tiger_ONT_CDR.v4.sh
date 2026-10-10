#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

PIPELINE_VERSION="2026-08-17.v4.1.1-sumatran-tiger-4khz"

# Sumatran tiger / PTI inputs prepared as four symlinked POD5 run directories.
PROJECT_ROOT="${PROJECT_ROOT:-/path/to/carnivora_project}"
SPECIES="${SPECIES:-PTI}"
POD5_ROOT="${POD5_ROOT:-${PROJECT_ROOT}/huabao/CDR/PTI/pod5}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/huabao/CDR/PTI/output}"

DORADO_HOME="${DORADO_HOME:-/path/to/software/dorado-0.8.2-linux-x64}"
DORADO_BIN="${DORADO_BIN:-${DORADO_HOME}/bin/dorado}"
DORADO_MODELS="${DORADO_MODELS:-/path/to/software/dorado_models_4khz_v4.1.0}"
CDR_SEARCH_ROOT="${CDR_SEARCH_ROOT:-/path/to/software/CDR-Finder}"
CDR_CONDA_PREFIX="${CDR_CONDA_PREFIX:-/path/to/software/CDR-Finder/conda_envs}"

THREADS="${THREADS:-48}"
SAMTOOLS_VIEW_THREADS="${SAMTOOLS_VIEW_THREADS:-8}"
SAMTOOLS_SORT_THREADS="${SAMTOOLS_SORT_THREADS:-16}"
SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-2G}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0,1,2,3,4,5}"
DORADO_DEVICE="${DORADO_DEVICE:-cuda:all}"
# All four PTI runs are R10.4.1/E8.2 400 bps, sampled at 4 kHz. Dorado 2.1.0's
# model catalogue starts at v4.2.0 (5 kHz), so use Dorado 0.8.2 and pin the
# last compatible 4 kHz simplex/modified-base model pair.
DORADO_MODEL="${DORADO_MODEL:-sup@v4.1.0,5mCG_5hmCG@v2}"
REQUIRED_DORADO_VERSION="${REQUIRED_DORADO_VERSION:-0.8.2}"
EXPECTED_SIMPLEX_MODEL="${EXPECTED_SIMPLEX_MODEL:-dna_r10.4.1_e8.2_400bps_sup@v4.1.0}"
EXPECTED_MODIFIED_MODEL="${EXPECTED_MODIFIED_MODEL:-dna_r10.4.1_e8.2_400bps_sup@v4.1.0_5mCG_5hmCG@v2}"

CENTROMERE_SLOP="${CENTROMERE_SLOP:-500000}"
FULL_BAM_QC="${FULL_BAM_QC:-0}"
FORCE_RUN_BASECALL="${FORCE_RUN_BASECALL:-0}"
FORCE_CDR="${FORCE_CDR:-0}"
# Refuse to proceed with a subset of p1-p4 unless explicitly authorized.
ALLOW_PARTIAL_RUNS="${ALLOW_PARTIAL_RUNS:-0}"
# Compare the current POD5 counts with the supplied tree listings.
STRICT_POD5_COUNTS="${STRICT_POD5_COUNTS:-0}"

CDR_SPECIES="${CDR_SPECIES:-Carnivora}"
WINDOW_SIZE="${WINDOW_SIZE:-5000}"
HEIGHT_THRESHOLD="${HEIGHT_THRESHOLD:-0.34}"
PROM_THRESHOLD="${PROM_THRESHOLD:-0.30}"
BASELINE_AVG_METHYL="${BASELINE_AVG_METHYL:-0.40}"
BP_MERGE="${BP_MERGE:-1}"
BP_EDGE="${BP_EDGE:-500000}"
EXTEND_EDGES_STD="${EXTEND_EDGES_STD:--1}"

usage(){
    cat <<'USAGE'
Usage:
  bash run_sumatran_tiger_ONT_CDR.v4.sh GENOME_FASTA [CENTROMERE_BED|AUTO]

Examples:
  bash run_sumatran_tiger_ONT_CDR.v4.sh /path/PTI.fa /path/PTI.centromere.bed
  bash run_sumatran_tiger_ONT_CDR.v4.sh /path/PTI.fa AUTO

Important environment overrides:
  POD5_ROOT=/path/to/carnivora_project/huabao/CDR/PTI/pod5
  OUTDIR=/path/to/carnivora_project/huabao/CDR/PTI/output
  CUDA_VISIBLE_DEVICES=0,1,2,3,4,5
  DORADO_DEVICE=cuda:all
  DORADO_HOME=/path/to/software/dorado-0.8.2-linux-x64
  DORADO_MODEL=sup@v4.1.0,5mCG_5hmCG@v2
  # The defaults above are intentionally pinned for these 4 kHz POD5 files.
  ALLOW_PARTIAL_RUNS=0       # set to 1 only to knowingly run with failed runs
  STRICT_POD5_COUNTS=0       # set to 1 to require p1=384,p2=1571,p3=237,p4=193
  FORCE_RUN_BASECALL=0
  FORCE_CDR=0
USAGE
}

[[ "${1:-}" != "-h" && "${1:-}" != "--help" ]] || { usage; exit 0; }
GENOME_ARG="${GENOME_FASTA:-${1:-}}"
CEN_ARG="${CENTROMERE_BED:-${2:-AUTO}}"

log(){ echo "[$(date '+%F %T')] $*" >&2; }
warn(){ echo "[$(date '+%F %T')] [WARNING] $*" >&2; }
die(){ echo "[$(date '+%F %T')] [ERROR] $*" >&2; exit 1; }
require_cmd(){ command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }
dorado(){ env LD_LIBRARY_PATH="${DORADO_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" "${DORADO_BIN}" "$@"; }

check_dorado_4khz_models(){
    local version catalog clean_catalog
    version="$(dorado --version 2>&1 | tail -1)"
    [[ "${version}" == *"${REQUIRED_DORADO_VERSION}"* ]] ||
      die "This 4 kHz dataset requires Dorado ${REQUIRED_DORADO_VERSION}; found: ${version}. Set DORADO_HOME to dorado-${REQUIRED_DORADO_VERSION}-linux-x64"

    catalog="$(dorado download --list 2>&1)" ||
      die "Cannot read the Dorado ${REQUIRED_DORADO_VERSION} model catalogue"
    # Some Dorado builds decorate spdlog output with terminal control bytes.
    # Strip those bytes before testing; a display-format mismatch must not block
    # a model that Dorado itself can resolve and download in run_pod5_run().
    clean_catalog="$(printf '%s\n' "${catalog}" |
      sed $'s/\033\\[[0-9;?]*[ -\\/]*[@-~]//g' | tr -d '\r')"
    if ! LC_ALL=C grep -aFq -- "${EXPECTED_SIMPLEX_MODEL}" <<< "${clean_catalog}"; then
        warn "Could not text-match ${EXPECTED_SIMPLEX_MODEL} in Dorado catalogue; continuing to Dorado's authoritative model resolution"
    fi
    if ! LC_ALL=C grep -aFq -- "${EXPECTED_MODIFIED_MODEL}" <<< "${clean_catalog}"; then
        warn "Could not text-match ${EXPECTED_MODIFIED_MODEL} in Dorado catalogue; continuing to Dorado's authoritative model resolution"
    fi
}

find_cdr_repo(){
    local x
    for x in "${CDR_SEARCH_ROOT}" "${CDR_SEARCH_ROOT}/CDR-Finder"; do
        [[ -s "${x}/workflow/Snakefile" ]] && { realpath "${x}"; return; }
    done
    x="$(find "${CDR_SEARCH_ROOT}" -maxdepth 3 -type f -path '*/workflow/Snakefile' -print 2>/dev/null | head -1 || true)"
    [[ -n "${x}" ]] && { dirname "$(dirname "${x}")"; return; }
    mkdir -p "${CDR_SEARCH_ROOT}"
    x="${CDR_SEARCH_ROOT}/CDR-Finder"
    [[ ! -e "${x}" ]] || die "${x} exists but has no workflow/Snakefile"
    git clone --depth 1 https://github.com/EichlerLab/CDR-Finder.git "${x}"
    realpath "${x}"
}

resolve_cenbed(){
    local requested="$1" search_root
    local -a hits
    if [[ "${requested}" != "AUTO" ]]; then
        [[ -s "${requested}" ]] || return 1
        realpath "${requested}"
        return
    fi
    search_root="${PROJECT_ROOT}/huabao"
    [[ -d "${search_root}" ]] || search_root="${PROJECT_ROOT}"
    mapfile -t hits < <(
        find "${search_root}" -type f \
          \( -iname 'PTI.centromere.region.bed' \
             -o -iname 'PTI.centromere.bed' \
             -o -iname 'sumatran_tiger.centromere.region.bed' \
             -o -iname 'sumatran_tiger.centromere.bed' \) \
          -print 2>/dev/null | sort -u
    )
    if [[ "${#hits[@]}" -eq 1 ]]; then
        realpath "${hits[0]}"
        return
    fi
    if [[ "${#hits[@]}" -eq 0 ]]; then
        echo "No PTI/Sumatran tiger centromere BED found automatically under ${search_root}" >&2
    else
        echo "Multiple candidate centromere BEDs found; provide the exact path:" >&2
        printf '  %s\n' "${hits[@]}" >&2
    fi
    return 1
}

sample_mm_ml(){
    local bam="$1" out="$2" rc
    set +o pipefail
    samtools view "${bam}" | head -n 1000 | awk -F'\t' '
    {mm=0;ml=0; for(i=12;i<=NF;i++){if($i~/^(MM|Mm):Z:/)mm=1;if($i~/^(ML|Ml):B:/)ml=1} n++;a+=mm;b+=ml}
    END{print "metric\tvalue";print "sampled_records\t"n+0;print "MM_records\t"a+0;print "ML_records\t"b+0;if(n==0||a==0||b==0)exit 7}' > "${out}"
    rc=${PIPESTATUS[2]}
    set -o pipefail
    return "${rc}"
}

prepare_regions(){
    local genome="$1" cen="$2"
    local d bed3 slop dup
    d="${OUTDIR}/04.regions"
    bed3="${d}/${SPECIES}.centromere.3col.bed"
    slop="${d}/${SPECIES}.centromere.plus${CENTROMERE_SLOP}.bed"
    mkdir -p "${d}"
    [[ -s "${genome}.fai" ]] || samtools faidx "${genome}"
    awk -v OFS='\t' 'NR==FNR{L[$1]=$2;next}/^#/||NF==0{next}{if(NF<3||!($1 in L)||$2<0||$3<=$2||$3>L[$1]){print "BAD BED:",$0>"/dev/stderr";bad=1;next}print $1,$2,$3}END{if(bad)exit 2}' "${genome}.fai" "${cen}" > "${bed3}"
    [[ -s "${bed3}" ]] || die "${SPECIES}: no valid centromere intervals"
    dup="$(cut -f1 "${bed3}" | sort | uniq -d | paste -sd, -)"
    [[ -z "${dup}" ]] || die "${SPECIES}: multiple target intervals on contig(s): ${dup}"
    bedtools slop -i "${bed3}" -g "${genome}.fai" -b "${CENTROMERE_SLOP}" |
      bedtools sort -faidx "${genome}.fai" -i - > "${slop}"
    echo "${slop}"
}

run_pod5_run(){
    local id="$1" poddir="$2" genome="$3"
    local bam marker first tmp pod5_count
    bam="${OUTDIR}/02.run_bams/${id}.primary.mapped.sorted.bam"
    marker="${OUTDIR}/markers/run.${id}.done"
    mkdir -p "${OUTDIR}"/{02.run_bams,03.run_qc,logs,markers,tmp} "${DORADO_MODELS}"

    if [[ "${FORCE_RUN_BASECALL}" -eq 0 && -s "${marker}" && -s "${bam}" ]] &&
       samtools quickcheck -v "${bam}" >/dev/null 2>&1; then
        log "${SPECIES}/${id}: already complete"
        return
    fi

    poddir="$(realpath "${poddir}")"
    [[ -d "${poddir}" ]] || die "${SPECIES}/${id}: POD5 directory missing: ${poddir}"
    pod5_count="$(find "${poddir}" -type f -name '*.pod5' | wc -l)"
    [[ "${pod5_count}" -gt 0 ]] || die "${SPECIES}/${id}: no POD5 files"
    first="$(find "${poddir}" -type f -name '*.pod5' -print -quit)"

    pod5 inspect summary "${first}" \
      > "${OUTDIR}/03.run_qc/${id}.pod5.summary.txt" \
      2> "${OUTDIR}/logs/${id}.pod5.inspect.log"

    log "${SPECIES}/${id}: resolving ${DORADO_MODEL} for ${pod5_count} POD5 files"
    dorado download \
      --model "${DORADO_MODEL}" \
      --data "${poddir}" \
      --recursive \
      --models-directory "${DORADO_MODELS}" \
      > "${OUTDIR}/logs/${id}.dorado.model_download.log" 2>&1 ||
        die "${SPECIES}/${id}: Dorado model resolution failed; inspect ${OUTDIR}/logs/${id}.dorado.model_download.log"

    tmp="${bam}.partial"
    rm -f "${tmp}" "${bam}" "${bam}.bai"
    log "${SPECIES}/${id}: Dorado ${DORADO_MODEL}; devices=${CUDA_VISIBLE_DEVICES}"
    dorado basecaller "${DORADO_MODEL}" "${poddir}" \
      --recursive \
      --device "${DORADO_DEVICE}" \
      --reference "${genome}" \
      --models-directory "${DORADO_MODELS}" \
      2> "${OUTDIR}/logs/${id}.dorado.log" |
      samtools view -@ "${SAMTOOLS_VIEW_THREADS}" -bh -F 2308 - \
        2> "${OUTDIR}/logs/${id}.samtools.view.log" |
      samtools sort -@ "${SAMTOOLS_SORT_THREADS}" -m "${SAMTOOLS_SORT_MEM}" \
        -T "${OUTDIR}/tmp/${id}.sort" -o "${tmp}" - \
        2> "${OUTDIR}/logs/${id}.samtools.sort.log"

    samtools quickcheck -v "${tmp}" || die "${SPECIES}/${id}: BAM quickcheck failed"
    mv "${tmp}" "${bam}"
    samtools index -@ "${THREADS}" "${bam}"
    sample_mm_ml "${bam}" "${OUTDIR}/03.run_qc/${id}.MM_ML.tsv" ||
      die "${SPECIES}/${id}: MM/ML tags not detected"
    samtools view -H "${bam}" | awk '$1=="@RG"' > "${OUTDIR}/03.run_qc/${id}.RG.header.txt"
    samtools idxstats "${bam}" > "${OUTDIR}/03.run_qc/${id}.idxstats.tsv"
    printf 'run_id\tpod5_dir\tpod5_count\tcompleted\n%s\t%s\t%s\t%s\n' \
      "${id}" "${poddir}" "${pod5_count}" "$(date '+%F %T')" \
      > "${OUTDIR}/03.run_qc/${id}.run_input.tsv"
    date '+%F %T' > "${marker}"
}

merge_runs(){
    local merged id bam marker
    local -a bams
    merged="${OUTDIR}/05.merged/${SPECIES}.ONT.5mCG.primary.mapped.sorted.bam"
    mkdir -p "${OUTDIR}/05.merged"
    bams=()
    # Only merge BAMs with a successful run marker. This prevents an old or
    # interrupted BAM from being included when ALLOW_PARTIAL_RUNS=1.
    for id in "${RUN_IDS[@]}"; do
        bam="${OUTDIR}/02.run_bams/${id}.primary.mapped.sorted.bam"
        marker="${OUTDIR}/markers/run.${id}.done"
        if [[ -s "${marker}" && -s "${bam}" ]] && samtools quickcheck -v "${bam}" >/dev/null 2>&1; then
            bams+=("${bam}")
        fi
    done
    [[ "${#bams[@]}" -gt 0 ]] || die "${SPECIES}: no completed run BAMs"
    printf '%s\n' "${bams[@]}" > "${OUTDIR}/05.merged/${SPECIES}.run_bams.list"
    if [[ "${#bams[@]}" -eq 1 ]]; then
        ln -sfn "${bams[0]}" "${merged}"
        ln -sfn "${bams[0]}.bai" "${merged}.bai"
    else
        rm -f "${merged}" "${merged}.bai"
        samtools merge -@ "${THREADS}" -f -o "${merged}" "${bams[@]}" \
          2> "${OUTDIR}/logs/${SPECIES}.samtools.merge.log"
        samtools quickcheck -v "${merged}" || die "${SPECIES}: merged BAM quickcheck failed"
        samtools index -@ "${THREADS}" "${merged}"
    fi
    sample_mm_ml "${merged}" "${OUTDIR}/05.merged/${SPECIES}.MM_ML.tsv" ||
      die "${SPECIES}: merged BAM lacks MM/ML tags"
    samtools idxstats "${merged}" > "${OUTDIR}/05.merged/${SPECIES}.idxstats.tsv"
    [[ "${FULL_BAM_QC}" -eq 0 ]] ||
      samtools flagstat -@ "${THREADS}" "${merged}" > "${OUTDIR}/05.merged/${SPECIES}.flagstat.txt"
    echo "${merged}"
}

run_cdr(){
    local genome="$1" regions="$2" bam="$3" repo="$4"
    local root res logs bench cfg cdr tools
    local -a dep args
    root="${OUTDIR}/06.cdr_finder"
    res="${root}/results"
    logs="${root}/logs"
    bench="${root}/benchmarks"
    cfg="${root}/config.${SPECIES}.yaml"
    cdr="${res}/bed/${SPECIES}_CDR.bed"
    mkdir -p "${root}" "${res}" "${logs}" "${bench}" "${CDR_CONDA_PREFIX}"
    cat > "${cfg}" <<CFG
samples:
  ${SPECIES}:
    fasta: "${genome}"
    regions: "${regions}"
    bam: "${bam}"
output_dir: "${res}"
log_dir: "${logs}"
benchmark_dir: "${bench}"
species: "${CDR_SPECIES}"
restrict_alr: false
window_size: ${WINDOW_SIZE}
alr_threshold: 100000
bp_merge: ${BP_MERGE}
height_perc_valley_threshold: ${HEIGHT_THRESHOLD}
prom_perc_valley_threshold: ${PROM_THRESHOLD}
baseline_avg_methyl: ${BASELINE_AVG_METHYL}
bp_alr_merge: 1000
bp_edge: ${BP_EDGE}
extend_edges_std: ${EXTEND_EDGES_STD}
add_hbar: true
overlap_opacity: 0.34
CFG
    tools="${repo}/workflow/envs/tools.yaml"
    grep -qE '^[[:space:]]*-[[:space:]]*ont-modkit' "${tools}" ||
      die "No ont-modkit entry in CDR-Finder tools.yaml"
    sed -i -E 's|^([[:space:]]*-[[:space:]]*ont-modkit).*$|\1==0.3.2|' "${tools}"
    if [[ "${FORCE_CDR}" -eq 0 && -e "${cdr}" && -e "${res}/plot/${SPECIES}.done" ]]; then
        echo "${cdr}"
        return
    fi
    if snakemake --help 2>&1 | grep -q -- '--software-deployment-method'; then
        dep=(--software-deployment-method conda)
    else
        dep=(--use-conda)
    fi
    args=(
      --snakefile "${repo}/workflow/Snakefile"
      --configfile "${cfg}"
      --workflow-profile none
      --cores "${THREADS}"
      "${dep[@]}"
      --conda-prefix "${CDR_CONDA_PREFIX}"
      --rerun-incomplete
      --latency-wait 120
      --printshellcmds
    )
    snakemake "${args[@]}" --dry-run 2>&1 |
      tee "${OUTDIR}/logs/${SPECIES}.CDR-Finder.dryrun.log" >&2
    snakemake "${args[@]}" 2>&1 |
      tee "${OUTDIR}/logs/${SPECIES}.CDR-Finder.log" >&2
    [[ -e "${cdr}" ]] || die "${SPECIES}: CDR BED not produced"
    echo "${cdr}"
}

finalize(){
    local bam="$1" regions="$2" cdr="$3"
    local f meth
    f="${OUTDIR}/07.final"
    mkdir -p "${f}"
    ln -sfn "${bam}" "${f}/${SPECIES}.ONT.5mCG.primary.mapped.sorted.bam"
    [[ ! -s "${bam}.bai" ]] ||
      ln -sfn "${bam}.bai" "${f}/${SPECIES}.ONT.5mCG.primary.mapped.sorted.bam.bai"
    ln -sfn "${regions}" "${f}/${SPECIES}.CDR.search_regions.bed"
    ln -sfn "${cdr}" "${f}/${SPECIES}.CDR.bed"
    meth="${OUTDIR}/06.cdr_finder/results/bed/${SPECIES}_methyl.bed"
    [[ ! -s "${meth}" ]] || ln -sfn "${meth}" "${f}/${SPECIES}.methyl.bed"
    awk 'BEGIN{OFS="\t";print "metric","value"}!/^#/&&NF>=3{n++;bp+=$3-$2}END{print "cdr_count",n+0;print "cdr_total_bp",bp+0}' \
      "${cdr}" > "${f}/${SPECIES}.CDR.summary.tsv"
}

for c in samtools bedtools pod5 snakemake conda git find awk sed sort uniq grep realpath nvidia-smi; do
    require_cmd "${c}"
done
[[ -n "${GENOME_ARG}" ]] || { usage >&2; die "GENOME_FASTA is required"; }
[[ -x "${DORADO_BIN}" ]] || die "Dorado not executable: ${DORADO_BIN}"
[[ -d "${DORADO_HOME}/lib" ]] || die "Dorado lib missing: ${DORADO_HOME}/lib"
[[ -d "${POD5_ROOT}" ]] || die "POD5 root missing: ${POD5_ROOT}"
[[ -s "${GENOME_ARG}" ]] || die "Genome FASTA missing or empty: ${GENOME_ARG}"
GENOME="$(realpath "${GENOME_ARG}")"
CENBED="$(resolve_cenbed "${CEN_ARG}")" ||
  die "Cannot resolve centromere BED; provide its exact path as argument 2"

mkdir -p "${OUTDIR}"/{00.input,02.run_bams,03.run_qc,04.regions,05.merged,06.cdr_finder,07.final,logs,markers,tmp} \
  "${DORADO_MODELS}" "${CDR_CONDA_PREFIX}"

if ! env -u LD_LIBRARY_PATH -u LD_PRELOAD nvidia-smi >/dev/null 2>&1; then
    die "nvidia-smi failed; fix NVIDIA driver/NVML before running Dorado"
fi
GPU_COUNT="$(env -u LD_LIBRARY_PATH -u LD_PRELOAD nvidia-smi -L | grep -c '^GPU ' || true)"
[[ "${GPU_COUNT}" -gt 0 ]] || die "No NVIDIA GPU detected"
check_dorado_4khz_models

declare -A EXPECTED_COUNTS=( [p1]=384 [p2]=1571 [p3]=237 [p4]=193 )
declare -a RUN_IDS=() RUN_DIRS=() FAILED_RUNS=()
RUN_MANIFEST="${OUTDIR}/00.input/${SPECIES}.pod5_runs.tsv"
printf 'run_id\tsymlink_path\tresolved_pod5_dir\tpod5_count\texpected_count\n' > "${RUN_MANIFEST}"

for link in "${POD5_ROOT}"/*; do
    [[ -d "${link}" ]] || continue
    id="$(basename "${link}")"
    [[ "${id}" =~ ^p[1-4]$ ]] || { warn "Ignoring unexpected directory: ${link}"; continue; }
    target="$(realpath "${link}")"
    count="$(find "${target}" -type f -name '*.pod5' | wc -l)"
    [[ "${count}" -gt 0 ]] || die "${id}: no POD5 files under ${target}"
    expected="${EXPECTED_COUNTS[${id}]}"
    printf '%s\t%s\t%s\t%s\t%s\n' "${id}" "${link}" "${target}" "${count}" "${expected}" >> "${RUN_MANIFEST}"
    if [[ "${count}" -ne "${expected}" ]]; then
        if [[ "${STRICT_POD5_COUNTS}" -eq 1 ]]; then
            die "${id}: POD5 count ${count}, expected ${expected} from supplied tree"
        fi
        warn "${id}: POD5 count ${count}, supplied tree listed ${expected}; continuing"
    fi
    RUN_IDS+=("${id}")
    RUN_DIRS+=("${target}")
done

[[ "${#RUN_IDS[@]}" -eq 4 ]] ||
  die "Expected four run links p1-p4 under ${POD5_ROOT}; found ${#RUN_IDS[@]}"

printf 'pipeline_version\t%s\nspecies\t%s\ngenome\t%s\ncentromere_bed\t%s\npod5_root\t%s\noutput\t%s\n' \
  "${PIPELINE_VERSION}" "${SPECIES}" "${GENOME}" "${CENBED}" "$(realpath "${POD5_ROOT}")" "${OUTDIR}" \
  > "${OUTDIR}/00.input/${SPECIES}.resolved_inputs.tsv"

log "Pipeline=${PIPELINE_VERSION}"
log "Species=${SPECIES}; POD5=${POD5_ROOT}; genome=${GENOME}; centromere=${CENBED}"
log "GPUs=${GPU_COUNT}; CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}; Dorado=$(dorado --version 2>&1 | tail -1)"
log "4 kHz models: ${EXPECTED_SIMPLEX_MODEL}; ${EXPECTED_MODIFIED_MODEL}"
CDR_REPO="$(find_cdr_repo)"
log "CDR-Finder=${CDR_REPO}; output=${OUTDIR}"

for i in "${!RUN_IDS[@]}"; do
    id="${RUN_IDS[${i}]}"
    poddir="${RUN_DIRS[${i}]}"
    set +e
    ( set -Eeuo pipefail; run_pod5_run "${id}" "${poddir}" "${GENOME}" )
    rc=$?
    set -e
    if [[ "${rc}" -ne 0 ]]; then
        FAILED_RUNS+=("${id}")
        warn "${id}: basecalling failed with exit ${rc}"
    fi
done

if [[ "${#FAILED_RUNS[@]}" -gt 0 ]]; then
    printf '%s\n' "${FAILED_RUNS[@]}" > "${OUTDIR}/00.input/${SPECIES}.failed_runs.list"
    if [[ "${ALLOW_PARTIAL_RUNS}" -ne 1 ]]; then
        die "Failed runs: ${FAILED_RUNS[*]}. Fix them and rerun, or explicitly set ALLOW_PARTIAL_RUNS=1"
    fi
    warn "Proceeding with partial data by explicit request; failed runs: ${FAILED_RUNS[*]}"
else
    rm -f "${OUTDIR}/00.input/${SPECIES}.failed_runs.list"
fi

MERGED="$(merge_runs)"
REGIONS="$(prepare_regions "${GENOME}" "${CENBED}")"
CDR="$(run_cdr "${GENOME}" "${REGIONS}" "${MERGED}" "${CDR_REPO}")"
finalize "${MERGED}" "${REGIONS}" "${CDR}"
date '+%F %T' > "${OUTDIR}/markers/SPECIES.COMPLETE"
log "${SPECIES}: COMPLETE -> ${OUTDIR}/07.final"
