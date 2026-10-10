#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

VERSION="2026-08-18.primary-CDR-v2"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
R_BIN="${R_BIN:-Rscript}"
MERGE_PY="${MERGE_PY:-${SCRIPT_DIR}/cdr_merge_rank.v1.py}"
SAT_R="${SAT_R:-${SCRIPT_DIR}/SAT113_CDR_enrichment.v1.R}"
RUN_SAT113_TEST="${RUN_SAT113_TEST:-1}"
NPERM="${NPERM:-10000}"
RANDOM_SEED="${RANDOM_SEED:-20260818}"
WRITE_WINDOWS="${WRITE_WINDOWS:-0}"

# Primary-CDR defaults. All are environment-overridable.
HARD_GAP="${HARD_GAP:-25000}"
CONDITIONAL_GAP="${CONDITIONAL_GAP:-50000}"
BRIDGE_MIN_FRACTION="${BRIDGE_MIN_FRACTION:-0.50}"
BRIDGE_MIN_DROP="${BRIDGE_MIN_DROP:-8}"
BRIDGE_MAX_RATIO="${BRIDGE_MAX_RATIO:-0.90}"
BRIDGE_MAX_METHYLATION="${BRIDGE_MAX_METHYLATION:-70}"
BASELINE_QUANTILE="${BASELINE_QUANTILE:-0.75}"
FLANK="${FLANK:-100000}"
MIN_BASELINE_WINDOWS="${MIN_BASELINE_WINDOWS:-4}"
MIN_SUPPORTED_BP="${MIN_SUPPORTED_BP:-50000}"
MIN_DENSITY="${MIN_DENSITY:-0.60}"
MIN_DEPLETION="${MIN_DEPLETION:-10}"
MIN_BROAD_SUPPORT="${MIN_BROAD_SUPPORT:-0.50}"
MIN_COVERAGE_FRACTION="${MIN_COVERAGE_FRACTION:-0.50}"
SECONDARY_RATIO="${SECONDARY_RATIO:-0.50}"
MIN_MEAN_SITE_COV="${MIN_MEAN_SITE_COV:-5}"
MODIFIED_MODE="${MODIFIED_MODE:-total}"

log(){ echo "[$(date '+%F %T')] $*" >&2; }
warn(){ echo "[$(date '+%F %T')] [WARNING] $*" >&2; }
die(){ echo "[$(date '+%F %T')] [ERROR] $*" >&2; exit 1; }

usage(){
cat <<'EOF'
Usage:
  Single species:
    bash run_carnivora_primary_CDR_pipeline.v2.sh single \
      SPECIES RAW_CDR METHYL_BED CENTROMERE_BED SAT113_GFF_OR_NA OUTDIR

  Batch manifest:
    bash run_carnivora_primary_CDR_pipeline.v2.sh manifest manifest.tsv

Manifest header/columns:
  species  raw_cdr  methyl_bed  centromere_bed  sat113  outdir

Purpose:
  CDR-Finder *_CDR.bed is treated ONLY as discovery seeds. It is never copied
  directly to the final CDR result. This wrapper performs 25-kb hard merging,
  25-50-kb methylation-aware bridging, 1/5/10-kb multi-scale scoring, confidence
  grading, and depletion-AUC ranking to select Primary/Secondary CDRs.

Important:
  CENTROMERE_BED must be the curated centromere BED (not the +/-500-kb CDR-Finder
  search BED). SAT113 is optional and is used only for annotation/randomization;
  it NEVER contributes to Primary CDR ranking.

Environment overrides include:
  HARD_GAP=25000 CONDITIONAL_GAP=50000 BRIDGE_MIN_FRACTION=0.50
  MIN_SUPPORTED_BP=50000 MIN_DENSITY=0.60 MIN_DEPLETION=10
  SECONDARY_RATIO=0.50 NPERM=10000 RUN_SAT113_TEST=1 WRITE_WINDOWS=0
EOF
}

run_one(){
    local sp="$1" raw="$2" methyl="$3" cen="$4" sat="$5" root="$6"
    [[ -s "${raw}" ]] || die "${sp}: raw CDR-Finder BED missing/empty: ${raw}"
    [[ -s "${methyl}" ]] || die "${sp}: methyl BED missing/empty: ${methyl}"
    [[ -s "${cen}" ]] || die "${sp}: curated centromere BED missing/empty: ${cen}"
    if [[ "${sat}" != "NA" && "${sat}" != "." && -n "${sat}" ]]; then
        [[ -s "${sat}" ]] || die "${sp}: SAT113 file missing/empty: ${sat}"
    else
        sat="NA"
    fi

    local post="${root}/01.primary_selection"
    local satout="${root}/02.sat113_test"
    local final="${root}/03.final"
    mkdir -p "${post}" "${satout}" "${final}" "${root}/logs"

    local -a cmd=(
      "${PYTHON_BIN}" "${MERGE_PY}"
      --raw-cdr "${raw}"
      --methyl-bed "${methyl}"
      --centromere-bed "${cen}"
      --outdir "${post}"
      --prefix "${sp}"
      --modified-mode "${MODIFIED_MODE}"
      --min-mean-site-cov "${MIN_MEAN_SITE_COV}"
      --hard-gap "${HARD_GAP}"
      --conditional-gap "${CONDITIONAL_GAP}"
      --bridge-min-fraction "${BRIDGE_MIN_FRACTION}"
      --bridge-min-drop "${BRIDGE_MIN_DROP}"
      --bridge-max-ratio "${BRIDGE_MAX_RATIO}"
      --bridge-max-methylation "${BRIDGE_MAX_METHYLATION}"
      --baseline-quantile "${BASELINE_QUANTILE}"
      --flank "${FLANK}"
      --min-baseline-windows "${MIN_BASELINE_WINDOWS}"
      --min-supported-bp "${MIN_SUPPORTED_BP}"
      --min-density "${MIN_DENSITY}"
      --min-depletion "${MIN_DEPLETION}"
      --min-broad-support "${MIN_BROAD_SUPPORT}"
      --min-coverage-fraction "${MIN_COVERAGE_FRACTION}"
      --secondary-ratio "${SECONDARY_RATIO}"
    )
    [[ "${sat}" == "NA" ]] || cmd+=(--sat113 "${sat}" --sat-format auto)
    [[ "${WRITE_WINDOWS}" == "1" ]] && cmd+=(--write-windows)

    log "${sp}: selecting final Primary CDRs from CDR-Finder discovery seeds"
    "${cmd[@]}" 2>&1 | tee "${root}/logs/${sp}.primary_selection.log"

    local primary="${post}/${sp}.primary_CDR.bed"
    local graded="${post}/${sp}.CDR.graded.tsv"
    [[ -e "${primary}" ]] || die "${sp}: primary CDR BED was not produced"
    [[ -s "${graded}" ]] || die "${sp}: graded CDR table was not produced"

    # Preserve raw discovery calls explicitly, but make the compatibility name
    # ${sp}.CDR.bed point to the newly selected Primary CDR set.
    ln -sfn "$(realpath "${raw}")" "${final}/${sp}.CDR-Finder.discovery.bed"
    ln -sfn "$(realpath "${methyl}")" "${final}/${sp}.methyl.bed"
    ln -sfn "$(realpath "${cen}")" "${final}/${sp}.centromere.bed"
    ln -sfn "$(realpath "${primary}")" "${final}/${sp}.primary_CDR.bed"
    ln -sfn "$(realpath "${primary}")" "${final}/${sp}.CDR.bed"
    ln -sfn "$(realpath "${post}/${sp}.primary_secondary_CDR.bed")" "${final}/${sp}.primary_secondary_CDR.bed"
    ln -sfn "$(realpath "${graded}")" "${final}/${sp}.CDR.graded.tsv"
    ln -sfn "$(realpath "${post}/${sp}.chromosome_summary.tsv")" "${final}/${sp}.chromosome_summary.tsv"
    ln -sfn "$(realpath "${post}/${sp}.METHODS.txt")" "${final}/${sp}.PRIMARY_CDR_METHODS.txt"
    [[ "${sat}" == "NA" ]] || ln -sfn "$(realpath "${sat}")" "${final}/${sp}.SAT113.annotation"

    {
      echo -e "metric\tvalue"
      echo -e "pipeline_version\t${VERSION}"
      echo -e "raw_discovery_cdr_count\t$(awk '!/^#/&&NF>=3{n++}END{print n+0}' "${raw}")"
      echo -e "merged_candidate_count\t$(awk 'NR>1{n++}END{print n+0}' "${graded}")"
      echo -e "primary_count\t$(awk '!/^#/&&NF>=3{n++}END{print n+0}' "${primary}")"
      awk -F'\t' 'NR>1{a[$6]++}END{for(k in a)print "role_"k"\t"a[k]}' "${graded}" | sort
      awk -F'\t' 'NR>1{a[$7]++}END{for(k in a)print "confidence_"k"\t"a[k]}' "${graded}" | sort
    } > "${final}/${sp}.CDR.summary.tsv"

    if [[ "${sat}" != "NA" && "${RUN_SAT113_TEST}" == "1" ]]; then
        if [[ ! -s "${primary}" ]]; then
            warn "${sp}: no Primary/Primary-uncertain CDRs; SAT113 randomization test skipped."
        elif command -v "${R_BIN}" >/dev/null 2>&1 && [[ -s "${SAT_R}" ]]; then
            log "${sp}: running independent SAT113 within-centromere randomization test"
            "${R_BIN}" "${SAT_R}" \
              --cdr "${primary}" \
              --sat113 "${sat}" \
              --centromere "${cen}" \
              --prefix "${sp}" \
              --outdir "${satout}" \
              --nperm "${NPERM}" \
              --seed "${RANDOM_SEED}" \
              2>&1 | tee "${root}/logs/${sp}.SAT113_test.log"
        else
            warn "${sp}: Rscript or SAT113 R worker unavailable; Primary CDR selection is complete, SAT113 test skipped."
        fi
    fi

    cat > "${final}/README.txt" <<EOF
${sp} final CDR outputs
=======================

${sp}.CDR-Finder.discovery.bed
  Raw CDR-Finder calls retained only as discovery seeds.

${sp}.CDR.bed / ${sp}.primary_CDR.bed
  Final Primary CDR calls after methylation-aware merging, multi-scale validation,
  confidence grading, and depletion-AUC ranking. These are the files to use in
  downstream comparative analyses.

${sp}.primary_secondary_CDR.bed
  Primary plus strong Secondary domains.

${sp}.CDR.graded.tsv
  Complete candidate/QC table. SAT113 columns are annotations only and do not
  affect ranking.
EOF

    log "${sp}: COMPLETE -> ${final}/${sp}.CDR.bed"
}

[[ -s "${MERGE_PY}" ]] || die "Primary CDR worker missing: ${MERGE_PY}"
"${PYTHON_BIN}" -m py_compile "${MERGE_PY}" || die "Python worker syntax check failed"

case "${1:-}" in
  single)
    [[ $# -eq 7 ]] || { usage; exit 2; }
    run_one "$2" "$3" "$4" "$5" "$6" "$7"
    ;;
  manifest)
    [[ $# -eq 2 ]] || { usage; exit 2; }
    manifest="$2"; [[ -s "${manifest}" ]] || die "Manifest missing/empty: ${manifest}"
    first=1
    while IFS=$'\t' read -r sp raw methyl cen sat out extra; do
      sp="${sp//$'\r'/}"; raw="${raw//$'\r'/}"; methyl="${methyl//$'\r'/}"; cen="${cen//$'\r'/}"; sat="${sat//$'\r'/}"; out="${out//$'\r'/}"
      [[ -z "${sp}" || "${sp}" == \#* ]] && continue
      if [[ "${first}" == "1" && "${sp}" == "species" ]]; then first=0; continue; fi
      first=0
      run_one "${sp}" "${raw}" "${methyl}" "${cen}" "${sat}" "${out}"
    done < "${manifest}"
    ;;
  *) usage; exit 2 ;;
esac
