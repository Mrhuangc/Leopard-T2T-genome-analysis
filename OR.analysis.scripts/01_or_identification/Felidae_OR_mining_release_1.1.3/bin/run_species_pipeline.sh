#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash bin/run_species_pipeline.sh --params params.Clu_fa.yml [options]

Options:
  --profile STR       Nextflow profiles (default: slurm,conda)
  --scheduler MODE    optimized or legacy (default: optimized)
  --from-stage N      First compatibility stage, 1-19 (default: 1)
  --to-stage N        Last compatibility stage, 1-19 (default: 19)
  --no-resume         Do not pass -resume
  --prepare-conda-only
                      Pre-create all process environments on this login node and stop
  --initial-only      Generate/validate the initial merged TBLASTN and stop
  --preflight-only    Validate inputs and stop
  --list-stages       Print biological task groups and compatibility stages
  -h, --help          Show this help

The optimized scheduler submits one Slurm allocation for each biological/software
task group. Only the naturally scattered Exonerate searches retain one Slurm task
per search cell. The legacy scheduler preserves the 1.0.1 task-per-process layout.
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$ROOT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
PARAMS=""
PROFILE="slurm,conda"
SCHEDULER="optimized"
FROM_STAGE=1
TO_STAGE=19
RESUME="-resume"
PREFLIGHT_ONLY=false
PREPARE_CONDA_ONLY=false
INITIAL_ONLY=false

STAGES=(
    "validate_merged.nf"
    "build_candidate_regions.nf"
    "orf_discovery.nf"
    "gpcr_filter.nf"
    "phylogenetic_filter.nf"
    "nc_r_parity.nf"
    "species_or_expansion.nf"
    "multiexon_regions.nf"
    "multiexon_extraction.nf"
    "exonerate_discovery.nf"
    "exonerate_model_selection.nf"
    "strict_refinement.nf"
    "multiexon_classification.nf"
    "pre7tm_catalog.nf"
    "interproscan_candidates.nf"
    "final_homology_filter.nf"
    "single_exon_nonfunctional_recovery.nf"
    "or_clade_classification.nf"
    "single_species_or_analysis.nf"
)

# name|first stage|last stage|execution|cpus|memory GiB|walltime hours
# 'bundle' means one Slurm allocation with sequential local Nextflow processes.
# 'distributed' preserves Slurm-per-cell execution for expensive scatter searches.
TASK_GROUPS=(
    "candidate_locus_discovery|1|3|bundle|8|32|12"
    "gpcr_homology_screen|4|4|bundle|8|32|12"
    "or_phylogenetic_validation|5|6|bundle|16|64|36"
    "species_specific_or_expansion|7|7|bundle|16|64|24"
    "multiexon_candidate_definition|8|9|bundle|16|64|24"
    "targeted_exonerate_discovery|10|10|distributed|0|0|0"
    "multiexon_model_selection|11|11|bundle|4|16|4"
    "strict_exonerate_refinement|12|12|distributed|0|0|0"
    "multiexon_catalog_integration|13|14|bundle|4|16|6"
    "seven_tm_domain_validation|15|15|bundle|8|32|8"
    "final_or_homology_validation|16|16|bundle|8|24|6"
    "single_exon_residual_recovery|17|17|bundle|8|32|8"
    "or_clade_and_locus_architecture|18|19|bundle|8|32|8"
)

list_stages() {
    printf 'GROUP\tSTAGES\tEXECUTION\tRESOURCE\n'
    printf 'initial_reference_or_search\tPRE\tdistributed\t%s shards x %s CPU, max %s forks\n' \
        "${INITIAL_QUERY_SHARDS_DISPLAY:-2}" \
        "${INITIAL_TBLASTN_CPUS_DISPLAY:-8}" \
        "${INITIAL_TBLASTN_FORKS_DISPLAY:-2}"
    local row name first last mode cpus memory hours
    for row in "${TASK_GROUPS[@]}"; do
        IFS='|' read -r name first last mode cpus memory hours <<<"$row"
        if [[ "$mode" == "distributed" ]]; then
            printf '%s\t%02d-%02d\t%s\tprocess-specific\n' "$name" "$first" "$last" "$mode"
        else
            printf '%s\t%02d-%02d\t%s\t%s CPU/%s GiB/%s h\n' \
                "$name" "$first" "$last" "$mode" "$cpus" "$memory" "$hours"
        fi
    done
    printf '\nCompatibility stage map:\n'
    local index
    for index in "${!STAGES[@]}"; do
        printf '%02d\t%s\n' "$((index + 1))" "${STAGES[$index]}"
    done
}

while (($#)); do
    case "$1" in
        --params) PARAMS="${2:-}"; shift 2 ;;
        --profile) PROFILE="${2:-}"; shift 2 ;;
        --scheduler) SCHEDULER="${2:-}"; shift 2 ;;
        --from-stage) FROM_STAGE="${2:-}"; shift 2 ;;
        --to-stage) TO_STAGE="${2:-}"; shift 2 ;;
        --no-resume) RESUME=""; shift ;;
        --prepare-conda-only) PREPARE_CONDA_ONLY=true; shift ;;
        --initial-only) INITIAL_ONLY=true; shift ;;
        --preflight-only) PREFLIGHT_ONLY=true; shift ;;
        --list-stages) list_stages; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$PARAMS" ]] || { echo "ERROR: --params is required" >&2; exit 2; }
[[ -f "$PARAMS" ]] || { echo "ERROR: parameter file not found: $PARAMS" >&2; exit 2; }
PARAMS="$(cd "$(dirname "$PARAMS")" && pwd)/$(basename "$PARAMS")"
[[ "$FROM_STAGE" =~ ^[0-9]+$ && "$TO_STAGE" =~ ^[0-9]+$ ]] || {
    echo "ERROR: stage values must be integers" >&2; exit 2;
}
((FROM_STAGE >= 1 && TO_STAGE <= 19 && FROM_STAGE <= TO_STAGE)) || {
    echo "ERROR: require 1 <= from-stage <= to-stage <= 19" >&2; exit 2;
}
[[ "$SCHEDULER" == "optimized" || "$SCHEDULER" == "legacy" ]] || {
    echo "ERROR: --scheduler must be optimized or legacy" >&2; exit 2;
}

cd "$ROOT_DIR"
mkdir -p pipeline_info
ORIGINAL_PARAMS="$PARAMS"
PREFLIGHT_REPORT="pipeline_info/$(basename "$ORIGINAL_PARAMS").preflight.json"
python -m or_mining preflight "$PARAMS" --report "$PREFLIGHT_REPORT"
$PREFLIGHT_ONLY && exit 0

# Resolve environments while this controller/login node still has network access.
# Every release and species reuses the same shared cache. Grouped compute jobs are
# switched to offline Nextflow mode by run_stage_bundle.sh.
if [[ ",$PROFILE," == *,conda,* ]]; then
    PARAMS_CONDA_CACHE="$(
        python - "$PARAMS" <<'PY'
from pathlib import Path
import sys
from or_mining.preflight import load_params

value = load_params(Path(sys.argv[1])).get("conda_cache_dir")
print("" if value is None else str(value))
PY
    )"
    CONDA_CACHE_DIR="${PARAMS_CONDA_CACHE:-${NXF_CONDA_CACHEDIR:-$HOME/.cache/felidae-or/conda}}"
    [[ "$CONDA_CACHE_DIR" = /* ]] || {
        echo "ERROR: conda_cache_dir must be an absolute shared-filesystem path" >&2
        exit 2
    }
    mkdir -p "$CONDA_CACHE_DIR"
    export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"
    echo
    echo "===== LOGIN-NODE CONDA PREWARM ====="
    echo "Shared cache: $NXF_CONDA_CACHEDIR"
    python -m or_mining prewarm-conda "$ROOT_DIR" "$NXF_CONDA_CACHEDIR"
elif $PREPARE_CONDA_ONLY; then
    echo "ERROR: --prepare-conda-only requires a profile containing conda" >&2
    exit 2
fi

$PREPARE_CONDA_ONLY && exit 0

# Resolve the optional initial search before any compatibility stage reads input_merged.
# Existing 1.1.2 parameter files default to reuse and therefore follow the old path unchanged.
mapfile -t INITIAL_VALUES < <(
    python - "$ORIGINAL_PARAMS" <<'PY'
from pathlib import Path
import sys

from or_mining.initial_search import resolve_merged_input
from or_mining.preflight import load_params

p = load_params(Path(sys.argv[1]))
merged, source = resolve_merged_input(p)
print(str(p.get("species") or ""))
print(str(p.get("outdir") or "results"))
print(source)
print(str(merged))
PY
)
INITIAL_SPECIES="${INITIAL_VALUES[0]}"
INITIAL_OUTDIR="${INITIAL_VALUES[1]}"
INITIAL_SOURCE="${INITIAL_VALUES[2]}"
INITIAL_MERGED="${INITIAL_VALUES[3]}"

if [[ "$INITIAL_SOURCE" == "generated" ]]; then
    echo
    echo "===== INITIAL REFERENCE-OR TBLASTN (stage 0) ====="
    echo "Output: $INITIAL_MERGED"
    INITIAL_INFO="$INITIAL_OUTDIR/00_initial_search/pipeline_info"
    mkdir -p "$INITIAL_INFO"
    initial_command=(
        nextflow -log "$INITIAL_INFO/nextflow.log"
        run ./tests/nextflow/initial_tblastn.nf
        -profile "$PROFILE"
        -params-file "$ORIGINAL_PARAMS"
        -with-report "$INITIAL_INFO/execution_report.html"
        -with-trace "$INITIAL_INFO/trace.tsv"
        -with-timeline "$INITIAL_INFO/timeline.html"
        -with-dag "$INITIAL_INFO/workflow_dag.html"
    )
    [[ -n "$RESUME" ]] && initial_command+=("$RESUME")
    "${initial_command[@]}"
    python -m or_mining check-stage "$ORIGINAL_PARAMS" 0
else
    echo
    echo "===== INITIAL REFERENCE-OR TBLASTN ====="
    echo "Reusing validated input: $INITIAL_MERGED"
fi

RESOLVED_PARAMS="$INITIAL_OUTDIR/pipeline_info/${INITIAL_SPECIES}.resolved.params.yml"
python -m or_mining write-resolved-params \
    "$ORIGINAL_PARAMS" \
    --output "$RESOLVED_PARAMS"
PARAMS="$(cd "$(dirname "$RESOLVED_PARAMS")" && pwd)/$(basename "$RESOLVED_PARAMS")"
RESOLVED_PREFLIGHT="$INITIAL_OUTDIR/pipeline_info/${INITIAL_SPECIES}.resolved.preflight.json"
python -m or_mining preflight "$PARAMS" --report "$RESOLVED_PREFLIGHT"

if $INITIAL_ONLY; then
    echo "INITIAL SEARCH COMPLETED"
    echo "Merged input: $INITIAL_MERGED"
    echo "Resolved parameters: $PARAMS"
    exit 0
fi

run_compatibility_stage() {
    local stage="$1"
    local script="tests/nextflow/${STAGES[$((stage - 1))]}"
    echo
    echo "===== COMPATIBILITY STAGE $(printf '%02d' "$stage")/19: ${STAGES[$((stage - 1))]} ====="
    local command=(nextflow run "$script" -profile "$PROFILE" -params-file "$PARAMS")
    [[ -n "$RESUME" ]] && command+=("$RESUME")
    "${command[@]}"
    python -m or_mining check-stage "$PARAMS" "$stage"
}

if [[ "$SCHEDULER" == "legacy" ]]; then
    for ((stage = FROM_STAGE; stage <= TO_STAGE; stage++)); do
        run_compatibility_stage "$stage"
    done
    echo
    echo "PIPELINE COMPLETED THROUGH STAGE $TO_STAGE (legacy scheduler)"
    exit 0
fi

[[ ",$PROFILE," == *,slurm,* ]] || {
    echo "ERROR: optimized scheduler requires the slurm profile" >&2
    echo "Use --scheduler legacy for a non-Slurm run." >&2
    exit 2
}
command -v sbatch >/dev/null || { echo "ERROR: sbatch is not available" >&2; exit 127; }

mapfile -t SCHEDULER_VALUES < <(
    python - "$PARAMS" <<'PY'
from pathlib import Path
import sys
from or_mining.preflight import load_params

p = load_params(Path(sys.argv[1]))
for key in (
    "species",
    "outdir",
    "slurm_partition",
    "slurm_account",
    "slurm_qos",
    "publish_mode",
    "max_cpus",
    "max_memory_gb",
    "max_time_h",
):
    value = p.get(key)
    print("" if value is None else str(value))
PY
)
SPECIES="${SCHEDULER_VALUES[0]}"
OUTDIR="${SCHEDULER_VALUES[1]}"
PARTITION="${SCHEDULER_VALUES[2]}"
ACCOUNT="${SCHEDULER_VALUES[3]}"
QOS="${SCHEDULER_VALUES[4]}"
PUBLISH_MODE="${SCHEDULER_VALUES[5]:-link}"
MAX_CPUS="${SCHEDULER_VALUES[6]:-16}"
MAX_MEMORY_GB="${SCHEDULER_VALUES[7]:-128}"
MAX_TIME_H="${SCHEDULER_VALUES[8]:-24}"
[[ "$ACCOUNT" == "REPLACE_ME" ]] && ACCOUNT=""
[[ "$QOS" == "REPLACE_ME" ]] && QOS=""
mkdir -p "$OUTDIR/pipeline_info/slurm"

if [[ "$PUBLISH_MODE" == "link" ]]; then
    python - "$ROOT_DIR" "$OUTDIR" <<'PY'
import os
from pathlib import Path
import sys
import tempfile

source_dir = Path(sys.argv[1]) / "pipeline_info"
target_dir = Path(sys.argv[2]) / "pipeline_info"
source_dir.mkdir(parents=True, exist_ok=True)
target_dir.mkdir(parents=True, exist_ok=True)
source_fd, source_name = tempfile.mkstemp(prefix=".hardlink_probe.", dir=source_dir)
os.close(source_fd)
target_name = target_dir / (Path(source_name).name + ".link")
try:
    os.link(source_name, target_name)
except OSError as exc:
    raise SystemExit(
        "ERROR: publish_mode=link requires work and outdir on one filesystem; "
        f"hard-link probe failed: {exc}. Set publish_mode: copy."
    )
finally:
    Path(source_name).unlink(missing_ok=True)
    target_name.unlink(missing_ok=True)
PY
else
    echo "WARNING: publish_mode=copy duplicates intermediate payload I/O; use link when possible." >&2
fi

# Replace only the executor profile. Environment profile(s) remain unchanged.
BUNDLE_PROFILE="$(python - "$PROFILE" <<'PY'
import sys
profiles = [item.strip() for item in sys.argv[1].split(",") if item.strip()]
profiles = ["batch" if item == "slurm" else item for item in profiles]
seen = []
for item in profiles:
    if item not in seen:
        seen.append(item)
print(",".join(seen))
PY
)"

submit_bundle() {
    local name="$1" first="$2" last="$3" cpus="$4" memory="$5" hours="$6"
    ((cpus > MAX_CPUS)) && cpus="$MAX_CPUS"
    ((memory > MAX_MEMORY_GB)) && memory="$MAX_MEMORY_GB"
    ((hours > MAX_TIME_H)) && hours="$MAX_TIME_H"
    local job_name="OR_${SPECIES}_${name}"
    local output_pattern="$OUTDIR/pipeline_info/slurm/%x-%j.out"
    local error_pattern="$OUTDIR/pipeline_info/slurm/%x-%j.err"
    local sbatch_args=(
        --wait
        --parsable
        --export="ALL,NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR:-},NXF_OFFLINE=true"
        --job-name="$job_name"
        --nodes=1
        --ntasks=1
        --cpus-per-task="$cpus"
        --mem="${memory}G"
        --time="${hours}:00:00"
        --output="$output_pattern"
        --error="$error_pattern"
    )
    [[ -n "$PARTITION" ]] && sbatch_args+=(--partition="$PARTITION")
    [[ -n "$ACCOUNT" ]] && sbatch_args+=(--account="$ACCOUNT")
    [[ -n "$QOS" ]] && sbatch_args+=(--qos="$QOS")

    echo
    echo "===== BIOLOGICAL TASK GROUP: $name (stages $first-$last) ====="
    echo "Slurm allocation: ${cpus} CPU, ${memory} GiB, ${hours} h"
    echo "Slurm logs: $OUTDIR/pipeline_info/slurm/${job_name}-<jobid>.{out,err}"
    local bundle_args=(
        "$ROOT_DIR/bin/run_stage_bundle.sh"
        --root "$ROOT_DIR"
        --params "$PARAMS"
        --profile "$BUNDLE_PROFILE"
        --group "$name"
        --from-stage "$first"
        --to-stage "$last"
    )
    [[ -z "$RESUME" ]] && bundle_args+=(--no-resume)
    sbatch "${sbatch_args[@]}" "${bundle_args[@]}"
}

SELECTED_GROUPS=0
SELECTED_STAGES=0
for row in "${TASK_GROUPS[@]}"; do
    IFS='|' read -r name first last mode cpus memory hours <<<"$row"
    ((last < FROM_STAGE || first > TO_STAGE)) && continue
    ((first < FROM_STAGE)) && first="$FROM_STAGE"
    ((last > TO_STAGE)) && last="$TO_STAGE"
    SELECTED_GROUPS=$((SELECTED_GROUPS + 1))
    SELECTED_STAGES=$((SELECTED_STAGES + last - first + 1))

    if [[ "$mode" == "bundle" ]]; then
        submit_bundle "$name" "$first" "$last" "$cpus" "$memory" "$hours"
    else
        for ((stage = first; stage <= last; stage++)); do
            run_compatibility_stage "$stage"
        done
    fi
done

EXPECTED_STAGES=$((TO_STAGE - FROM_STAGE + 1))
((SELECTED_GROUPS > 0)) || {
    echo "ERROR: optimized scheduler selected no task groups" >&2
    exit 70
}
((SELECTED_STAGES == EXPECTED_STAGES)) || {
    echo "ERROR: optimized scheduler covered $SELECTED_STAGES stages; expected $EXPECTED_STAGES" >&2
    exit 70
}

echo
echo "===== FINAL STAGE-GATE AUDIT: $FROM_STAGE-$TO_STAGE ====="
for ((stage = FROM_STAGE; stage <= TO_STAGE; stage++)); do
    python -m or_mining check-stage "$PARAMS" "$stage"
done

echo
echo "PIPELINE COMPLETED THROUGH STAGE $TO_STAGE (optimized scheduler)"
echo "Results: $OUTDIR"
echo "Grouped Slurm logs: $OUTDIR/pipeline_info/slurm"
