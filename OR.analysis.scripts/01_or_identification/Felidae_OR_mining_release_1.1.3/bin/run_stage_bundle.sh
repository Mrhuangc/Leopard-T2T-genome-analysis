#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=""
PARAMS=""
PROFILE="batch,conda"
GROUP="unnamed_group"
FROM_STAGE=""
TO_STAGE=""
RESUME="-resume"

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

while (($#)); do
    case "$1" in
        --root) ROOT_DIR="${2:-}"; shift 2 ;;
        --params) PARAMS="${2:-}"; shift 2 ;;
        --profile) PROFILE="${2:-}"; shift 2 ;;
        --group) GROUP="${2:-}"; shift 2 ;;
        --from-stage) FROM_STAGE="${2:-}"; shift 2 ;;
        --to-stage) TO_STAGE="${2:-}"; shift 2 ;;
        --no-resume) RESUME=""; shift ;;
        *) echo "ERROR: unknown bundle option: $1" >&2; exit 2 ;;
    esac
done

[[ -d "$ROOT_DIR" ]] || { echo "ERROR: invalid --root: $ROOT_DIR" >&2; exit 2; }
[[ -f "$PARAMS" ]] || { echo "ERROR: invalid --params: $PARAMS" >&2; exit 2; }
[[ "$FROM_STAGE" =~ ^[0-9]+$ && "$TO_STAGE" =~ ^[0-9]+$ ]] || {
    echo "ERROR: bundle stage range must be numeric" >&2; exit 2;
}

cd "$ROOT_DIR"
export PYTHONPATH="$ROOT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
echo "GROUP=$GROUP"
echo "HOST=$(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-unset}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
echo "PROFILE=$PROFILE"

# Environments were materialized by the login-node controller. Avoid network and
# Nextflow version checks inside the compute allocation.
export NXF_OFFLINE=true
if [[ ",$PROFILE," == *,conda,* ]]; then
    [[ -n "${NXF_CONDA_CACHEDIR:-}" && -d "$NXF_CONDA_CACHEDIR" ]] || {
        echo "ERROR: shared NXF_CONDA_CACHEDIR was not prepared before submission" >&2
        exit 2
    }
    echo "NXF_CONDA_CACHEDIR=$NXF_CONDA_CACHEDIR"
fi

for ((stage = FROM_STAGE; stage <= TO_STAGE; stage++)); do
    script="tests/nextflow/${STAGES[$((stage - 1))]}"
    echo
    echo "----- $GROUP: compatibility stage $(printf '%02d' "$stage") -----"
    command=(nextflow run "$script" -profile "$PROFILE" -params-file "$PARAMS")
    [[ -n "$RESUME" ]] && command+=("$RESUME")
    "${command[@]}"
    python -m or_mining check-stage "$PARAMS" "$stage"
done

echo "GROUP_COMPLETED=$GROUP"
