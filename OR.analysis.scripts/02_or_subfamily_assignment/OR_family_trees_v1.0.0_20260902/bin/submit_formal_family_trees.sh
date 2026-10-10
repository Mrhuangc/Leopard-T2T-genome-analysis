#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash submit_formal_family_trees.sh /abs/config.sh" >&2
    exit 2
fi
CONFIG=$(readlink -f "$1")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
# shellcheck source=/dev/null
source "$CONFIG"

bash "$SCRIPT_DIR/preflight.sh" "$CONFIG"

SLURM_PARTITION=${SLURM_PARTITION:-debug}
MAX_PARALLEL_LARGE=${MAX_PARALLEL_LARGE:-4}
MAX_PARALLEL_MEDIUM=${MAX_PARALLEL_MEDIUM:-4}
MAX_PARALLEL_SMALL=${MAX_PARALLEL_SMALL:-4}
LOG_DIR="$OUTPUT_ROOT/logs"
mkdir -p "$LOG_DIR"

submit_array() {
    local tier=$1 list=$2 cpus=$3 memory=$4 time_limit=$5 parallel=$6
    local count
    count=$(grep -Evc '^\s*(#|$)' "$list")
    local last=$((count - 1))
    sbatch --parsable \
        --partition="$SLURM_PARTITION" \
        --job-name="ORtree_${tier}" \
        --array="0-${last}%${parallel}" \
        --cpus-per-task="$cpus" \
        --mem="$memory" \
        --time="$time_limit" \
        --output="$LOG_DIR/%x-%A_%a.out" \
        --error="$LOG_DIR/%x-%A_%a.err" \
        --export=ALL \
        --wrap="bash '$SCRIPT_DIR/run_family_from_list.sh' '$CONFIG' '$list'"
}

LARGE_ID=$(submit_array large "$PACKAGE_ROOT/lists/large.txt" 16 64G 36:00:00 "$MAX_PARALLEL_LARGE")
MEDIUM_ID=$(submit_array medium "$PACKAGE_ROOT/lists/medium.txt" 12 32G 24:00:00 "$MAX_PARALLEL_MEDIUM")
SMALL_ID=$(submit_array small "$PACKAGE_ROOT/lists/small.txt" 8 16G 12:00:00 "$MAX_PARALLEL_SMALL")

# --parsable can append a cluster name after a semicolon.
LARGE_ID=${LARGE_ID%%;*}
MEDIUM_ID=${MEDIUM_ID%%;*}
SMALL_ID=${SMALL_ID%%;*}

FINAL_ID=$(sbatch --parsable \
    --partition="$SLURM_PARTITION" \
    --job-name="ORtree_finalize" \
    --dependency="afterany:${LARGE_ID}:${MEDIUM_ID}:${SMALL_ID}" \
    --cpus-per-task=1 --mem=4G --time=01:00:00 \
    --output="$LOG_DIR/%x-%j.out" \
    --error="$LOG_DIR/%x-%j.err" \
    --export=ALL \
    --wrap="bash '$SCRIPT_DIR/finalize.sh' '$CONFIG'")
FINAL_ID=${FINAL_ID%%;*}

cat <<EOF
FORMAL_FAMILY_TREE_JOBS_SUBMITTED
large_array=$LARGE_ID
medium_array=$MEDIUM_ID
small_array=$SMALL_ID
final_validator=$FINAL_ID
maximum_simultaneous_tree_tasks=$((MAX_PARALLEL_LARGE + MAX_PARALLEL_MEDIUM + MAX_PARALLEL_SMALL))
final_marker=$OUTPUT_ROOT/state/FORMAL_FAMILY_TREES_COMPLETE.json
EOF

