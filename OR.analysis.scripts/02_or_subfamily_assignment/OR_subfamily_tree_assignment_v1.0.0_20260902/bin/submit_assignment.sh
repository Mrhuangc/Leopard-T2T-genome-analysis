#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash submit_assignment.sh /abs/config.sh" >&2
    exit 2
fi
CONFIG=$(readlink -f "$1")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$CONFIG"
bash "$SCRIPT_DIR/preflight.sh" "$CONFIG"

JOB_ID=$(sbatch --parsable \
    --partition="${SLURM_PARTITION:-debug}" \
    --job-name=OR_subfamily_assign \
    --cpus-per-task=4 --mem=16G --time=04:00:00 \
    --output="$OUTPUT_ROOT/logs/%x-%j.out" \
    --error="$OUTPUT_ROOT/logs/%x-%j.err" \
    --export=ALL \
    --wrap="bash '$SCRIPT_DIR/run_assignment.sh' '$CONFIG'")
JOB_ID=${JOB_ID%%;*}
printf 'ASSIGNMENT_JOB_SUBMITTED\njob_id=%s\nfinal_marker=%s\n' \
    "$JOB_ID" "$OUTPUT_ROOT/ASSIGNMENT_COMPLETE.json"

