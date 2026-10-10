#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: bash run_family_from_list.sh /abs/config.sh /abs/family_list.txt" >&2
    exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
INDEX=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
mapfile -t FAMILIES < <(grep -Ev '^\s*(#|$)' "$2")
if (( INDEX < 0 || INDEX >= ${#FAMILIES[@]} )); then
    echo "Invalid array index $INDEX for $2" >&2
    exit 2
fi
exec bash "$SCRIPT_DIR/run_family_tree.sh" "$1" "${FAMILIES[$INDEX]}"

