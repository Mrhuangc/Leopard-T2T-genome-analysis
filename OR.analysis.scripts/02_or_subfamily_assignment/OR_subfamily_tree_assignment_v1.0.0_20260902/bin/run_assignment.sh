#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash run_assignment.sh /abs/config.sh" >&2
    exit 2
fi
CONFIG=$(readlink -f "$1")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$CONFIG"
PYTHON_BIN=${PYTHON_BIN:-python}

exec "$PYTHON_BIN" "$SCRIPT_DIR/assign_subfamilies.py" \
    --family-tree-root "$FAMILY_TREE_ROOT" \
    --candidate-assignments "$CANDIDATE_ASSIGNMENTS" \
    --human-metadata "$HUMAN_REFERENCE_METADATA" \
    --identical-sequence-map "$IDENTICAL_SEQUENCE_MAP" \
    --species-metadata "$SPECIES_METADATA" \
    --output-root "$OUTPUT_ROOT" \
    --high-sh-alrt "${HIGH_SH_ALRT:-80}" \
    --high-ufboot "${HIGH_UFBOOT:-95}" \
    --high-distance-ratio "${HIGH_DISTANCE_RATIO:-1.10}" \
    --medium-distance-ratio "${MEDIUM_DISTANCE_RATIO:-1.02}"

