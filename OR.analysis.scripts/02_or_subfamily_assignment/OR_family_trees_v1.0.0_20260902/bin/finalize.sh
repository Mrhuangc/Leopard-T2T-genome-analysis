#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash finalize.sh /abs/config.sh" >&2
    exit 2
fi
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$(readlink -f "$1")"
PYTHON_BIN=${PYTHON_BIN:-python}
exec "$PYTHON_BIN" "$SCRIPT_DIR/validate_all_family_trees.py" --output-root "$OUTPUT_ROOT"

