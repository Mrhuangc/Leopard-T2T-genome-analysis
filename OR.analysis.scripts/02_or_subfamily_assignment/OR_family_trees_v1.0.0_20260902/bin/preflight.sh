#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash preflight.sh /abs/config.sh" >&2
    exit 2
fi
CONFIG=$(readlink -f "$1")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
# shellcheck source=/dev/null
source "$CONFIG"

for name in MANIFEST BRIDGE_ROOT OUTPUT_ROOT; do
    [[ -n ${!name:-} ]] || { echo "FATAL: $name is unset" >&2; exit 1; }
done
for path in "$MANIFEST"; do
    [[ -s "$path" ]] || { echo "FATAL: required file missing or empty: $path" >&2; exit 1; }
done
[[ -d "$BRIDGE_ROOT" ]] || { echo "FATAL: BRIDGE_ROOT missing: $BRIDGE_ROOT" >&2; exit 1; }

PYTHON_BIN=${PYTHON_BIN:-python}
MAFFT_BIN=${MAFFT_BIN:-mafft}
TRIMAL_BIN=${TRIMAL_BIN:-trimal}
IQTREE_BIN=${IQTREE_BIN:-iqtree3}
for tool in "$PYTHON_BIN" "$MAFFT_BIN" "$TRIMAL_BIN" "$IQTREE_BIN" sbatch; do
    command -v "$tool" >/dev/null || { echo "FATAL: tool not found: $tool" >&2; exit 1; }
done
"$PYTHON_BIN" -c 'import Bio'
"$PYTHON_BIN" - "${PREPARATION_MARKER:-}" "${DOG_BRIDGE_MARKER:-}" "${SMOKE_MARKER:-}" <<'PY'
import json, pathlib, sys
for raw in sys.argv[1:]:
    if not raw:
        continue
    path = pathlib.Path(raw)
    payload = json.loads(path.read_text())
    if payload.get("status") != "pass":
        raise SystemExit(f"FATAL: upstream marker is not pass: {path}")
PY

for list in "$PACKAGE_ROOT/lists/large.txt" "$PACKAGE_ROOT/lists/medium.txt" "$PACKAGE_ROOT/lists/small.txt"; do
    while IFS= read -r family; do
        [[ -n "$family" ]] || continue
        "$PYTHON_BIN" "$SCRIPT_DIR/resolve_family_input.py" \
            --manifest "$MANIFEST" --bridge-root "$BRIDGE_ROOT" --family "$family" >/dev/null
    done < "$list"
done

mkdir -p "$OUTPUT_ROOT/logs" "$OUTPUT_ROOT/state"
echo "PREFLIGHT_PASS families=17 output_root=$OUTPUT_ROOT"

