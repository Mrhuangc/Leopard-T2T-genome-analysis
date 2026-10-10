#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash preflight.sh /abs/config.sh" >&2
    exit 2
fi
CONFIG=$(readlink -f "$1")
# shellcheck source=/dev/null
source "$CONFIG"

for name in FAMILY_TREE_ROOT CANDIDATE_ASSIGNMENTS HUMAN_REFERENCE_METADATA \
            IDENTICAL_SEQUENCE_MAP SPECIES_METADATA OUTPUT_ROOT; do
    [[ -n ${!name:-} ]] || { echo "FATAL: $name is unset" >&2; exit 1; }
done
for path in "$CANDIDATE_ASSIGNMENTS" "$HUMAN_REFERENCE_METADATA" \
            "$IDENTICAL_SEQUENCE_MAP" "$SPECIES_METADATA" \
            "$FAMILY_TREE_ROOT/state/FORMAL_FAMILY_TREES_COMPLETE.json"; do
    [[ -s "$path" ]] || { echo "FATAL: required input missing or empty: $path" >&2; exit 1; }
done

PYTHON_BIN=${PYTHON_BIN:-python}
command -v "$PYTHON_BIN" >/dev/null || { echo "FATAL: Python not found: $PYTHON_BIN" >&2; exit 1; }
command -v sbatch >/dev/null || { echo "FATAL: sbatch not found" >&2; exit 1; }

"$PYTHON_BIN" - "$FAMILY_TREE_ROOT/state/FORMAL_FAMILY_TREES_COMPLETE.json" <<'PY'
import json, pathlib, sys
p=pathlib.Path(sys.argv[1])
d=json.loads(p.read_text())
checks=d.get("checks", {})
if d.get("status") != "pass" or d.get("counts", {}).get("completed_families") != 17 or not all(checks.values()):
    raise SystemExit(f"FATAL: formal family-tree marker is not a complete pass: {p}")
PY

for family in OR1 OR2 OR3 OR4 OR5 OR6 OR7 OR8 OR9 OR10 OR11 OR12 OR13 OR14 OR51 OR52 OR56; do
    for path in \
        "$FAMILY_TREE_ROOT/families/$family/input.selected.prot.fa" \
        "$FAMILY_TREE_ROOT/families/$family/tree/$family.treefile" \
        "$FAMILY_TREE_ROOT/families/$family/COMPLETE.json"; do
        [[ -s "$path" ]] || { echo "FATAL: family-tree input missing: $path" >&2; exit 1; }
    done
done

mkdir -p "$OUTPUT_ROOT" "$OUTPUT_ROOT/logs"
echo "PREFLIGHT_PASS family_trees=17 output_root=$OUTPUT_ROOT"

