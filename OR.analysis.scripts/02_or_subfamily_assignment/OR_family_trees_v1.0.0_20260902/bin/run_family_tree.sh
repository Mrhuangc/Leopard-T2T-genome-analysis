#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: bash run_family_tree.sh /abs/config.sh OR_FAMILY" >&2
    exit 2
fi

CONFIG=$(readlink -f "$1")
FAMILY=$2
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

# shellcheck source=/dev/null
source "$CONFIG"

: "${MANIFEST:?MANIFEST is required}"
: "${BRIDGE_ROOT:?BRIDGE_ROOT is required}"
: "${OUTPUT_ROOT:?OUTPUT_ROOT is required}"

PYTHON_BIN=${PYTHON_BIN:-python}
MAFFT_BIN=${MAFFT_BIN:-mafft}
TRIMAL_BIN=${TRIMAL_BIN:-trimal}
IQTREE_BIN=${IQTREE_BIN:-iqtree3}
IQTREE_MODEL=${IQTREE_MODEL:-MFP}
IQTREE_BOOTSTRAP=${IQTREE_BOOTSTRAP:-1000}
IQTREE_ALRT=${IQTREE_ALRT:-1000}
IQTREE_SEED=${IQTREE_SEED:-20260902}
CPUS=${SLURM_CPUS_PER_TASK:-8}

INPUT=$($PYTHON_BIN "$SCRIPT_DIR/resolve_family_input.py" \
    --manifest "$MANIFEST" --bridge-root "$BRIDGE_ROOT" --family "$FAMILY")
INPUT=$(readlink -f "$INPUT")
FAMILY_ROOT="$OUTPUT_ROOT/families/$FAMILY"
ALIGN_DIR="$FAMILY_ROOT/alignment"
TREE_DIR="$FAMILY_ROOT/tree"
LOG_DIR="$FAMILY_ROOT/logs"
MARKER="$FAMILY_ROOT/COMPLETE.json"
mkdir -p "$ALIGN_DIR" "$TREE_DIR" "$LOG_DIR"

if [[ -s "$MARKER" ]]; then
    if "$PYTHON_BIN" - "$MARKER" "$INPUT" <<'PY'
import hashlib, json, pathlib, sys
marker = json.loads(pathlib.Path(sys.argv[1]).read_text())
digest = hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest()
raise SystemExit(0 if marker.get("status") == "pass" and marker.get("provenance", {}).get("input_sha256") == digest else 1)
PY
    then
        echo "CHECKPOINT_REUSED family=$FAMILY marker=$MARKER"
        exit 0
    fi
fi

cp -f "$INPUT" "$FAMILY_ROOT/input.selected.prot.fa"
INPUT_LOCAL="$FAMILY_ROOT/input.selected.prot.fa"
NSEQ=$(grep -c '^>' "$INPUT_LOCAL")
if (( NSEQ >= 500 )); then
    MAFFT_ARGS=(--retree 2 --maxiterate 0)
    MAFFT_METHOD="FFT-NS-2"
elif (( NSEQ >= 200 )); then
    MAFFT_ARGS=(--retree 2 --maxiterate 100)
    MAFFT_METHOD="FFT-NS-i-100"
else
    MAFFT_ARGS=(--localpair --maxiterate 1000)
    MAFFT_METHOD="L-INS-i"
fi

ALIGNED="$ALIGN_DIR/${FAMILY}.mafft.fa"
TRIMMED="$ALIGN_DIR/${FAMILY}.trimmed.fa"
PREFIX="$TREE_DIR/${FAMILY}"

"$MAFFT_BIN" --thread "$CPUS" "${MAFFT_ARGS[@]}" "$INPUT_LOCAL" \
    > "$ALIGNED" 2> "$LOG_DIR/mafft.log"
"$TRIMAL_BIN" -in "$ALIGNED" -out "$TRIMMED" -automated1 \
    > "$LOG_DIR/trimal.stdout.log" 2> "$LOG_DIR/trimal.stderr.log"

"$IQTREE_BIN" -s "$TRIMMED" -m "$IQTREE_MODEL" \
    -B "$IQTREE_BOOTSTRAP" --alrt "$IQTREE_ALRT" \
    -T "$CPUS" --seed "$IQTREE_SEED" --redo --prefix "$PREFIX" \
    > "$LOG_DIR/iqtree.stdout.log" 2> "$LOG_DIR/iqtree.stderr.log"

"$PYTHON_BIN" "$SCRIPT_DIR/validate_family_tree.py" \
    --family "$FAMILY" \
    --input "$INPUT_LOCAL" \
    --aligned "$ALIGNED" \
    --trimmed "$TRIMMED" \
    --tree "${PREFIX}.treefile" \
    --iqtree-report "${PREFIX}.iqtree" \
    --mafft-method "$MAFFT_METHOD" \
    --output "$MARKER"

echo "FAMILY_TREE_COMPLETED family=$FAMILY marker=$MARKER"

