#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
export PYTHONPATH="$ROOT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"

nextflow run ./tests/nextflow/initial_tblastn.nf \
    -profile standard,conda \
    -params-file tests/data/initial_search/params.yml \
    -resume

python -m or_mining check-stage tests/data/initial_search/params.yml 0
echo "INITIAL SEARCH FUNCTIONAL SMOKE PASSED"
