#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
export PYTHONPATH="$ROOT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"

python -m pytest -q
python -m ruff check \
    bin/build_release_archive.py \
    src/or_mining/conda_prewarm.py \
    src/or_mining/initial_search.py \
    src/or_mining/run_initial_search.py \
    src/or_mining/species_code.py \
    src/or_mining/preflight.py \
    src/or_mining/stage_gate.py \
    src/or_mining/cli.py \
    tests/python/test_species_code.py \
    tests/python/test_conda_prewarm.py \
    tests/python/test_initial_search.py \
    tests/python/test_preflight.py \
    tests/python/test_stage_gate.py

bash -n bin/run_species_pipeline.sh
bash -n bin/run_stage_bundle.sh
bash -n bin/run_initial_search_smoke.sh
bash -n bin/release_self_check.sh
python -m or_mining prewarm-conda "$ROOT_DIR" /tmp/felidae_or_conda_selfcheck --generate-only
bash bin/run_species_pipeline.sh --list-stages > /tmp/felidae_or_task_groups.txt
grep -q '^candidate_locus_discovery' /tmp/felidae_or_task_groups.txt
grep -q '^initial_reference_or_search' /tmp/felidae_or_task_groups.txt
grep -q '^targeted_exonerate_discovery' /tmp/felidae_or_task_groups.txt
grep -q '^or_clade_and_locus_architecture' /tmp/felidae_or_task_groups.txt

nextflow config -profile batch >/dev/null

while IFS= read -r workflow; do
    nextflow lint "$workflow"
done < <(find tests/nextflow -maxdepth 1 -name '*.nf' -type f | sort)

echo "RELEASE SELF-CHECK PASSED"
