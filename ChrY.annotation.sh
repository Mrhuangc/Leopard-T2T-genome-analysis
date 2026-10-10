#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YTOOLS="${REPO_ROOT}/ChrY.analysis.scripts"

#protein coding gene annotation
python3 "${YTOOLS}/annotation.scripts/run_annotation.py" \
  --config "${YTOOLS}/annotation.scripts/config/FCA.json"
python3 "${YTOOLS}/annotation.scripts/run_annotation.py" \
  --config "${YTOOLS}/annotation.scripts/config/PPA.json"

#pseudogenes annotation
# Install mappy in your environment before running (see envs/requirements.txt).
# Copy/edit the supplied example configs first; do not overwrite a configured file.

python3 "${YTOOLS}/annotation.scripts/annotate_pseudogenes.py" \
  --config "${YTOOLS}/annotation.scripts/config/config.json" \
  --out path/results \
  --miniprot miniprot --threads 8
