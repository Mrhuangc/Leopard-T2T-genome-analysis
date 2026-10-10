#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash /absolute/path/run_candidate_stage.sh /absolute/path/config.sh" >&2
  exit 2
fi

if [[ ! -s "$1" ]]; then
  echo "FATAL: config is missing or empty before path resolution: $1" >&2
  echo "Current directory: $PWD" >&2
  exit 2
fi
CONFIG_PATH=$(readlink -f "$1")
if [[ ! -s "$CONFIG_PATH" ]]; then
  echo "FATAL: missing config: $CONFIG_PATH" >&2
  exit 2
fi
# shellcheck source=/dev/null
source "$CONFIG_PATH"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
for variable in QUERY_FASTA BROAD_CLADE_TSV REFERENCE_FASTA REFERENCE_METADATA OUTPUT_ROOT; do
  value=${!variable:-}
  if [[ -z "$value" ]]; then
    echo "FATAL: config variable is empty: $variable" >&2
    exit 2
  fi
done
for file in "$QUERY_FASTA" "$BROAD_CLADE_TSV" "$REFERENCE_FASTA" "$REFERENCE_METADATA"; do
  if [[ ! -s "$file" ]]; then
    echo "FATAL: required input is missing or empty: $file" >&2
    exit 2
  fi
done
for tool in python diamond; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FATAL: tool not found: $tool" >&2; exit 2; }
done

mkdir -p "$OUTPUT_ROOT/logs" "$OUTPUT_ROOT/search" "$OUTPUT_ROOT/candidate_families"
exec > >(tee "$OUTPUT_ROOT/logs/candidate_stage.stdout.log") \
  2> >(tee "$OUTPUT_ROOT/logs/candidate_stage.stderr.log" >&2)

echo "started_at=$(date --iso-8601=seconds)"
echo "hostname=$(hostname)"
echo "slurm_job_id=${SLURM_JOB_ID:-none}"
python --version
diamond version

QUERY_COUNT=$(grep -c '^>' "$QUERY_FASTA")
REFERENCE_COUNT=$(grep -c '^>' "$REFERENCE_FASTA")
[[ "$QUERY_COUNT" -eq 5338 ]] || {
  echo "FATAL: expected 5338 non-redundant query proteins, observed $QUERY_COUNT" >&2
  exit 3
}
[[ "$REFERENCE_COUNT" -eq 417 ]] || {
  echo "FATAL: expected 417 human reference proteins, observed $REFERENCE_COUNT" >&2
  exit 3
}

DB_PREFIX="$OUTPUT_ROOT/search/human_OR_reference"
DIAMOND_OUT="$OUTPUT_ROOT/search/felid_vs_human_OR.diamond.tsv"
diamond makedb --in "$REFERENCE_FASTA" --db "$DB_PREFIX"
diamond blastp \
  --query "$QUERY_FASTA" \
  --db "${DB_PREFIX}.dmnd" \
  --out "$DIAMOND_OUT" \
  --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore \
  --very-sensitive --max-target-seqs 25 --evalue 1e-10 \
  --threads "${SLURM_CPUS_PER_TASK:-8}"

[[ -s "$DIAMOND_OUT" ]] || { echo "FATAL: DIAMOND output is empty" >&2; exit 4; }
python "$SCRIPT_DIR/prepare_family_candidates.py" \
  --query-fasta "$QUERY_FASTA" \
  --reference-fasta "$REFERENCE_FASTA" \
  --reference-metadata "$REFERENCE_METADATA" \
  --broad-clade-tsv "$BROAD_CLADE_TSV" \
  --diamond-tsv "$DIAMOND_OUT" \
  --output-dir "$OUTPUT_ROOT/candidate_families"

python - "$OUTPUT_ROOT/candidate_families/PREPARATION_COMPLETE.json" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
if payload.get("status") != "pass":
    raise SystemExit(f"FATAL: completion marker is not pass: {path}")
print(json.dumps(payload["counts"], indent=2))
PY

echo "finished_at=$(date --iso-8601=seconds)"
echo "CANDIDATE_STAGE_COMPLETE=$OUTPUT_ROOT/candidate_families/PREPARATION_COMPLETE.json"
