#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash submit_candidate_stage.sh /absolute/path/config.sh" >&2
  exit 2
fi
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
RUNNER="$SCRIPT_DIR/run_candidate_stage.sh"
if [[ ! -s "$RUNNER" ]]; then
  echo "FATAL: runner is missing or empty: $RUNNER" >&2
  exit 2
fi
if [[ ! -s "$1" ]]; then
  echo "FATAL: config is missing or empty: $1" >&2
  echo "Current directory: $PWD" >&2
  exit 2
fi
CONFIG=$(readlink -f "$1")
# shellcheck source=/dev/null
source "$CONFIG"
for variable in QUERY_FASTA BROAD_CLADE_TSV REFERENCE_FASTA REFERENCE_METADATA OUTPUT_ROOT; do
  value=${!variable:-}
  if [[ -z "$value" ]]; then
    echo "FATAL: config variable is empty: $variable" >&2
    exit 2
  fi
  if [[ "$value" == *'/USER/'* ]]; then
    echo "FATAL: replace USER in config variable $variable: $value" >&2
    exit 2
  fi
done
for file in "$QUERY_FASTA" "$BROAD_CLADE_TSV" "$REFERENCE_FASTA" "$REFERENCE_METADATA"; do
  if [[ ! -s "$file" ]]; then
    echo "FATAL: required input is missing or empty: $file" >&2
    exit 2
  fi
done
command -v sbatch >/dev/null 2>&1 || {
  echo "FATAL: sbatch is not available in PATH" >&2
  exit 2
}
mkdir -p "$OUTPUT_ROOT/slurm"

echo "SUBMISSION_PREFLIGHT_PASS"
echo "config=$CONFIG"
echo "runner=$RUNNER"
echo "output_root=$OUTPUT_ROOT"

JOB_ID=$(sbatch --parsable \
  --job-name=OR_subfamily_candidates \
  --partition=debug \
  --cpus-per-task=8 \
  --mem=16G \
  --time=02:00:00 \
  --output="$OUTPUT_ROOT/slurm/OR_subfamily_candidates-%j.out" \
  --error="$OUTPUT_ROOT/slurm/OR_subfamily_candidates-%j.err" \
  --wrap="bash '$RUNNER' '$CONFIG'")

echo "SUBMITTED job_id=$JOB_ID"
echo "stdout=$OUTPUT_ROOT/slurm/OR_subfamily_candidates-${JOB_ID}.out"
echo "stderr=$OUTPUT_ROOT/slurm/OR_subfamily_candidates-${JOB_ID}.err"
