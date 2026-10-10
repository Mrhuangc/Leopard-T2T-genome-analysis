#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"
RUN_SCRIPT="$PACKAGE_ROOT/bin/run_macro_synteny.sh"
OUT_ROOT="${OUT_ROOT:?Set OUT_ROOT to a writable results directory before submission}"
mkdir -p "$OUT_ROOT/logs"

SAMTOOLS="$(command -v samtools || true)"
MINIMAP2="$(command -v minimap2 || true)"
PYTHON="$(command -v python || true)"
[[ -x "$SAMTOOLS" ]] || { printf 'FATAL: samtools not found\n' >&2; exit 1; }
[[ -x "$MINIMAP2" ]] || { printf 'FATAL: minimap2 not found\n' >&2; exit 1; }
[[ -x "$PYTHON" ]] || { printf 'FATAL: python not found\n' >&2; exit 1; }
"$PYTHON" "$PACKAGE_ROOT/bin/validate_package.py"

job_id="$(sbatch --parsable \
  --job-name=OR_ChrD1_macro \
  --partition=debug \
  --cpus-per-task=32 \
  --mem=64G \
  --time=06:00:00 \
  --output="$OUT_ROOT/logs/slurm-%j.out" \
  --error="$OUT_ROOT/logs/slurm-%j.err" \
  --export=ALL,PACKAGE_ROOT="$PACKAGE_ROOT",OUT_ROOT="$OUT_ROOT",SAMTOOLS="$SAMTOOLS",MINIMAP2="$MINIMAP2",PYTHON="$PYTHON" \
  "$RUN_SCRIPT")"
printf 'job_id=%s\noutput_root=%s\n' "$job_id" "$OUT_ROOT"
