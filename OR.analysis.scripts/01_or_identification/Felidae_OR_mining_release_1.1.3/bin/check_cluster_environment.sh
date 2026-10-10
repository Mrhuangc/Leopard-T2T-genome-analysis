#!/usr/bin/env bash
set -u

status=0
environment_engine=""

check_command() {
    local name="$1"
    if path="$(command -v "$name" 2>/dev/null)"; then
        printf '[OK] %-14s %s\n' "$name" "$path"
    else
        printf '[MISSING] %s\n' "$name"
        status=1
    fi
}

report_optional_command() {
    local name="$1"
    if path="$(command -v "$name" 2>/dev/null)"; then
        printf '[AVAILABLE] %-9s %s\n' "$name" "$path"
    else
        printf '[DEFERRED] %s (provided by a per-process environment)\n' "$name"
    fi
}

echo "== Host =="
hostname || true
uname -a || true

echo
echo "== Required control commands =="
for cmd in bash java nextflow python sbatch squeue; do
    check_command "$cmd"
done

echo
echo "== Environment engines =="
if command -v apptainer >/dev/null 2>&1; then
    check_command apptainer
    environment_engine="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    check_command singularity
    environment_engine="singularity"
elif command -v mamba >/dev/null 2>&1; then
    check_command mamba
    environment_engine="mamba"
elif command -v conda >/dev/null 2>&1; then
    check_command conda
    environment_engine="conda"
else
    echo "[MISSING] apptainer/singularity/mamba/conda"
    status=1
fi

echo
echo "== Scientific tools (optional in the control environment) =="
for cmd in tblastn blastp blastx samtools getorf transeq miniprot \
           genewise exonerate mafft iqtree cd-hit; do
    report_optional_command "$cmd"
done

echo
echo "== Versions =="
java -version 2>&1 | head -3 || true
nextflow -version 2>&1 || true
python --version 2>&1 || true
sbatch --version 2>&1 || true
case "$environment_engine" in
    apptainer) apptainer --version 2>&1 || true ;;
    singularity) singularity --version 2>&1 || true ;;
    mamba) mamba --version 2>&1 || true ;;
    conda) conda --version 2>&1 || true ;;
esac

echo
if (( status == 0 )); then
    echo "ENVIRONMENT CHECK PASSED"
else
    echo "ENVIRONMENT CHECK INCOMPLETE"
fi

exit "$status"
