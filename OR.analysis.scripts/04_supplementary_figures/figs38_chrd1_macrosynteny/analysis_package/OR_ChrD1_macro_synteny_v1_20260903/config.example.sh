#!/usr/bin/env bash
PACKAGE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
export REGIONS_TSV="$PACKAGE_ROOT/input/regions.tsv"
export ARCHITECTURE_DIR=/ABSOLUTE/PATH/OR_data/architecture
export OUT_ROOT=/ABSOLUTE/PATH/OR_results/macrosynteny
export THREADS=32
