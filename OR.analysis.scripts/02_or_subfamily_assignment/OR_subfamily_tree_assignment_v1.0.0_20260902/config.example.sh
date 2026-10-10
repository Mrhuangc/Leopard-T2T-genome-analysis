#!/usr/bin/env bash
OR_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
DATA_ROOT="${OR_DATA_ROOT:-/ABSOLUTE/PATH/OR_data}"
RESULT_ROOT="${OR_RESULT_ROOT:-$OR_ROOT/results}"
FAMILY_TREE_ROOT="$RESULT_ROOT/OR_family_trees_v1"
CANDIDATE_ASSIGNMENTS="$RESULT_ROOT/subfamily_candidates/candidate_families/candidate_assignments.tsv"
HUMAN_REFERENCE_METADATA="$DATA_ROOT/reference/human_OR_reference.metadata.tsv"
IDENTICAL_SEQUENCE_MAP="$DATA_ROOT/joint/functional_OR_identical_sequence_map.tsv"
SPECIES_METADATA="$OR_ROOT/config/species_metadata.tsv"
OUTPUT_ROOT="$RESULT_ROOT/OR_subfamily_tree_assignment_v1"
PYTHON_BIN="python"
SLURM_PARTITION="debug"
HIGH_SH_ALRT=80
HIGH_UFBOOT=95
HIGH_DISTANCE_RATIO=1.10
MEDIUM_DISTANCE_RATIO=1.02
