"""CLI for strict multi-exon OR sequence reconstruction and classification."""

from __future__ import annotations

import argparse
from pathlib import Path

from or_mining.multiexon_classification import classify_multiexon_models


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("strict_raw", type=Path)
    parser.add_argument("strict_manifest", type=Path)
    parser.add_argument("regions_tsv", type=Path)
    parser.add_argument("regions_fasta", type=Path)
    parser.add_argument("reference_or", type=Path)
    parser.add_argument("species_or", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--minimum-orf-nt", type=int, default=750)
    parser.add_argument("--minimum-model-overlap", type=float, default=0.80)
    parser.add_argument("--maximum-exons", type=int, default=2)
    parser.add_argument("--edge-distance-bp", type=int, default=100)
    parser.add_argument("--edge-n-window-bp", type=int, default=200)
    parser.add_argument("--edge-n-run-bp", type=int, default=50)
    args = parser.parse_args()
    classify_multiexon_models(
        args.strict_raw,
        args.strict_manifest,
        args.regions_tsv,
        args.regions_fasta,
        args.reference_or,
        args.species_or,
        args.output_prefix,
        args.minimum_orf_nt,
        args.minimum_model_overlap,
        args.maximum_exons,
        args.edge_distance_bp,
        args.edge_n_window_bp,
        args.edge_n_run_bp,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
