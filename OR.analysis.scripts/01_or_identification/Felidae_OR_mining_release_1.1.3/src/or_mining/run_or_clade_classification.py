"""CLI for OR evolutionary-clade assignment."""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from or_mining.or_clade_classification import classify_or_clades


def main() -> int:
    parser = ArgumentParser()
    parser.add_argument("tree", type=Path)
    parser.add_argument("template_alignment", type=Path)
    parser.add_argument("type1_ids", type=Path)
    parser.add_argument("type2_ids", type=Path)
    parser.add_argument("catalog", type=Path)
    parser.add_argument("functional_protein", type=Path)
    parser.add_argument("functional_blastp", type=Path)
    parser.add_argument("pseudogene_blastx", type=Path)
    parser.add_argument("truncated_blastx", type=Path)
    parser.add_argument("ambiguous_blastx", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--minimum-tree-support", type=float, default=0.80)
    args = parser.parse_args()
    classify_or_clades(
        args.species, args.tree, args.template_alignment, args.type1_ids, args.type2_ids,
        args.catalog, args.functional_protein, args.functional_blastp,
        args.pseudogene_blastx, args.truncated_blastx, args.ambiguous_blastx,
        args.output_prefix, args.minimum_tree_support,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
