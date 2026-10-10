#!/usr/bin/env python3
"""Validate one family alignment/tree and emit a guarded completion marker."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from Bio import Phylo


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in records:
                    raise ValueError(f"empty or duplicate FASTA ID in {path}: {current!r}")
                records[current] = []
            elif current is None:
                raise ValueError(f"sequence before first header in {path}")
            else:
                records[current].append(line)
    return {key: "".join(value) for key, value in records.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--aligned", required=True, type=Path)
    parser.add_argument("--trimmed", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--iqtree-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mafft-method", required=True)
    args = parser.parse_args()

    for path in (args.input, args.aligned, args.trimmed, args.tree, args.iqtree_report):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"required output missing or empty: {path}")

    input_records = read_fasta(args.input)
    aligned_records = read_fasta(args.aligned)
    trimmed_records = read_fasta(args.trimmed)
    aligned_lengths = {len(seq) for seq in aligned_records.values()}
    trimmed_lengths = {len(seq) for seq in trimmed_records.values()}
    tree = Phylo.read(str(args.tree), "newick")
    tree_ids = [tip.name for tip in tree.get_terminals()]

    checks = {
        "minimum_three_input_records": len(input_records) >= 3,
        "input_and_alignment_ids_match": set(input_records) == set(aligned_records),
        "input_and_trimmed_ids_match": set(input_records) == set(trimmed_records),
        "alignment_is_rectangular": len(aligned_lengths) == 1,
        "trimmed_alignment_is_rectangular": len(trimmed_lengths) == 1,
        "trimmed_alignment_has_at_least_100_sites": bool(trimmed_lengths)
        and next(iter(trimmed_lengths)) >= 100,
        "tree_tip_ids_are_unique": len(tree_ids) == len(set(tree_ids)),
        "tree_and_input_ids_match": set(tree_ids) == set(input_records),
    }
    status = "pass" if all(checks.values()) else "fail"
    payload = {
        "schema_version": "1.0",
        "status": status,
        "family": args.family,
        "method": {
            "alignment": f"MAFFT {args.mafft_method}",
            "trimming": "trimAl automated1",
            "tree": "IQ-TREE model selection with UFBoot and SH-aLRT",
        },
        "counts": {
            "input_records": len(input_records),
            "aligned_records": len(aligned_records),
            "trimmed_records": len(trimmed_records),
            "tree_leaves": len(tree_ids),
            "aligned_sites": next(iter(aligned_lengths)) if len(aligned_lengths) == 1 else None,
            "trimmed_sites": next(iter(trimmed_lengths)) if len(trimmed_lengths) == 1 else None,
        },
        "checks": checks,
        "provenance": {
            "input_path": str(args.input.resolve()),
            "input_sha256": sha256(args.input),
            "aligned_sha256": sha256(args.aligned),
            "trimmed_sha256": sha256(args.trimmed),
            "tree_sha256": sha256(args.tree),
        },
        "outputs": {
            "aligned": str(args.aligned.resolve()),
            "trimmed": str(args.trimmed.resolve()),
            "tree": str(args.tree.resolve()),
            "iqtree_report": str(args.iqtree_report.resolve()),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if status != "pass":
        failed = [key for key, value in checks.items() if not value]
        raise ValueError(f"family-tree validation failed for {args.family}: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

