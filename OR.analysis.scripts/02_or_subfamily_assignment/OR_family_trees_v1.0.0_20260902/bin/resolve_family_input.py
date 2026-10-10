#!/usr/bin/env python3
"""Resolve and validate one family FASTA from a bridge manifest."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".faa", ".fasta", ".fas", ".prot"}
FAMILY_COLUMNS = ("family", "or_family", "family_id", "candidate_family")
COUNT_COLUMNS = ("total_sequences", "total_records", "sequence_count", "records")


def fasta_count(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def normalized(field: str) -> str:
    return field.strip().lower().replace(" ", "_")


def family_token_match(name: str, family: str) -> bool:
    return bool(re.search(rf"(^|[._-]){re.escape(family)}([._-]|$)", name, re.I))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--bridge-root", required=True, type=Path)
    parser.add_argument("--family", required=True)
    args = parser.parse_args()

    manifest = args.manifest.resolve()
    bridge_root = args.bridge_root.resolve()
    if not manifest.is_file() or manifest.stat().st_size == 0:
        raise SystemExit(f"manifest missing or empty: {manifest}")
    if not bridge_root.is_dir():
        raise SystemExit(f"bridge root missing: {bridge_root}")

    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("manifest has no header")
        original = {normalized(x): x for x in reader.fieldnames}
        family_key = next((original[x] for x in FAMILY_COLUMNS if x in original), None)
        if family_key is None:
            raise SystemExit(f"cannot identify family column: {reader.fieldnames}")
        rows = [row for row in reader if row[family_key].strip().upper() == args.family.upper()]

    if len(rows) != 1:
        raise SystemExit(f"expected one manifest row for {args.family}, observed {len(rows)}")
    row = rows[0]

    candidates: set[Path] = set()
    for value in row.values():
        value = (value or "").strip()
        if not value:
            continue
        path = Path(value)
        if path.suffix.lower() not in FASTA_SUFFIXES:
            continue
        probes = [path] if path.is_absolute() else [manifest.parent / path, bridge_root / path]
        for probe in probes:
            if probe.is_file() and probe.stat().st_size > 0:
                candidates.add(probe.resolve())
        if not path.is_absolute():
            for probe in bridge_root.rglob(path.name):
                if probe.is_file() and probe.stat().st_size > 0:
                    candidates.add(probe.resolve())

    if not candidates:
        for probe in bridge_root.rglob("*"):
            if (
                probe.is_file()
                and probe.suffix.lower() in FASTA_SUFFIXES
                and probe.stat().st_size > 0
                and family_token_match(probe.name, args.family)
            ):
                candidates.add(probe.resolve())

    def score(path: Path) -> tuple[int, int, str]:
        name = path.name.lower()
        preferred = sum(token in name for token in ("with_dog", "merged", "family"))
        return (preferred, path.stat().st_size, str(path))

    if not candidates:
        raise SystemExit(f"no FASTA resolved for {args.family}")
    ranked = sorted(candidates, key=score, reverse=True)
    best = ranked[0]
    if len(ranked) > 1 and score(ranked[0])[:2] == score(ranked[1])[:2]:
        raise SystemExit(f"ambiguous FASTA resolution for {args.family}: {ranked[:5]}")

    observed = fasta_count(best)
    if observed == 0:
        raise SystemExit(f"resolved FASTA is empty: {best}")

    normalized_row = {normalized(key): (value or "").strip() for key, value in row.items()}
    expected_text = next((normalized_row[x] for x in COUNT_COLUMNS if normalized_row.get(x)), "")
    if expected_text:
        expected = int(expected_text)
        if observed != expected:
            raise SystemExit(
                f"record-count mismatch for {args.family}: manifest={expected}, fasta={observed}, {best}"
            )

    print(best)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

