#!/usr/bin/env python3
"""Validate package indexes, entrypoint paths, and root SHA-256 records."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    problems: list[str] = []
    manifest_path = ROOT / "scripts_manifest.tsv"
    rows = read_tsv(manifest_path)
    required = {
        "module_id",
        "entrypoint_or_package",
        "role",
        "required_inputs",
        "outputs_or_result",
        "dependencies_and_scope",
    }
    if not rows or not required.issubset(rows[0]):
        problems.append("scripts_manifest.tsv is empty or has an invalid header")
    module_ids = [row.get("module_id", "") for row in rows]
    if len(module_ids) != len(set(module_ids)):
        problems.append("scripts_manifest.tsv contains duplicate module_id values")
    for row in rows:
        entry = ROOT / row.get("entrypoint_or_package", "")
        if not entry.exists():
            problems.append(f"missing manifest entry: {row.get('module_id')} -> {entry}")

    source_rows = read_tsv(ROOT / "source_map.tsv")
    if not source_rows or len(source_rows[0]) != 5:
        problems.append("source_map.tsv is empty or does not have five tab-separated columns")

    checksum_path = ROOT / "SHA256SUMS.txt"
    indexed: set[str] = set()
    if not checksum_path.is_file():
        problems.append("SHA256SUMS.txt is missing")
    else:
        for line_number, line in enumerate(
            checksum_path.read_text(encoding="utf-8-sig").splitlines(), start=1
        ):
            if not line.strip():
                continue
            try:
                expected, relative = line.split(maxsplit=1)
            except ValueError:
                problems.append(f"malformed checksum row {line_number}")
                continue
            relative = relative.strip().lstrip("*").replace("\\", "/")
            indexed.add(relative)
            target = ROOT / Path(relative)
            if not target.is_file():
                problems.append(f"checksum target is missing: {relative}")
                continue
            observed = hashlib.sha256(target.read_bytes()).hexdigest()
            if observed.lower() != expected.lower():
                problems.append(f"checksum mismatch: {relative}")

    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path != checksum_path
    }
    if indexed != actual:
        missing = sorted(actual - indexed)
        stale = sorted(indexed - actual)
        if missing:
            problems.append(f"files absent from root checksum index: {missing[:5]}")
        if stale:
            problems.append(f"stale root checksum entries: {stale[:5]}")

    if problems:
        for problem in problems:
            print(f"FAIL: {problem}")
        return 1
    print(
        f"PACKAGE_INDEX_PASS: {len(rows)} manifest entries, "
        f"{len(source_rows)} source-map rows, {len(indexed)} checksums"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
