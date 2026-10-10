"""Deterministic preparation and validation for the initial genome-wide TBLASTN."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from or_mining.blast import validate_merged_file, write_validation_report

INITIAL_OUTDIR = "00_initial_search/merged"
INITIAL_SUFFIX = "initial_merged.tblastn.tsv"


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """A normalized FASTA record with stable source ordering."""

    index: int
    identifier: str
    description: str
    sequence: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    seen: set[str] = set()
    header: str | None = None
    sequence: list[str] = []

    def append_record() -> None:
        if header is None:
            return
        identifier = header.split(maxsplit=1)[0]
        if not identifier:
            raise ValueError(f"empty FASTA identifier in {path}")
        if identifier in seen:
            raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
        normalized = "".join(sequence).replace(" ", "").upper()
        if not normalized:
            raise ValueError(f"empty FASTA sequence for {identifier!r} in {path}")
        if any(character.isspace() for character in normalized):
            raise ValueError(f"embedded whitespace in FASTA sequence {identifier!r}")
        seen.add(identifier)
        records.append(
            FastaRecord(
                index=len(records),
                identifier=identifier,
                description=header,
                sequence=normalized,
            )
        )

    with path.open("rt", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                append_record()
                header = line[1:].strip()
                sequence = []
            elif header is None:
                raise ValueError(f"{path}:{line_number}: sequence occurs before the first header")
            else:
                sequence.append(line)
    append_record()
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def split_initial_queries(
    query_fasta: Path,
    output_dir: Path,
    output_prefix: str,
    requested_shards: int,
) -> dict[str, Any]:
    """Length-balance protein records into deterministic query shards."""

    if requested_shards < 1:
        raise ValueError("requested_shards must be positive")
    records = _read_fasta(query_fasta)
    shard_count = min(requested_shards, len(records))
    shards: list[list[FastaRecord]] = [[] for _ in range(shard_count)]
    residues = [0] * shard_count

    # Longest-processing-time allocation provides stable, well-balanced TBLASTN tasks.
    for record in sorted(records, key=lambda value: (-len(value.sequence), value.index)):
        target = min(
            range(shard_count),
            key=lambda index: (residues[index], len(shards[index]), index),
        )
        shards[target].append(record)
        residues[target] += len(record.sequence)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(f"{output_prefix}.initial_query_manifest.tsv")
    summary_path = Path(f"{output_prefix}.initial_query_shards.summary.json")
    fields = ("shard_id", "fasta", "records", "residues", "first_query", "last_query")
    manifest_rows: list[dict[str, Any]] = []

    for index, values in enumerate(shards, start=1):
        ordered = sorted(values, key=lambda value: value.index)
        shard_id = f"shard_{index:03d}"
        shard_path = output_dir / f"{shard_id}.prot.fa"
        with shard_path.open("wt", encoding="utf-8", newline="\n") as handle:
            for record in ordered:
                handle.write(f">{record.description}\n{record.sequence}\n")
        manifest_rows.append(
            {
                "shard_id": shard_id,
                "fasta": shard_path.name,
                "records": len(ordered),
                "residues": sum(len(record.sequence) for record in ordered),
                "first_query": ordered[0].identifier,
                "last_query": ordered[-1].identifier,
            }
        )

    with manifest_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "method": "deterministic longest-processing-time protein query sharding",
        "input": {
            "query_fasta": str(query_fasta.resolve()),
            "query_sha256": _sha256(query_fasta),
        },
        "parameters": {"requested_shards": requested_shards},
        "counts": {
            "query_records": len(records),
            "query_residues": sum(len(record.sequence) for record in records),
            "query_shards": shard_count,
            "minimum_shard_records": min(len(values) for values in shards),
            "maximum_shard_records": max(len(values) for values in shards),
            "minimum_shard_residues": min(residues),
            "maximum_shard_residues": max(residues),
        },
        "checks": {
            "all_query_ids_unique": True,
            "all_queries_assigned_once": sum(len(values) for values in shards) == len(records),
            "no_empty_shards": all(shards),
        },
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def merge_initial_tblastn(
    part_files: list[Path],
    query_manifest: Path,
    output_merged: Path,
    validation_path: Path,
    summary_path: Path,
    *,
    species: str,
    evalue: str,
) -> dict[str, Any]:
    """Merge completed query-shard outputs in stable order and validate all rows."""

    ordered = sorted(part_files, key=lambda path: path.name)
    if not ordered:
        raise ValueError("no initial TBLASTN part files were supplied")
    if len({path.name for path in ordered}) != len(ordered):
        raise ValueError("duplicate initial TBLASTN part filenames")
    missing = [str(path) for path in ordered if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"initial TBLASTN parts are missing: {missing[:5]}")

    with query_manifest.open("rt", encoding="utf-8", newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_shards = [row.get("shard_id", "") for row in manifest_rows]
    if not expected_shards or any(not value for value in expected_shards):
        raise ValueError("initial query manifest is empty or malformed")
    if len(set(expected_shards)) != len(expected_shards):
        raise ValueError("initial query manifest contains duplicate shard IDs")
    actual_shards: list[str] = []
    for path in ordered:
        match = re.search(r"\.(shard_\d+)\.initial\.tblastn\.tsv$", path.name)
        if match is None:
            raise ValueError(f"cannot infer shard ID from TBLASTN part filename: {path.name}")
        actual_shards.append(match.group(1))
    missing_shards = sorted(set(expected_shards) - set(actual_shards))
    unexpected_shards = sorted(set(actual_shards) - set(expected_shards))
    if missing_shards or unexpected_shards or len(actual_shards) != len(expected_shards):
        raise ValueError(
            "initial TBLASTN shard set differs from the query manifest: "
            f"missing={missing_shards}; unexpected={unexpected_shards}"
        )

    output_merged.parent.mkdir(parents=True, exist_ok=True)
    total_bytes = 0
    with output_merged.open("wb") as output:
        for path in ordered:
            payload = path.read_bytes()
            total_bytes += len(payload)
            output.write(payload)
            if payload and not payload.endswith(b"\n"):
                output.write(b"\n")
    if output_merged.stat().st_size == 0:
        raise ValueError("merged initial TBLASTN output is empty")

    validation = validate_merged_file(output_merged)
    write_validation_report(validation, validation_path)
    if validation["status"] != "pass":
        raise ValueError("merged initial TBLASTN output failed the canonical 16-column contract")

    counts = validation["counts"]
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "species": species,
        "stage": "initial reference-OR protein versus genome TBLASTN",
        "search_completeness": "complete",
        "method": "query-sharded whole-genome TBLASTN followed by ordered byte merge",
        "parameters": {
            "evalue": str(evalue),
            "outfmt": (
                "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
                "evalue bitscore qframe sframe qlen slen"
            ),
        },
        "counts": {
            "input_query_shards": len(expected_shards),
            "completed_query_shards": len(ordered),
            "merged_source_bytes": total_bytes,
            "merged_rows": counts["total_rows"],
            "unique_query_ids_with_hits": counts["unique_query_ids"],
            "unique_subject_ids_with_hits": counts["unique_subject_ids"],
            "plus_strand_rows": counts["plus_strand_rows"],
            "minus_strand_rows": counts["minus_strand_rows"],
        },
        "checks": {
            "all_query_shards_completed": True,
            "query_manifest_matches_parts": set(expected_shards) == set(actual_shards),
            "merged_output_nonempty": counts["total_rows"] > 0,
            "canonical_16_column_contract": validation["status"] == "pass",
            "invalid_rows": counts["invalid_rows"],
        },
        "outputs": {
            "merged": output_merged.name,
            "validation": validation_path.name,
            "summary": summary_path.name,
        },
        "query_manifest_sha256": _sha256(query_manifest),
        "merged_sha256": _sha256(output_merged),
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def generated_merged_path(params: dict[str, Any]) -> Path:
    """Return the deterministic published merged path for a generated search."""

    outdir = Path(str(params.get("outdir") or "results"))
    species = str(params.get("species") or "")
    return outdir / INITIAL_OUTDIR / f"{species}.{INITIAL_SUFFIX}"


def resolve_merged_input(params: dict[str, Any]) -> tuple[Path, str]:
    """Resolve reuse/generate/auto mode without silently changing scientific inputs."""

    mode = str(params.get("initial_search_mode") or "reuse").strip().lower()
    if mode not in {"reuse", "generate", "auto"}:
        raise ValueError("initial_search_mode must be reuse, generate, or auto")
    supplied = str(params.get("input_merged") or "").strip()
    supplied_path = Path(supplied).expanduser() if supplied else None
    generated = generated_merged_path(params)
    if mode == "reuse":
        if supplied_path is None:
            raise ValueError("input_merged is required when initial_search_mode=reuse")
        return supplied_path, "supplied"
    if mode == "generate":
        return generated, "generated"
    if supplied_path is not None and supplied_path.is_file():
        return supplied_path, "supplied"
    return generated, "generated"


def write_resolved_params(params_file: Path, output: Path) -> dict[str, Any]:
    """Write an immutable run snapshot with the effective merged input made explicit."""

    import yaml

    from or_mining.preflight import load_params

    params = load_params(params_file)
    merged, source = resolve_merged_input(params)
    if not merged.is_file():
        raise FileNotFoundError(f"resolved merged TBLASTN file does not exist: {merged}")
    params["input_merged"] = str(merged.resolve())
    params["resolved_initial_input_source"] = source
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        yaml.safe_dump(params, sort_keys=False, allow_unicode=True), encoding="utf-8", newline="\n"
    )
    return params
