"""Streaming validation for the canonical 16-column TBLASTN table."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO
import gzip
import hashlib
import json


BLAST_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qframe",
    "sframe",
    "qlen",
    "slen",
)


class BlastFormatError(ValueError):
    """Raised when a TBLASTN row violates the declared data contract."""


@dataclass(frozen=True, slots=True)
class BlastHit:
    """A validated TBLASTN hit with explicit coordinate semantics."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qframe: int
    sframe: int
    qlen: int
    slen: int

    @property
    def strand(self) -> str:
        """Return the genomic strand encoded by ``sframe``."""

        return "+" if self.sframe > 0 else "-"

    @property
    def genome_start_0(self) -> int:
        """Return the zero-based, inclusive genomic start."""

        return min(self.sstart, self.send) - 1

    @property
    def genome_end_0(self) -> int:
        """Return the zero-based, exclusive genomic end."""

        return max(self.sstart, self.send)


def parse_blast_row(line: str) -> BlastHit:
    """Parse and validate one headerless canonical TBLASTN row."""

    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != len(BLAST_COLUMNS):
        raise BlastFormatError(f"expected 16 tab-separated columns, found {len(fields)}")

    if not fields[0] or not fields[1]:
        raise BlastFormatError("qseqid and sseqid must be non-empty")

    try:
        hit = BlastHit(
            qseqid=fields[0],
            sseqid=fields[1],
            pident=float(fields[2]),
            length=int(fields[3]),
            mismatch=int(fields[4]),
            gapopen=int(fields[5]),
            qstart=int(fields[6]),
            qend=int(fields[7]),
            sstart=int(fields[8]),
            send=int(fields[9]),
            evalue=float(fields[10]),
            bitscore=float(fields[11]),
            qframe=int(fields[12]),
            sframe=int(fields[13]),
            qlen=int(fields[14]),
            slen=int(fields[15]),
        )
    except ValueError as exc:
        raise BlastFormatError(f"invalid numeric value: {exc}") from exc

    if not 0.0 <= hit.pident <= 100.0:
        raise BlastFormatError("pident must be between 0 and 100")
    if hit.length <= 0 or hit.qlen <= 0 or hit.slen <= 0:
        raise BlastFormatError("length, qlen and slen must be positive")
    if hit.mismatch < 0 or hit.gapopen < 0 or hit.evalue < 0 or hit.bitscore < 0:
        raise BlastFormatError("mismatch, gapopen, evalue and bitscore must be non-negative")
    if not (1 <= hit.qstart <= hit.qend <= hit.qlen):
        raise BlastFormatError("query coordinates must satisfy 1 <= qstart <= qend <= qlen")
    if not (1 <= hit.sstart <= hit.slen and 1 <= hit.send <= hit.slen):
        raise BlastFormatError("subject coordinates must fall within 1..slen")
    if hit.qframe != 0:
        raise BlastFormatError("qframe must be 0 for a protein-query TBLASTN result")
    if hit.sframe not in {-3, -2, -1, 1, 2, 3}:
        raise BlastFormatError("sframe must be one of -3,-2,-1,1,2,3")
    if hit.sframe > 0 and hit.sstart > hit.send:
        raise BlastFormatError("positive sframe requires sstart <= send")
    if hit.sframe < 0 and hit.sstart < hit.send:
        raise BlastFormatError("negative sframe requires sstart >= send")

    return hit


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", newline="")
    return path.open(mode="rt", encoding="utf-8", newline="")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open(mode="rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_merged_file(path: Path, max_reported_errors: int = 20) -> dict[str, Any]:
    """Validate a merged TBLASTN file in constant row memory."""

    if max_reported_errors < 1:
        raise ValueError("max_reported_errors must be at least 1")

    total_rows = 0
    valid_rows = 0
    invalid_rows = 0
    plus_rows = 0
    minus_rows = 0
    query_ids: set[str] = set()
    subject_ids: set[str] = set()
    query_lengths: dict[str, int] = {}
    subject_lengths: dict[str, int] = {}
    errors: list[dict[str, Any]] = []

    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            total_rows += 1
            try:
                if not line.strip():
                    raise BlastFormatError("blank lines are not allowed")
                hit = parse_blast_row(line)
                known_qlen = query_lengths.get(hit.qseqid)
                if known_qlen is not None and known_qlen != hit.qlen:
                    raise BlastFormatError(
                        f"inconsistent qlen for {hit.qseqid!r}: {known_qlen} versus {hit.qlen}"
                    )
                known_slen = subject_lengths.get(hit.sseqid)
                if known_slen is not None and known_slen != hit.slen:
                    raise BlastFormatError(
                        f"inconsistent slen for {hit.sseqid!r}: {known_slen} versus {hit.slen}"
                    )
            except BlastFormatError as exc:
                invalid_rows += 1
                if len(errors) < max_reported_errors:
                    errors.append({"line": line_number, "message": str(exc)})
                continue

            valid_rows += 1
            plus_rows += hit.strand == "+"
            minus_rows += hit.strand == "-"
            query_ids.add(hit.qseqid)
            subject_ids.add(hit.sseqid)
            query_lengths[hit.qseqid] = hit.qlen
            subject_lengths[hit.sseqid] = hit.slen

    if total_rows == 0:
        invalid_rows = 1
        errors.append({"line": 0, "message": "input file is empty"})

    return {
        "schema_version": "1.0",
        "status": "pass" if invalid_rows == 0 else "fail",
        "input": str(path.resolve()),
        "input_size_bytes": path.stat().st_size,
        "input_sha256": _sha256(path),
        "format": "headerless canonical tblastn outfmt 6",
        "columns": list(BLAST_COLUMNS),
        "coordinate_contract": {
            "source": "BLAST 1-based inclusive",
            "internal": "0-based half-open",
            "strand_field": "sframe",
        },
        "counts": {
            "total_rows": total_rows,
            "valid_rows": valid_rows,
            "invalid_rows": invalid_rows,
            "plus_strand_rows": plus_rows,
            "minus_strand_rows": minus_rows,
            "unique_query_ids": len(query_ids),
            "unique_subject_ids": len(subject_ids),
        },
        "errors": errors,
        "errors_truncated": invalid_rows > len(errors),
    }


def write_validation_report(report: dict[str, Any], output: Path) -> None:
    """Write a deterministic, human-readable JSON validation report."""

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open(mode="wt", encoding="utf-8", newline="\n") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
