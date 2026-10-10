"""Minimal strict FASTA I/O used by workflow boundary adapters."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header_without_>, sequence)`` records from a FASTA file."""

    header: str | None = None
    parts: list[str] = []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                if not header:
                    raise ValueError(f"{path}: empty FASTA header at line {line_number}")
                parts = []
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence before first header at line {line_number}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def write_fasta_record(handle: TextIO, header: str, sequence: str, width: int = 60) -> None:
    """Write one wrapped FASTA record to an open text handle."""

    if width < 1:
        raise ValueError("FASTA line width must be positive")
    handle.write(f">{header}\n")
    for offset in range(0, len(sequence), width):
        handle.write(sequence[offset : offset + width] + "\n")
