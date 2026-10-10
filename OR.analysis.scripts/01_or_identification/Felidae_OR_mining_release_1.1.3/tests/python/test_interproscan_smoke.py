import json
from pathlib import Path

import pytest

from or_mining.fasta import read_fasta
from or_mining.interproscan_smoke import (
    prepare_interproscan_smoke,
    validate_interproscan_smoke,
)


def _write_reference(path: Path, count: int = 9) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for index in range(count):
            handle.write(f">OR{index:03d} description\n")
            handle.write("M" + "ACDEFGHIKLMNPQRSTVWY" * 13 + "*\n")


def test_prepare_interproscan_smoke_is_small_and_deterministic(tmp_path: Path) -> None:
    reference = tmp_path / "known_or.fa"
    _write_reference(reference)
    output = tmp_path / "smoke.fa"
    summary = prepare_interproscan_smoke(reference, output, tmp_path / "summary.json", 5)
    records = list(read_fasta(output))
    assert [header for header, _ in records] == ["OR000", "OR002", "OR004", "OR006", "OR008"]
    assert all(not sequence.endswith("*") for _, sequence in records)
    assert summary["counts"]["selected_records"] == 5


def test_validate_interproscan_smoke_requires_ipr017452(tmp_path: Path) -> None:
    input_fasta = tmp_path / "smoke.fa"
    input_fasta.write_text(">OR001\nMACDEFG\n>OR002\nMACDEFG\n", encoding="utf-8")
    result = tmp_path / "result.tsv"
    result.write_text(
        "OR001\tmd5\t7\tPfam\tPF00001\t7tm_1\t1\t7\t1e-10\tT\t2026-08-12\t"
        "IPR017452\tGPCR, rhodopsin-like, 7TM\n",
        encoding="utf-8",
    )
    version = tmp_path / "version.txt"
    version.write_text("InterProScan version 5.78-109.0\n", encoding="utf-8")
    gate = tmp_path / "gate.ok"
    summary = validate_interproscan_smoke(
        input_fasta,
        result,
        version,
        tmp_path / "summary.json",
        gate,
        "5.78-109.0",
    )
    assert summary["counts"]["sequences_with_expected_interpro_accession"] == 1
    assert gate.read_text(encoding="utf-8").startswith("status=pass\n")

    result.write_text(
        "OR001\tmd5\t7\tPfam\tPF99999\tother\t1\t7\t1e-10\tT\t2026-08-12\t-\t-\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="expected_interpro_accession_present"):
        validate_interproscan_smoke(
            input_fasta,
            result,
            version,
            tmp_path / "failed.json",
            tmp_path / "failed.ok",
            "5.78-109.0",
        )
