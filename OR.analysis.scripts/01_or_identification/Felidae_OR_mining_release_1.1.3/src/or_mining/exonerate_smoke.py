"""Generate and validate a tiny positive-control Exonerate protein2genome test."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any

from or_mining.fasta import write_fasta_record


CODONS = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}


def write_smoke_inputs(
    query_path: Path,
    target_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    rng = random.Random(20260810)
    protein_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    protein = "M" + "".join(rng.choice(protein_alphabet) for _ in range(89))
    cds = "".join(CODONS[residue] for residue in protein)
    split_nt = 45 * 3
    intron = "GT" + "".join(rng.choice("ACGT") for _ in range(56)) + "AG"
    left_flank = "".join(rng.choice("ACGT") for _ in range(40))
    right_flank = "".join(rng.choice("ACGT") for _ in range(40))
    target = left_flank + cds[:split_nt] + intron + cds[split_nt:] + right_flank

    with query_path.open("wt", encoding="utf-8", newline="\n") as handle:
        write_fasta_record(handle, "synthetic_spliced_or", protein)
    with target_path.open("wt", encoding="utf-8", newline="\n") as handle:
        write_fasta_record(handle, "synthetic_target", target)

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "test_type": "synthetic_positive_control",
        "counts": {
            "query_records": 1,
            "query_residues": len(protein),
            "target_records": 1,
            "target_bases": len(target),
            "expected_exons": 2,
            "synthetic_intron_bases": len(intron),
        },
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def validate_smoke_result(
    result_path: Path,
    chunk_summary_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    chunk = json.loads(chunk_summary_path.read_text(encoding="utf-8"))
    counts = chunk.get("counts", {})
    text = result_path.read_text(encoding="utf-8", errors="replace")
    checks = {
        "chunk_status_pass": chunk.get("status") == "pass",
        "completed_before_timeout": chunk.get("termination") == "completed",
        "complete_alignment_present": counts.get("retained_complete_alignment_blocks", 0) > 0,
        "vulgar_present": counts.get("retained_vulgar_alignments", 0) > 0,
        "gff_present": counts.get("retained_gff_dumps", 0) > 0,
        "spliced_intron_present": any(
            line.startswith("vulgar:") and " I " in line for line in text.splitlines()
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("Exonerate functional smoke test failed: " + ", ".join(failed))

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "test_type": "exonerate_functional_smoke",
        "checks": checks,
        "source_chunk_summary": chunk_summary_path.name,
        "counts": {
            "retained_complete_alignment_blocks": counts[
                "retained_complete_alignment_blocks"
            ],
            "retained_vulgar_alignments": counts["retained_vulgar_alignments"],
            "retained_gff_dumps": counts["retained_gff_dumps"],
        },
    }
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def check_smoke_gate(summary_path: Path, output_path: Path) -> None:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "pass":
        raise ValueError("functional smoke summary is not pass")
    if summary.get("test_type") != "exonerate_functional_smoke":
        raise ValueError("unexpected functional smoke summary type")
    checks = summary.get("checks")
    if not isinstance(checks, dict) or not checks or not all(checks.values()):
        raise ValueError("functional smoke summary contains failed checks")
    output_path.write_text("status=pass\n", encoding="utf-8")
