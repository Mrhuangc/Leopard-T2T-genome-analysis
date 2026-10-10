import json
from pathlib import Path

from or_mining.exonerate_smoke import check_smoke_gate, validate_smoke_result, write_smoke_inputs
from or_mining.fasta import read_fasta


def test_smoke_inputs_are_small_spliced_positive_control(tmp_path: Path) -> None:
    query = tmp_path / "query.fa"
    target = tmp_path / "target.fa"
    summary = write_smoke_inputs(query, target, tmp_path / "inputs.json")
    assert len(list(read_fasta(query))) == 1
    assert len(list(read_fasta(target))) == 1
    assert summary["counts"]["expected_exons"] == 2
    assert summary["counts"]["target_bases"] < 1000


def test_smoke_validation_and_gate_require_all_positive_signals(tmp_path: Path) -> None:
    result = tmp_path / "result.txt"
    result.write_text(
        "C4 Alignment:\n"
        "vulgar: q 0 121 + t 100 563 + 1000 M 61 183 I 0 100 M 60 180\n"
        "# --- START OF GFF DUMP ---\n"
        "t\texonerate\tgene\t101\t563\t.\t+\t.\tID=g1\n"
        "# --- END OF GFF DUMP ---\n",
        encoding="utf-8",
    )
    chunk = tmp_path / "chunk.json"
    chunk.write_text(
        json.dumps(
            {
                "status": "pass",
                "termination": "completed",
                "counts": {
                    "retained_complete_alignment_blocks": 1,
                    "retained_vulgar_alignments": 1,
                    "retained_gff_dumps": 1,
                },
            }
        ),
        encoding="utf-8",
    )
    audit = tmp_path / "audit.json"
    validate_smoke_result(result, chunk, audit)
    gate = tmp_path / "gate.ok"
    check_smoke_gate(audit, gate)
    assert gate.read_text(encoding="utf-8") == "status=pass\n"
