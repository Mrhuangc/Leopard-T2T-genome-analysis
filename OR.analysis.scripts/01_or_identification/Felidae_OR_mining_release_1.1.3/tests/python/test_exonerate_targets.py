from pathlib import Path

from or_mining.exonerate_targets import prepare_exonerate_targets
from or_mining.fasta import read_fasta


def test_prepare_exonerate_targets_balances_without_losing_records(tmp_path: Path) -> None:
    target = tmp_path / "targets.fa"
    target.write_text(
        ">r1\n" + "A" * 100 + "\n"
        ">r2\n" + "C" * 80 + "\n"
        ">r3\n" + "G" * 60 + "\n"
        ">r4\n" + "T" * 40 + "\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "chunks"
    summary = prepare_exonerate_targets(
        target,
        output_dir,
        tmp_path / "summary.json",
        requested_chunks=2,
    )
    records = [record for path in sorted(output_dir.iterdir()) for record in read_fasta(path)]
    assert {header for header, _ in records} == {"r1", "r2", "r3", "r4"}
    assert summary["counts"]["target_records"] == 4
    assert summary["counts"]["target_bases"] == 280
    assert summary["counts"]["target_chunks"] == 2
    assert summary["counts"]["maximum_chunk_bases"] == 140
    assert summary["counts"]["minimum_chunk_bases"] == 140
