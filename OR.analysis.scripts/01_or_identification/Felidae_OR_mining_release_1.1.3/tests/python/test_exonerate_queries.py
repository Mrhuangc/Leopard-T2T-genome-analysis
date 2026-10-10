from pathlib import Path

from or_mining.exonerate_queries import prepare_exonerate_queries
from or_mining.fasta import read_fasta


def test_prepare_exonerate_queries_combines_and_balances(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    species = tmp_path / "species.fa"
    reference.write_text(">r1\nAAAA\n>r2\nAAAAAA\n", encoding="utf-8")
    species.write_text(">s1\nAAA\n>s2\nAAAAA\n", encoding="utf-8")
    output_dir = tmp_path / "chunks"
    summary = prepare_exonerate_queries(
        reference, species, output_dir, tmp_path / "summary.json", requested_chunks=2
    )
    records = [record for path in sorted(output_dir.iterdir()) for record in read_fasta(path)]
    assert {header for header, _ in records} == {"r1", "r2", "s1", "s2"}
    assert summary["counts"]["combined_query_records"] == 4
    assert summary["counts"]["query_chunks"] == 2
