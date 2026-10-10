import csv
from pathlib import Path

from or_mining.fasta import read_fasta
from or_mining.targeted_exonerate import prepare_targeted_pairs


def _blast_row(
    query: str,
    subject: str,
    qstart: int,
    qend: int,
    sstart: int,
    send: int,
    bitscore: float,
) -> str:
    length = qend - qstart + 1
    return "\t".join(
        map(
            str,
            (
                query,
                subject,
                80.0,
                length,
                0,
                0,
                qstart,
                qend,
                sstart,
                send,
                1e-20,
                bitscore,
                0,
                1,
                100,
                2000,
            ),
        )
    )


def test_prepare_targeted_pairs_replaces_all_by_all_search(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    species = tmp_path / "species.fa"
    targets = tmp_path / "targets.fa"
    regions = tmp_path / "regions.tsv"
    hits = tmp_path / "hits.tsv"
    reference.write_text(
        ">q1\n" + "A" * 100 + "\n>q2\n" + "C" * 100 + "\n>q3\n" + "G" * 100 + "\n",
        encoding="utf-8",
    )
    species.write_text(">s1\n" + "T" * 100 + "\n", encoding="utf-8")
    targets.write_text(">r1\n" + "A" * 1000 + "\n>r2\n" + "C" * 1000 + "\n", encoding="utf-8")
    regions.write_text(
        "region_id\tsseqid\texpanded_start_1\texpanded_end_1\n"
        "r1\tscaffold1\t1\t1000\n"
        "r2\tscaffold1\t1001\t2000\n",
        encoding="utf-8",
    )
    hits.write_text(
        _blast_row("q1", "scaffold1", 1, 50, 100, 249, 100.0)
        + "\n"
        + _blast_row("q2", "scaffold1", 1, 60, 150, 329, 120.0)
        + "\n"
        + _blast_row("q3", "scaffold1", 1, 80, 1200, 1439, 150.0)
        + "\n",
        encoding="utf-8",
    )

    pair_dir = tmp_path / "pairs"
    manifest = tmp_path / "manifest.tsv"
    summary = prepare_targeted_pairs(
        reference,
        species,
        targets,
        regions,
        hits,
        pair_dir,
        manifest,
        tmp_path / "summary.json",
        top_queries_per_region=2,
    )
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["region_id"] == "r2"
    assert rows[0]["best_query_id"] == "q3"
    assert rows[1]["region_id"] == "r1"
    assert rows[1]["selected_query_ids"] == "q2,q1"
    assert len(list(pair_dir.glob("*.query.prot.fa"))) == 2
    assert len(list(pair_dir.glob("*.target.nucl.fa"))) == 2
    assert len(list(read_fasta(pair_dir / "pair_000002.query.prot.fa"))) == 2
    assert summary["counts"]["pair_tasks"] == 2
    assert summary["counts"]["selected_query_records"] == 3
    assert summary["counts"]["all_by_all_comparisons_avoided"] == 5
