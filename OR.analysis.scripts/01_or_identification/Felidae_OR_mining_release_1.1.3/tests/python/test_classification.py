"""Tests for translation normalization and traceable GPCR classification."""

from pathlib import Path

from or_mining.classification import classify_best_gpcr_hits
from or_mining.translations import normalize_transeq


def test_normalize_transeq_removes_terminal_stop(tmp_path: Path) -> None:
    nucleotide = tmp_path / "orfs.nt.fa"
    raw_protein = tmp_path / "raw.prot.fa"
    protein = tmp_path / "orfs.prot.fa"
    summary = tmp_path / "translation.json"
    nucleotide.write_text(">sp_ORF000001\nATGAAATAA\n", encoding="utf-8")
    raw_protein.write_text(">sp_ORF000001_1\nMK*\n", encoding="utf-8")

    result = normalize_transeq(nucleotide, raw_protein, protein, summary)

    assert protein.read_text(encoding="utf-8") == ">sp_ORF000001\nMK\n"
    assert result["counts"]["terminal_stops_removed"] == 1


def test_classifies_every_orf_and_retains_only_or_best_hits(tmp_path: Path) -> None:
    proteins = tmp_path / "orfs.prot.fa"
    nucleotides = tmp_path / "orfs.nt.fa"
    coordinates = tmp_path / "orfs.tsv"
    blastp = tmp_path / "blastp.tsv"
    prefix = tmp_path / "Ppa"
    ids = ["Ppa_ORF000001", "Ppa_ORF000002", "Ppa_ORF000003", "Ppa_ORF000004"]
    proteins.write_text("".join(f">{identifier}\nMK\n" for identifier in ids), encoding="utf-8")
    nucleotides.write_text(
        "".join(f">{identifier}\nATGAAATAA\n" for identifier in ids), encoding="utf-8"
    )
    coordinates.write_text(
        "orf_id\tsseqid\tstart_1\tend_1\tstrand\n"
        + "".join(f"{identifier}\tChr1\t1\t9\t+\n" for identifier in ids),
        encoding="utf-8",
    )

    def row(query: str, subject: str, title: str) -> str:
        return (
            f"{query}\t{subject}\t90\t2\t0\t0\t1\t2\t1\t2\t1e-20\t100\t2\t300\t{title}\n"
        )

    blastp.write_text(
        row(ids[0], "OR-Receptor-cat-1", "OR-Receptor cat")
        + row(ids[0], "OR-Receptor-cat-1", "OR-Receptor cat secondary HSP")
        + row(ids[1], "rhodopsin", "rhodopsin receptor")
        + row(
            ids[2],
            "OR-Receptor-growth-factor-receptor",
            "OR-Receptor growth-factor-receptor",
        ),
        encoding="utf-8",
    )

    result = classify_best_gpcr_hits(
        blastp, proteins, nucleotides, coordinates, output_prefix=prefix
    )

    assert result["counts"] == {
        "input_orfs": 4,
        "blastp_rows": 4,
        "blastp_best_hits": 3,
        "queries_with_multiple_hsps": 1,
        "retained_or": 1,
        "no_hit": 1,
        "non_or_best_hit": 1,
        "growth_factor": 1,
    }
    assert Path(f"{prefix}.OR_candidates.ids.txt").read_text(encoding="utf-8") == (
        "Ppa_ORF000001\n"
    )
