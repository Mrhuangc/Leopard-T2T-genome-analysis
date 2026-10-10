import csv
import json
from pathlib import Path

import pytest

from or_mining.fasta import read_fasta
from or_mining.final_homology import build_nc_final_catalog


def _blast_row(query: str, subject: str, title: str) -> str:
    return (
        f"{query}\t{subject}\t80.0\t3\t0\t0\t1\t3\t1\t3\t1e-20\t100\t3\t3\t{title}\n"
    )


def _write_inputs(tmp_path: Path) -> dict[str, Path]:
    names = (
        "catalog.tsv", "catalog.gff3", "functional.nt.fa", "functional.prot.fa",
        "pseudogene.fa", "truncated.fa", "ambiguous.fa", "functional.blastp.tsv",
        "pseudogene.blastx.tsv", "truncated.blastx.tsv", "ambiguous.blastx.tsv",
    )
    paths = {name: tmp_path / name for name in names}
    fields = (
        "catalog_id", "source_id", "classification", "evidence", "sseqid",
        "start_1", "end_1", "strand", "exon_count",
    )
    rows = [
        ("Felid_ORC1", "s1", "functional", "7tm", "chr1", "1", "12", "+", "1"),
        ("Felid_ORC2", "s2", "functional", "7tm", "chr1", "21", "32", "+", "1"),
        ("Felid_ORP1", "p1", "pseudogene", "stop", "chr1", "41", "52", "+", "2"),
        ("Felid_ORT1", "t1", "truncated", "edge", "chr1", "61", "72", "+", "2"),
        ("Felid_ORA1", "a1", "ambiguous", "N", "chr1", "81", "92", "+", "2"),
    ]
    with paths["catalog.tsv"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    paths["catalog.gff3"].write_text(
        "##gff-version 3\n" + "".join(
            f"chr1\tfelidae_or\tgene\t{row[5]}\t{row[6]}\t.\t+\t.\t"
            f"ID={row[0]};classification={row[2]}\n" for row in rows
        ),
        encoding="utf-8",
    )
    paths["functional.nt.fa"].write_text(
        ">Felid_ORC1\nATGGCTGCTTAA\n>Felid_ORC2\nATGGCTGCTTAA\n", encoding="utf-8"
    )
    paths["functional.prot.fa"].write_text(
        ">Felid_ORC1\nMAA\n>Felid_ORC2\nMAA\n", encoding="utf-8"
    )
    paths["pseudogene.fa"].write_text(">Felid_ORP1\nATGTAATGGTAA\n", encoding="utf-8")
    paths["truncated.fa"].write_text(">Felid_ORT1\nATGGCTGCTTAA\n", encoding="utf-8")
    paths["ambiguous.fa"].write_text(">Felid_ORA1\nATGNNNGCTTAA\n", encoding="utf-8")
    paths["functional.blastp.tsv"].write_text(
        _blast_row("Felid_ORC1", "OR1", "OR-Receptor family protein")
        + _blast_row("Felid_ORC2", "OPS1", "opsin receptor"),
        encoding="utf-8",
    )
    paths["pseudogene.blastx.tsv"].write_text(
        _blast_row("Felid_ORP1", "OR2", "OR-Receptor pseudogene homolog"), encoding="utf-8"
    )
    paths["truncated.blastx.tsv"].write_text("", encoding="utf-8")
    paths["ambiguous.blastx.tsv"].write_text(
        _blast_row("Felid_ORA1", "GF1", "growth-factor-receptor-like protein"),
        encoding="utf-8",
    )
    return paths


def _run(paths: dict[str, Path], prefix: Path):
    return build_nc_final_catalog(
        "Felid", paths["catalog.tsv"], paths["catalog.gff3"],
        paths["functional.nt.fa"], paths["functional.prot.fa"],
        paths["pseudogene.fa"], paths["truncated.fa"], paths["ambiguous.fa"],
        paths["functional.blastp.tsv"], paths["pseudogene.blastx.tsv"],
        paths["truncated.blastx.tsv"], paths["ambiguous.blastx.tsv"],
        prefix, "1e-5",
    )


def test_final_homology_filter_conserves_every_locus(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    prefix = tmp_path / "Felid"
    summary = _run(paths, prefix)
    assert summary["status"] == "pass"
    assert summary["reportable_final_catalog"] is True
    assert summary["counts"]["input_catalog_loci"] == 5
    assert summary["counts"]["retained_catalog_loci"] == 2
    assert summary["counts"]["excluded_non_or_loci"] == 3
    assert summary["counts"]["retained_functional"] == 1
    assert summary["counts"]["retained_pseudogene"] == 1

    with Path(f"{prefix}.NC_final_OR_catalog.tsv").open(
        "rt", encoding="utf-8", newline=""
    ) as handle:
        retained = {row["catalog_id"] for row in csv.DictReader(handle, delimiter="\t")}
    assert retained == {"Felid_ORC1", "Felid_ORP1"}
    with Path(f"{prefix}.NC_excluded_nonOR_catalog.tsv").open(
        "rt", encoding="utf-8", newline=""
    ) as handle:
        rejected = {row["catalog_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rejected["Felid_ORC2"]["final_homology_reason"] == "non_or_best_hit"
    assert rejected["Felid_ORT1"]["final_homology_reason"] == "no_hit"
    assert rejected["Felid_ORA1"]["final_homology_reason"] == "growth_factor_receptor"

    proteins = list(read_fasta(Path(f"{prefix}.NC_final_functional_OR.prot.fa")))
    excluded_nt = list(read_fasta(Path(f"{prefix}.NC_excluded_nonOR.nucl.fa")))
    assert [header.split()[0] for header, _ in proteins] == ["Felid_ORC1"]
    assert {header.split()[0] for header, _ in excluded_nt} == {
        "Felid_ORC2", "Felid_ORT1", "Felid_ORA1"
    }
    gff = Path(f"{prefix}.NC_final_OR_catalog.gff3").read_text(encoding="utf-8")
    assert "ID=Felid_ORC1" in gff
    assert "ID=Felid_ORP1" in gff
    assert "ID=Felid_ORC2" not in gff
    assert "final_homology=pass" in gff


def test_final_homology_summary_is_serializable(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    prefix = tmp_path / "Felid"
    _run(paths, prefix)
    summary = json.loads(
        Path(f"{prefix}.NC_final_OR_catalog.summary.json").read_text(encoding="utf-8")
    )
    assert summary["checks"]["catalog_partition_is_complete"] is True


def test_final_homology_refuses_zero_functional_retention(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    paths["functional.blastp.tsv"].write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="retained no functional ORs"):
        _run(paths, tmp_path / "failed")
