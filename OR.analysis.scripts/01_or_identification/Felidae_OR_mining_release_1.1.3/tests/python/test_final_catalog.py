import csv
import json
from pathlib import Path

import pytest

from or_mining.fasta import read_fasta
from or_mining.final_catalog import build_final_catalog


def _write_inputs(tmp_path: Path) -> dict[str, Path]:
    paths = {name: tmp_path / name for name in (
        "catalog.tsv", "catalog.gff3", "functional.nt.fa", "functional.prot.fa",
        "pseudogene.fa", "truncated.fa", "ambiguous.fa", "interpro.tsv", "version.txt",
    )}
    fields = (
        "catalog_id", "source_id", "source_stage", "classification", "evidence",
        "query_id", "sseqid", "start_1", "end_1", "strand", "genomic_span_bp",
        "spliced_cds_bp", "protein_aa", "exon_count", "requires_7tm_gate",
    )
    rows = [
        ("Felid_ORC000001", "s1", "tree_filtered_single_exon", "functional_candidate", "tree", "", "chr1", "1", "12", "+", "12", "12", "3", "1", "true"),
        ("Felid_ORC000002", "s2", "tree_filtered_single_exon", "functional_candidate", "tree", "", "chr1", "21", "32", "+", "12", "12", "3", "1", "true"),
        ("Felid_ORP000001", "p1", "strict_multiexon", "pseudogene", "stop", "q1", "chr1", "41", "60", "+", "20", "12", "", "2", "false"),
        ("Felid_ORT000001", "t1", "strict_multiexon", "truncated", "edge", "q2", "chr1", "71", "90", "+", "20", "12", "", "2", "false"),
        ("Felid_ORA000001", "a1", "strict_multiexon", "ambiguous", "N", "q3", "chr1", "101", "120", "+", "20", "12", "", "2", "false"),
    ]
    with paths["catalog.tsv"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    paths["catalog.gff3"].write_text(
        "##gff-version 3\n"
        + "".join(
            f"chr1\tfelidae_or\tgene\t{row[7]}\t{row[8]}\t.\t+\t.\t"
            f"ID={row[0]};source_id={row[1]};classification={row[3]};source_stage={row[2]}\n"
            for row in rows
        ),
        encoding="utf-8",
    )
    paths["functional.nt.fa"].write_text(
        ">Felid_ORC000001\nATGGCTGCTTAA\n>Felid_ORC000002\nATGGCCGCCTAA\n",
        encoding="utf-8",
    )
    paths["functional.prot.fa"].write_text(
        ">Felid_ORC000001\nMAA\n>Felid_ORC000002\nMAA\n", encoding="utf-8"
    )
    paths["pseudogene.fa"].write_text(">Felid_ORP000001\nATGTAATGGTAA\n", encoding="utf-8")
    paths["truncated.fa"].write_text(">Felid_ORT000001\nATGGCTGCTTAA\n", encoding="utf-8")
    paths["ambiguous.fa"].write_text(">Felid_ORA000001\nATGNNNGCTTAA\n", encoding="utf-8")
    paths["interpro.tsv"].write_text(
        "Felid_ORC000001\tmd5\t3\tProSiteProfiles\tPS50262\tGPCR profile\t1\t3\t1e-8\tT\t2026-08-12\tIPR017452\tGPCR, rhodopsin-like, 7TM\n"
        "Felid_ORC000002\tmd5\t3\tProSiteProfiles\tPS00001\tother\t1\t3\t1e-2\tT\t2026-08-12\tIPR000001\tother\n",
        encoding="utf-8",
    )
    paths["version.txt"].write_text(
        "InterProScan version 5.78-109.0\n", encoding="utf-8"
    )
    return paths


def _run(paths: dict[str, Path], prefix: Path):
    return build_final_catalog(
        "Felid", paths["catalog.tsv"], paths["catalog.gff3"],
        paths["functional.nt.fa"], paths["functional.prot.fa"],
        paths["pseudogene.fa"], paths["truncated.fa"], paths["ambiguous.fa"],
        paths["interpro.tsv"], paths["version.txt"], "5.78-109.0", prefix,
    )


def test_final_catalog_applies_ipr_gate_without_dropping_loci(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    prefix = tmp_path / "Felid"
    summary = _run(paths, prefix)
    assert summary["status"] == "pass"
    assert summary["reportable_final_catalog"] is True
    assert summary["counts"]["functional"] == 1
    assert summary["counts"]["truncated"] == 2
    assert summary["counts"]["catalog_loci"] == 5

    with Path(f"{prefix}.final_OR_catalog.tsv").open(
        "rt", encoding="utf-8", newline=""
    ) as handle:
        rows = {row["catalog_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["Felid_ORC000001"]["classification"] == "functional"
    assert rows["Felid_ORC000001"]["interpro_7tm_status"] == "pass"
    assert rows["Felid_ORC000002"]["classification"] == "truncated"
    assert rows["Felid_ORC000002"]["pre7tm_classification"] == "functional_candidate"
    assert rows["Felid_ORP000001"]["interpro_7tm_status"] == "not_applicable"

    functional = list(read_fasta(Path(f"{prefix}.final_functional_OR.prot.fa")))
    rejected = list(read_fasta(Path(f"{prefix}.non7tm_rejected.prot.fa")))
    truncated = list(read_fasta(Path(f"{prefix}.final_truncated_OR.nucl.fa")))
    assert [header.split()[0] for header, _ in functional] == ["Felid_ORC000001"]
    assert [header.split()[0] for header, _ in rejected] == ["Felid_ORC000002"]
    assert {header.split()[0] for header, _ in truncated} == {
        "Felid_ORC000002", "Felid_ORT000001"
    }
    gff = Path(f"{prefix}.final_OR_catalog.gff3").read_text(encoding="utf-8")
    assert "ID=Felid_ORC000001;source_id=s1;classification=functional" in gff
    assert "ID=Felid_ORC000002;source_id=s2;classification=truncated" in gff
    assert "interpro_7tm=fail" in gff


def test_final_catalog_refuses_scan_without_expected_ipr(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    paths["interpro.tsv"].write_text(
        "Felid_ORC000002\tmd5\t3\tProSiteProfiles\tPS00001\tother\t1\t3\t1e-2\tT\t2026-08-12\tIPR000001\tother\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="contains no IPR017452"):
        _run(paths, tmp_path / "failed")


def test_final_catalog_summary_is_serializable(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    prefix = tmp_path / "Felid"
    _run(paths, prefix)
    summary = json.loads(
        Path(f"{prefix}.final_OR_catalog.summary.json").read_text(encoding="utf-8")
    )
    assert summary["checks"]["catalog_partition_is_complete"] is True
