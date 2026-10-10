"""Tests for audited pre-7TM locus consolidation."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from or_mining.fasta import read_fasta
from or_mining.preliminary_catalog import build_preliminary_catalog


def _inputs(tmp_path: Path, overlapping: bool = False) -> dict[str, Path]:
    paths = {name: tmp_path / name for name in (
        "genome.fa.fai", "single.nt.fa", "single.prot.fa", "single.tsv", "multi.tsv",
        "multi.gff3", "multi.functional.nt.fa", "multi.functional.prot.fa",
        "multi.pseudo.nt.fa", "multi.truncated.nt.fa", "multi.ambiguous.nt.fa",
    )}
    paths["genome.fa.fai"].write_text("chr1\t10000\t0\t0\t0\n", encoding="utf-8")
    # Identical sequences at distinct genomic loci represent distinct paralogous loci.
    paths["single.nt.fa"].write_text(
        ">s1\nATGGCTGCTTAA\n>s2\nATGGCTGCTTAA\n", encoding="utf-8"
    )
    paths["single.prot.fa"].write_text(">s1\nMAA\n>s2\nMAA\n", encoding="utf-8")
    second_start = 205 if overlapping else 501
    paths["single.tsv"].write_text(
        "orf_id\tsseqid\tstart_1\tend_1\tstrand\tlength_nt\n"
        "s1\tchr1\t1\t12\t+\t12\n"
        f"s2\tchr1\t{second_start}\t{second_start + 11}\t+\t12\n",
        encoding="utf-8",
    )
    paths["multi.tsv"].write_text(
        "candidate_id\tquery_id\tsseqid\tstrand\tclassification\tevidence\t"
        "genomic_start_1\tgenomic_end_1\n"
        "m1\tq1\tchr1\t+\tfunctional\tintact_start_stop_orf\t101\t212\n"
        "m2\tq2\tchr1\t-\tpseudogene\tpremature_stop\t301\t412\n",
        encoding="utf-8",
    )
    paths["multi.gff3"].write_text(
        "##gff-version 3\n"
        "chr1\tfelidae_or\tgene\t101\t212\t.\t+\t.\tID=m1;classification=functional\n"
        "chr1\tfelidae_or\texon\t101\t106\t.\t+\t.\tID=m1.exon1;Parent=m1\n"
        "chr1\tfelidae_or\tCDS\t101\t106\t.\t+\t0\tParent=m1\n"
        "chr1\tfelidae_or\texon\t207\t212\t.\t+\t.\tID=m1.exon2;Parent=m1\n"
        "chr1\tfelidae_or\tCDS\t207\t212\t.\t+\t0\tParent=m1\n"
        "chr1\tfelidae_or\tgene\t301\t412\t.\t-\t.\tID=m2;classification=pseudogene\n"
        "chr1\tfelidae_or\texon\t407\t412\t.\t-\t.\tID=m2.exon1;Parent=m2\n"
        "chr1\tfelidae_or\tCDS\t407\t412\t.\t-\t.\tParent=m2\n"
        "chr1\tfelidae_or\texon\t301\t306\t.\t-\t.\tID=m2.exon2;Parent=m2\n"
        "chr1\tfelidae_or\tCDS\t301\t306\t.\t-\t.\tParent=m2\n",
        encoding="utf-8",
    )
    paths["multi.functional.nt.fa"].write_text(">m1\nATGGCTGCTTAA\n", encoding="utf-8")
    paths["multi.functional.prot.fa"].write_text(">m1\nMAA\n", encoding="utf-8")
    paths["multi.pseudo.nt.fa"].write_text(">m2\nATGGCTGCTTAG\n", encoding="utf-8")
    paths["multi.truncated.nt.fa"].write_text("", encoding="utf-8")
    paths["multi.ambiguous.nt.fa"].write_text("", encoding="utf-8")
    return paths


def _run(paths: dict[str, Path], prefix: Path) -> dict[str, object]:
    return build_preliminary_catalog(
        "Felid", paths["genome.fa.fai"], paths["single.nt.fa"],
        paths["single.prot.fa"], paths["single.tsv"], paths["multi.tsv"],
        paths["multi.gff3"], paths["multi.functional.nt.fa"],
        paths["multi.functional.prot.fa"], paths["multi.pseudo.nt.fa"],
        paths["multi.truncated.nt.fa"], paths["multi.ambiguous.nt.fa"], prefix,
    )


def test_pre7tm_catalog_retains_distinct_loci_and_preserves_classes(tmp_path: Path) -> None:
    prefix = tmp_path / "Felid"
    summary = _run(_inputs(tmp_path), prefix)
    assert summary["reportable_final_catalog"] is False
    assert summary["counts"] == {
        "input_single_exon_functional_candidates": 2,
        "input_multiexon_models": 2,
        "functional_candidates_requiring_7tm": 3,
        "pseudogene": 1,
        "truncated": 0,
        "ambiguous": 0,
        "catalog_loci": 4,
        "input_overlapping_locus_pairs": 0,
        "single_exon_redundancy_components": 0,
        "redundant_single_exon_models_removed": 0,
        "invalid_cross_stage_or_multiexon_overlaps": 0,
        "remaining_overlapping_locus_pairs": 0,
        "sequence_only_deduplications": 0,
    }
    with Path(f"{prefix}.pre7tm_OR_catalog.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["source_id"] for row in rows] == ["s1", "m1", "m2", "s2"]
    assert [row["classification"] for row in rows].count("functional_candidate") == 3
    proteins = list(read_fasta(Path(f"{prefix}.pre7tm_functional_candidates.prot.fa")))
    assert len(proteins) == 3
    assert len({sequence for _, sequence in proteins}) == 1
    gff = Path(f"{prefix}.pre7tm_OR_catalog.gff3").read_text(encoding="utf-8")
    assert "\tgene\t101\t212\t.\t+\t." in gff
    assert "\tgene\t301\t412\t.\t-\t." in gff


def test_pre7tm_catalog_refuses_silent_overlap_deduplication(tmp_path: Path) -> None:
    prefix = tmp_path / "Felid"
    with pytest.raises(ValueError, match="overlap outside"):
        _run(_inputs(tmp_path, overlapping=True), prefix)
    audit = Path(f"{prefix}.pre7tm_overlap_audit.tsv").read_text(encoding="utf-8")
    assert "m1\ts2\tchr1\t101\t212\t205\t216\t8" in audit
    assert "fatal_cross_stage_or_multiexon_conflict" in audit


def test_pre7tm_catalog_collapses_overlapping_single_exon_component_by_longest(
    tmp_path: Path,
) -> None:
    paths = _inputs(tmp_path)
    paths["single.nt.fa"].write_text(
        ">short\nATGGCTGCTTAA\n>long\nATGGCTGCTGCTTAA\n", encoding="utf-8"
    )
    paths["single.prot.fa"].write_text(">short\nMAA\n>long\nMAAA\n", encoding="utf-8")
    paths["single.tsv"].write_text(
        "orf_id\tsseqid\tstart_1\tend_1\tstrand\tlength_nt\n"
        "short\tchr1\t1\t12\t+\t12\n"
        "long\tchr1\t1\t15\t+\t15\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "Felid"
    summary = _run(paths, prefix)
    assert summary["counts"]["input_overlapping_locus_pairs"] == 1
    assert summary["counts"]["single_exon_redundancy_components"] == 1
    assert summary["counts"]["redundant_single_exon_models_removed"] == 1
    assert summary["counts"]["functional_candidates_requiring_7tm"] == 2
    audit = Path(f"{prefix}.pre7tm_overlap_audit.tsv").read_text(encoding="utf-8")
    assert "collapse_single_exon_component_keep_longest\tlong\tshort" in audit
    proteins = list(read_fasta(Path(f"{prefix}.pre7tm_functional_candidates.prot.fa")))
    assert {header.split()[0] for header, _ in proteins} == {"Felid_ORC000001", "Felid_ORC000002"}
    assert {sequence for _, sequence in proteins} == {"MAAA", "MAA"}
