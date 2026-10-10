import csv
from pathlib import Path

from or_mining.fasta import read_fasta
from or_mining.single_exon_recovery import (
    CANDIDATE_FIELDS,
    build_complete_catalog,
    build_single_exon_candidates,
    classify_recovered_candidates,
)


def _merged_row(query: str, start: int, end: int, length: int, frame: int = 1) -> str:
    sstart, send = (start, end) if frame > 0 else (end, start)
    return (
        f"{query}\tchr1\t80\t{length}\t0\t0\t1\t{length}\t{sstart}\t{send}\t"
        f"1e-30\t100\t0\t{frame}\t100\t2000\n"
    )


def test_prepare_residual_candidates_excludes_catalog_and_reduces_hits(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "catalog_id\tsseqid\tstart_1\tend_1\nOR1\tchr1\t100\t300\n", encoding="utf-8"
    )
    merged = tmp_path / "merged.tsv"
    merged.write_text(
        _merged_row("q1", 150, 300, 50)
        + _merged_row("q1", 500, 700, 60)
        + _merged_row("q2", 650, 800, 80)
        + _merged_row("q1", 1200, 1400, 65, -1),
        encoding="utf-8",
    )
    reference = tmp_path / "reference.fa"
    species = tmp_path / "species.fa"
    reference.write_text(">q1\n" + "M" * 100 + "\n", encoding="utf-8")
    species.write_text(">q2\n" + "A" * 100 + "\n", encoding="utf-8")
    prefix = tmp_path / "Felid"
    summary = build_single_exon_candidates(
        merged, catalog, reference, species, "Felid", prefix, clearance=100, flank=100
    )
    assert summary["counts"]["input_hits"] == 4
    assert summary["counts"]["hits_excluded_by_confirmed_catalog"] == 1
    assert summary["counts"]["single_exon_residual_candidates"] == 2
    with Path(f"{prefix}.single_exon_residual_candidates.tsv").open(
        "rt", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["query_id"] == "q2"
    assert rows[0]["true_start_1"] == "650"
    assert rows[1]["strand"] == "-"


def _write_classification_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    candidates = tmp_path / "candidates.tsv"
    rows = []
    for index, start in enumerate((500, 900, 1300), start=1):
        identifier = f"Felid_ORSE{index:06d}"
        rows.append({
            "candidate_id": identifier, "sseqid": "chr1", "strand": "+",
            "query_id": "q1", "true_start_1": str(start), "true_end_1": str(start + 8),
            "expanded_start_1": str(start - 3), "expanded_end_1": str(start + 11),
            "slen": "2000", "alignment_aa": "3", "pident": "90", "evalue": "1e-20",
            "bitscore": "50", "faidx_region": f"chr1:{start - 3}-{start + 11}",
            "source_hit_count": "1",
        })
    with candidates.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    targets = tmp_path / "targets.fa"
    targets.write_text(
        ">Felid_ORSE000001\nAAAATGAAATAAAAA\n"
        ">Felid_ORSE000002\nAAAATGAAATAAAAA\n"
        ">Felid_ORSE000003\nAAAATGNNNTAAAAA\n",
        encoding="utf-8",
    )
    standard = tmp_path / "standard.tsv"
    standard.write_text(
        "Felid_ORSE000001\tFelid_ORSE000001\t90\t3\t0\t0\t1\t3\t4\t12\t1e-20\t1\tM*K\n"
        "Felid_ORSE000002\tFelid_ORSE000002\t90\t3\t0\t0\t1\t3\t4\t12\t1e-20\t1\tMKK\n"
        "Felid_ORSE000003\tFelid_ORSE000003\t90\t3\t0\t0\t1\t3\t4\t12\t1e-20\t1\tMXX\n",
        encoding="utf-8",
    )
    gapless = tmp_path / "gapless.tsv"
    gapless.write_text(standard.read_text(encoding="utf-8"), encoding="utf-8")
    return candidates, targets, standard, gapless


def test_classification_calls_pseudogene_truncated_and_ambiguous(tmp_path: Path) -> None:
    candidates, targets, standard, gapless = _write_classification_inputs(tmp_path)
    prefix = tmp_path / "Felid"
    summary = classify_recovered_candidates(candidates, targets, standard, gapless, prefix)
    assert summary["counts"]["pseudogene"] == 1
    assert summary["counts"]["truncated"] == 1
    assert summary["counts"]["ambiguous"] == 1
    assert len(list(read_fasta(Path(f"{prefix}.single_exon_recovery_candidates.nucl.fa")))) == 3


def test_complete_catalog_preserves_existing_ids_and_adds_recovered(tmp_path: Path) -> None:
    candidates, targets, standard, gapless = _write_classification_inputs(tmp_path)
    prefix = tmp_path / "classified"
    classify_recovered_candidates(candidates, targets, standard, gapless, prefix)

    fields = (
        "catalog_id", "source_id", "source_stage", "classification", "evidence",
        "query_id", "sseqid", "start_1", "end_1", "strand", "genomic_span_bp",
        "spliced_cds_bp", "protein_aa", "exon_count", "requires_7tm_gate",
        "pre7tm_classification", "interpro_7tm_status", "final_homology_decision",
        "final_homology_reason", "final_homology_program", "final_homology_sseqid",
        "final_homology_stitle", "final_homology_pident",
        "final_homology_alignment_length", "final_homology_evalue",
        "final_homology_bitscore", "final_homology_query_coverage",
    )
    catalog = tmp_path / "existing.tsv"
    with catalog.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            **{field: "" for field in fields}, "catalog_id": "Felid_ORC000001",
            "source_id": "old", "source_stage": "old", "classification": "functional",
            "evidence": "old", "sseqid": "chr1", "start_1": "100", "end_1": "300",
            "strand": "+", "genomic_span_bp": "201", "spliced_cds_bp": "9",
            "protein_aa": "3", "exon_count": "1", "requires_7tm_gate": "false",
        })
    gff = tmp_path / "existing.gff3"
    gff.write_text(
        "##gff-version 3\nchr1\tx\tgene\t100\t300\t.\t+\t.\tID=Felid_ORC000001\n",
        encoding="utf-8",
    )
    existing = {
        "fnt": tmp_path / "fnt.fa", "fprot": tmp_path / "fprot.fa",
        "p": tmp_path / "p.fa", "t": tmp_path / "t.fa", "a": tmp_path / "a.fa",
    }
    existing["fnt"].write_text(">Felid_ORC000001\nATGAAAAAA\n", encoding="utf-8")
    existing["fprot"].write_text(">Felid_ORC000001\nMKK\n", encoding="utf-8")
    for key in ("p", "t", "a"):
        existing[key].write_text("", encoding="utf-8")
    blastx = tmp_path / "blastx.tsv"
    blastx.write_text("".join(
        f"Felid_ORSE{index:06d}\tOR{index}\t80\t3\t0\t0\t1\t9\t1\t3\t1e-20\t50\t9\t3\tOR-Receptor\n"
        for index in range(1, 4)
    ), encoding="utf-8")
    complete_prefix = tmp_path / "Felid"
    summary = build_complete_catalog(
        "Felid", catalog, gff, existing["fnt"], existing["fprot"], existing["p"],
        existing["t"], existing["a"], Path(f"{prefix}.single_exon_recovery_calls.tsv"),
        Path(f"{prefix}.single_exon_pseudogene.nucl.fa"),
        Path(f"{prefix}.single_exon_truncated.nucl.fa"),
        Path(f"{prefix}.single_exon_ambiguous.nucl.fa"), blastx, complete_prefix,
    )
    assert summary["counts"]["existing_catalog_loci"] == 1
    assert summary["counts"]["new_catalog_loci"] == 3
    assert summary["counts"]["complete_catalog_loci"] == 4
    with Path(f"{complete_prefix}.complete_OR_catalog.tsv").open(
        "rt", encoding="utf-8", newline=""
    ) as handle:
        ids = [row["catalog_id"] for row in csv.DictReader(handle, delimiter="\t")]
    assert ids == [
        "Felid_ORC000001", "Felid_ORP000001", "Felid_ORT000001", "Felid_ORA000001"
    ]
