import csv
from pathlib import Path

from or_mining.exonerate_model_selection import ExonerateModel, GffFeature, Region
from or_mining.fasta import read_fasta
from or_mining.multiexon_classification import (
    StrictCandidate,
    _find_intact_orf,
    _transcript_path,
    classify_multiexon_models,
    reverse_complement,
)


def _block(candidate_id: str, query_id: str, coding: str, frameshifts: int, query_end: int) -> str:
    frame_attribute = f" ; frameshifts {frameshifts}" if frameshifts else ""
    return (
        "C4 Alignment:\n"
        f"vulgar: {query_id} 0 {query_end} . {candidate_id} 100 1100 + 900 M {query_end} 900\n"
        "# --- START OF GFF DUMP ---\n"
        f"{candidate_id}\texonerate:protein2genome:local\tgene\t101\t1100\t900\t+\t.\t"
        f"gene_id 1 ; sequence {query_id} ; identity 90 ; similarity 95\n"
        f"{candidate_id}\texonerate:protein2genome:local\tcds\t101\t500\t.\t+\t.\t\n"
        f"{candidate_id}\texonerate:protein2genome:local\texon\t101\t500\t.\t+\t.\t"
        f"identity 90 ; similarity 95{frame_attribute}\n"
        f"{candidate_id}\texonerate:protein2genome:local\tintron\t501\t600\t.\t+\t.\tintron_id 1\n"
        f"{candidate_id}\texonerate:protein2genome:local\tcds\t601\t1100\t.\t+\t.\t\n"
        f"{candidate_id}\texonerate:protein2genome:local\texon\t601\t1100\t.\t+\t.\t"
        "identity 90 ; similarity 95\n"
        "# --- END OF GFF DUMP ---\n#\n"
        + coding[:400]
        + "\n"
        + coding[400:]
        + "\n# FELIDAE_OR_CHUNK_STATUS completed\n"
    )


def test_multiexon_classification_partitions_models_and_restores_coordinates(
    tmp_path: Path,
) -> None:
    intact = "ATG" + "GCT" * 298 + "TAA"
    disrupted = "ATG" + "GCT" * 98 + "TAA" + "GCT" * 200
    partial = "GCT" * 300
    assert len(intact) == len(disrupted) == len(partial) == 900

    region_sequences = {}
    for region_id, coding in (("r1", intact), ("r2", disrupted), ("r3", partial)):
        region_sequences[region_id] = "C" * 100 + coding[:400] + "A" * 100 + coding[400:] + "C" * 100

    regions_fasta = tmp_path / "regions.fa"
    regions_fasta.write_text(
        "".join(f">{region_id}\n{sequence}\n" for region_id, sequence in region_sequences.items()),
        encoding="utf-8",
    )
    regions_tsv = tmp_path / "regions.tsv"
    regions_tsv.write_text(
        "region_id\tsseqid\texpanded_start_1\texpanded_end_1\tslen\n"
        "r1\tchr1\t1001\t2200\t10000\n"
        "r2\tchr1\t3001\t4200\t10000\n"
        "r3\tchr1\t5001\t6200\t10000\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "strict.tsv"
    manifest.write_text(
        "candidate_id\tregion_id\tquery_id\tsseqid\tstrand\trefinement_start_1\t"
        "refinement_end_1\ttarget_bases\n"
        "c1\tr1\tq1\tchr1\t+\t1001\t2200\t1200\n"
        "c2\tr2\tq2\tchr1\t+\t3001\t4200\t1200\n"
        "c3\tr3\tq3\tchr1\t+\t5001\t6200\t1200\n",
        encoding="utf-8",
    )
    raw = tmp_path / "strict.raw.txt"
    raw.write_text(
        _block("c1", "q1", intact, 0, 299)
        + _block("c2", "q2", disrupted, 1, 299)
        + _block("c3", "q3", partial, 0, 150),
        encoding="utf-8",
    )
    reference = tmp_path / "reference.prot.fa"
    reference.write_text(
        ">q1\n" + "M" * 299 + "\n>q2\n" + "M" * 299 + "\n>q3\n" + "M" * 299 + "\n",
        encoding="utf-8",
    )
    species = tmp_path / "species.prot.fa"
    species.write_text(">species_query\nM\n", encoding="utf-8")

    prefix = tmp_path / "sample"
    summary = classify_multiexon_models(
        raw,
        manifest,
        regions_tsv,
        regions_fasta,
        reference,
        species,
        prefix,
    )
    assert summary["counts"] == {
        "input_candidates": 3,
        "classified_candidates": 3,
        "functional": 1,
        "pseudogene": 1,
        "truncated": 1,
        "ambiguous": 0,
        "output_exons": 6,
    }

    calls_path = tmp_path / "sample.multiexon_calls.tsv"
    with calls_path.open(encoding="utf-8", newline="") as handle:
        calls = {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert calls["c1"]["classification"] == "functional"
    assert calls["c1"]["genomic_start_1"] == "1101"
    assert calls["c1"]["genomic_end_1"] == "2100"
    assert calls["c1"]["orf_exons"] == "2"
    assert calls["c2"]["classification"] == "pseudogene"
    assert calls["c3"]["classification"] == "truncated"

    functional_nt = list(read_fasta(tmp_path / "sample.multiexon_functional.nucl.fa"))
    functional_protein = list(read_fasta(tmp_path / "sample.multiexon_functional.prot.fa"))
    assert len(functional_nt) == len(functional_protein) == 1
    assert len(functional_nt[0][1]) == 900
    assert functional_nt[0][1].startswith("ATG") and functional_nt[0][1].endswith("TAA")
    assert len(functional_protein[0][1]) == 299

    gff = (tmp_path / "sample.multiexon_calls.gff3").read_text(encoding="utf-8")
    assert "chr1\tfelidae_or\tCDS\t1101\t1500\t.\t+\t0\tParent=c1" in gff
    assert "chr1\tfelidae_or\tCDS\t1601\t2100\t.\t+\t2\tParent=c1" in gff


def test_minus_strand_transcript_path_is_reverse_complemented_in_transcription_order() -> None:
    intact = "ATG" + "GCT" * 298 + "TAA"
    high_exon = reverse_complement(intact[:500])
    low_exon = reverse_complement(intact[500:])
    target = "C" * 100 + low_exon + "A" * 100 + high_exon + "C" * 100
    region = Region("cminus", "chr2", 1001, 2200, 10000)
    gene = GffFeature("gene", 101, 1100, "-", "identity 90 ; similarity 95")
    exons = (
        GffFeature("exon", 101, 500, "-", "identity 90 ; similarity 95"),
        GffFeature("exon", 601, 1100, "-", "identity 90 ; similarity 95"),
    )
    model = ExonerateModel(
        block_index=1,
        query_id="qminus",
        query_start_0=0,
        query_end_0=299,
        target_id="cminus",
        target_start_0=1100,
        target_end_0=100,
        strand="-",
        score=900,
        gene=gene,
        exons=exons,
        introns=(GffFeature("intron", 501, 600, "-", "intron_id 1"),),
        target_coding_sequence=intact,
        identity=90,
        similarity=95,
        frameshifts=0,
        frameshift_bases=0,
        region=region,
    )
    candidate = StrictCandidate(
        "cminus", "rminus", "qminus", "chr2", "-", 1001, 2200, target, 10000
    )
    model_sequence, model_mapping, aligned = _transcript_path(
        model, candidate, include_outer_flanks=False
    )
    assert model_sequence == intact
    assert model_mapping[0] == 2100
    assert model_mapping[-1] == 1101

    transcript, mapping, aligned = _transcript_path(model, candidate, include_outer_flanks=True)
    orf = _find_intact_orf(transcript, mapping, aligned, 299, 750, 0.80, 2)
    assert orf is not None
    assert orf.nucleotide == intact
    assert orf.exon_intervals == ((1601, 2100), (1101, 1500))
