"""Single-species genomic architecture summaries for a completed OR catalog."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any
import csv
import json


def _read_fai(path: Path) -> tuple[list[str], dict[str, int]]:
    order: list[str] = []
    lengths: dict[str, int] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed FAI row")
            if fields[0] in lengths:
                raise ValueError(f"duplicate FAI sequence {fields[0]!r}")
            order.append(fields[0])
            lengths[fields[0]] = int(fields[1])
    if not order:
        raise ValueError("genome FAI is empty")
    return order, lengths


def analyze_or_locus_architecture(
    species: str,
    catalog_path: Path,
    genome_fai: Path,
    output_prefix: Path,
    cluster_gap_bp: int = 1_000_000,
) -> dict[str, Any]:
    """Report coordinate facts and explicitly operational OR clusters."""

    if cluster_gap_bp < 0:
        raise ValueError("cluster gap must be non-negative")
    contig_order, contig_lengths = _read_fai(genome_fai)
    with catalog_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("annotated catalog has no header")
        fields = list(reader.fieldnames)
        rows = list(reader)
    required = {
        "catalog_id", "classification", "sseqid", "start_1", "end_1", "strand",
        "exon_count", "or_clade", "or_class",
    }
    if not required <= set(fields):
        raise ValueError(f"annotated catalog lacks fields: {sorted(required - set(fields))}")
    identifiers = [row["catalog_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("annotated catalog contains duplicate IDs")

    rank = {contig: index for index, contig in enumerate(contig_order)}
    loci_by_contig: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        contig = row["sseqid"]
        if contig not in contig_lengths:
            raise ValueError(f"catalog sequence {contig!r} is absent from genome FAI")
        start = int(row["start_1"])
        end = int(row["end_1"])
        if not (1 <= start <= end <= contig_lengths[contig]):
            raise ValueError(f"invalid coordinates for {row['catalog_id']}")
        loci_by_contig.setdefault(contig, []).append(row)
    for values in loci_by_contig.values():
        values.sort(key=lambda row: (int(row["start_1"]), int(row["end_1"]), row["catalog_id"]))

    cluster_counter = 0
    clusters: list[dict[str, Any]] = []
    locus_metrics: dict[str, dict[str, Any]] = {}
    for contig in sorted(loci_by_contig, key=lambda value: rank[value]):
        values = loci_by_contig[contig]
        gaps: list[int] = []
        for left, right in zip(values, values[1:]):
            gap = int(right["start_1"]) - int(left["end_1"]) - 1
            if gap < 0:
                raise ValueError(
                    f"catalog loci overlap: {left['catalog_id']}/{right['catalog_id']} ({-gap} bp)"
                )
            gaps.append(gap)

        components: list[list[dict[str, str]]] = []
        current: list[dict[str, str]] = []
        for index, row in enumerate(values):
            if current and gaps[index - 1] > cluster_gap_bp:
                components.append(current)
                current = []
            current.append(row)
        if current:
            components.append(current)

        for component in components:
            cluster_counter += 1
            cluster_id = f"{species}_ORCL{cluster_counter:06d}"
            start = int(component[0]["start_1"])
            end = max(int(row["end_1"]) for row in component)
            class_counts = Counter(row["classification"] for row in component)
            clade_counts = Counter(row["or_clade"] for row in component)
            clusters.append({
                "cluster_id": cluster_id,
                "sseqid": contig,
                "start_1": start,
                "end_1": end,
                "span_bp": end - start + 1,
                "member_count": len(component),
                "functional_count": class_counts["functional"],
                "pseudogene_count": class_counts["pseudogene"],
                "truncated_count": class_counts["truncated"],
                "ambiguous_count": class_counts["ambiguous"],
                "clade_count": len(clade_counts),
                "clade_composition": ";".join(
                    f"{clade}:{count}" for clade, count in sorted(clade_counts.items())
                ),
                "member_ids": ";".join(row["catalog_id"] for row in component),
            })
            for row in component:
                locus_metrics[row["catalog_id"]] = {
                    "operational_cluster_id": cluster_id,
                    "operational_cluster_size": len(component),
                }

        for index, row in enumerate(values):
            previous_gap = "" if index == 0 else gaps[index - 1]
            next_gap = "" if index == len(values) - 1 else gaps[index]
            numeric = [gap for gap in (previous_gap, next_gap) if gap != ""]
            locus_metrics[row["catalog_id"]].update({
                "previous_or_gap_bp": previous_gap,
                "next_or_gap_bp": next_gap,
                "nearest_or_gap_bp": min(numeric) if numeric else "",
                "left_contig_edge_bp": int(row["start_1"]) - 1,
                "right_contig_edge_bp": contig_lengths[contig] - int(row["end_1"]),
            })

    locus_path = Path(f"{output_prefix}.OR_locus_architecture.tsv")
    cluster_path = Path(f"{output_prefix}.OR_operational_clusters.tsv")
    contig_path = Path(f"{output_prefix}.OR_contig_counts.tsv")
    composition_path = Path(f"{output_prefix}.OR_clade_composition.tsv")
    summary_path = Path(f"{output_prefix}.OR_locus_architecture.summary.json")
    metric_fields = (
        "operational_cluster_id", "operational_cluster_size", "previous_or_gap_bp",
        "next_or_gap_bp", "nearest_or_gap_bp", "left_contig_edge_bp", "right_contig_edge_bp",
    )
    with locus_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields + list(metric_fields), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda value: (
            rank[value["sseqid"]], int(value["start_1"]), int(value["end_1"]), value["catalog_id"]
        )):
            writer.writerow({**row, **locus_metrics[row["catalog_id"]]})

    cluster_fields = list(clusters[0]) if clusters else []
    with cluster_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cluster_fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(clusters)

    with contig_path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "sseqid\tsequence_length_bp\ttotal_or\tfunctional\tpseudogene\ttruncated"
            "\tambiguous\toperational_clusters\tlargest_cluster\n"
        )
        for contig in sorted(loci_by_contig, key=lambda value: rank[value]):
            values = loci_by_contig[contig]
            counts = Counter(row["classification"] for row in values)
            local_clusters = [row for row in clusters if row["sseqid"] == contig]
            handle.write(
                f"{contig}\t{contig_lengths[contig]}\t{len(values)}\t{counts['functional']}\t"
                f"{counts['pseudogene']}\t{counts['truncated']}\t{counts['ambiguous']}\t"
                f"{len(local_clusters)}\t{max((row['member_count'] for row in local_clusters), default=0)}\n"
            )

    composition = Counter(
        (row["or_clade"], row["or_class"], row["classification"]) for row in rows
    )
    with composition_path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("or_clade\tor_class\tclassification\tcount\n")
        for key, count in sorted(composition.items()):
            handle.write("\t".join((*key, str(count))) + "\n")

    class_counts = Counter(row["classification"] for row in rows)
    multi_member = [cluster for cluster in clusters if cluster["member_count"] > 1]
    functional_plus_pseudo = class_counts["functional"] + class_counts["pseudogene"]
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass",
        "species": species,
        "method": "coordinate-ordered single-species OR architecture audit",
        "parameters": {
            "operational_cluster_rule": "adjacent OR loci separated by at most the threshold",
            "operational_cluster_gap_bp": cluster_gap_bp,
            "cluster_interpretation": "physical neighborhood, not an evolutionary family",
        },
        "counts": {
            "catalog_loci": len(rows),
            "sequences_with_or": len(loci_by_contig),
            "functional": class_counts["functional"],
            "pseudogene": class_counts["pseudogene"],
            "truncated": class_counts["truncated"],
            "ambiguous": class_counts["ambiguous"],
            "operational_clusters_including_singletons": len(clusters),
            "multi_locus_operational_clusters": len(multi_member),
            "loci_in_multi_locus_clusters": sum(cluster["member_count"] for cluster in multi_member),
            "largest_operational_cluster": max((cluster["member_count"] for cluster in clusters), default=0),
        },
        "fractions": {
            "pseudogene_of_catalog": class_counts["pseudogene"] / len(rows) if rows else 0.0,
            "pseudogene_of_functional_plus_pseudogene": (
                class_counts["pseudogene"] / functional_plus_pseudo if functional_plus_pseudo else 0.0
            ),
        },
        "checks": {
            "all_loci_assigned_to_one_operational_cluster": len(locus_metrics) == len(rows),
            "all_catalog_sequences_present_in_fai": all(row["sseqid"] in contig_lengths for row in rows),
            "no_overlapping_catalog_loci": True,
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary
