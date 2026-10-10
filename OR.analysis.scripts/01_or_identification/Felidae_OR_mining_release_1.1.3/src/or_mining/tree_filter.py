"""Audited phylogenetic filtering of OR candidates."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any
import csv
import json

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from or_mining.fasta import read_fasta, write_fasta_record


def _read_ids(path: Path) -> list[str]:
    identifiers = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    identifiers = [identifier for identifier in identifiers if identifier]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"duplicate identifiers in {path}")
    return identifiers


def _fasta_records(path: Path) -> tuple[list[str], dict[str, tuple[str, str]]]:
    order: list[str] = []
    records: dict[str, tuple[str, str]] = {}
    for header, sequence in read_fasta(path):
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
        order.append(identifier)
        records[identifier] = (header, sequence)
    return order, records


def _descendant_names(clade: Clade) -> set[str]:
    return {tip.name for tip in clade.get_terminals() if tip.name}


def filter_or_clade(
    tree_file: Path,
    template_alignment: Path,
    type1_ids_file: Path,
    type2_ids_file: Path,
    candidate_proteins: Path,
    candidate_nucleotides: Path,
    candidate_coordinates: Path,
    output_prefix: Path,
) -> dict[str, Any]:
    """Apply the published NC MRCA rule and report an explicit-outgroup audit."""

    type1 = _read_ids(type1_ids_file)
    type2 = _read_ids(type2_ids_file)
    known = type1 + type2
    if len(known) != len(set(known)):
        raise ValueError("Type I and Type II reference ID sets overlap")

    template_order, _ = _fasta_records(template_alignment)
    template_ids = set(template_order)
    known_set = set(known)
    missing_known_template = sorted(known_set - template_ids)
    if missing_known_template:
        raise ValueError(
            f"{len(missing_known_template)} known OR IDs are absent from the template alignment"
        )
    outgroups = sorted(template_ids - known_set)
    if not outgroups:
        raise ValueError("template alignment contains no explicit non-OR outgroup sequences")

    candidate_order, proteins = _fasta_records(candidate_proteins)
    _, nucleotides = _fasta_records(candidate_nucleotides)
    candidate_set = set(candidate_order)
    if candidate_set != set(nucleotides):
        raise ValueError("candidate protein and nucleotide FASTA identifiers differ")
    if candidate_set & template_ids:
        raise ValueError("candidate identifiers overlap template identifiers")

    tree = Phylo.read(str(tree_file), "newick")
    leaf_names = [tip.name for tip in tree.get_terminals()]
    if any(name is None for name in leaf_names):
        raise ValueError("tree contains unnamed terminal nodes")
    if len(leaf_names) != len(set(leaf_names)):
        raise ValueError("tree contains duplicate terminal identifiers")
    expected_leaves = template_ids | candidate_set
    if set(leaf_names) != expected_leaves:
        missing = sorted(expected_leaves - set(leaf_names))
        unexpected = sorted(set(leaf_names) - expected_leaves)
        raise ValueError(
            f"tree leaf mismatch: missing={missing[:5]} unexpected={unexpected[:5]}"
        )
    legacy_tree = deepcopy(tree)
    initial_or_mrca = legacy_tree.common_ancestor(known)
    legacy_tree.root_with_outgroup(initial_or_mrca)
    legacy_or_mrca = legacy_tree.common_ancestor(known)
    legacy_retained = _descendant_names(legacy_or_mrca) & candidate_set

    audited_tree = deepcopy(tree)
    audited_terminals_by_name = {
        tip.name: tip for tip in audited_tree.get_terminals() if tip.name
    }
    outgroup_terminals = [audited_terminals_by_name[seq_id] for seq_id in outgroups]
    outgroup_clade = audited_tree.is_monophyletic(outgroup_terminals)
    if outgroup_clade is False:
        audit_status = "not_applicable_non_monophyletic_outgroups"
        audited_retained: set[str] | None = None
        legacy_only: list[str] | None = None
        explicit_outgroup_only: list[str] | None = None
        algorithms_agree: bool | None = None
    else:
        audited_tree.root_with_outgroup(outgroup_clade)
        audited_or_mrca = audited_tree.common_ancestor(known)
        audited_retained = _descendant_names(audited_or_mrca) & candidate_set
        legacy_only = sorted(legacy_retained - audited_retained)
        explicit_outgroup_only = sorted(audited_retained - legacy_retained)
        algorithms_agree = not legacy_only and not explicit_outgroup_only
        audit_status = "completed"

    # The primary outputs reproduce Tree_parser_2.R from the NC workflow. The
    # explicit-outgroup result is a sensitivity audit and is never substituted
    # silently for the published rule.
    retained_ids = [identifier for identifier in candidate_order if identifier in legacy_retained]
    rejected_ids = [identifier for identifier in candidate_order if identifier not in legacy_retained]
    retained_set = set(retained_ids)

    retained_protein_path = Path(f"{output_prefix}.tree_filtered.prot.fa")
    retained_nucleotide_path = Path(f"{output_prefix}.tree_filtered.nt.fa")
    retained_coordinates_path = Path(f"{output_prefix}.tree_filtered.tsv")
    decisions_path = Path(f"{output_prefix}.tree_filter.decisions.tsv")
    retained_ids_path = Path(f"{output_prefix}.tree_filtered.ids.txt")
    rooted_tree_path = Path(f"{output_prefix}.nc_legacy.rooted.treefile")
    summary_path = Path(f"{output_prefix}.tree_filter.summary.json")

    with (
        retained_protein_path.open("wt", encoding="utf-8", newline="\n") as protein_out,
        retained_nucleotide_path.open("wt", encoding="utf-8", newline="\n") as nucleotide_out,
        retained_ids_path.open("wt", encoding="utf-8", newline="\n") as ids_out,
        decisions_path.open("wt", encoding="utf-8", newline="\n") as decisions,
    ):
        decisions.write(
            "orf_id\tnc_legacy_decision\texplicit_outgroup_decision\taudit_status\n"
        )
        for identifier in candidate_order:
            if identifier in retained_set:
                write_fasta_record(protein_out, proteins[identifier][0], proteins[identifier][1])
                write_fasta_record(
                    nucleotide_out, nucleotides[identifier][0], nucleotides[identifier][1]
                )
                ids_out.write(f"{identifier}\n")
                nc_decision = "retain"
            else:
                nc_decision = "reject"
            if audited_retained is None:
                explicit_decision = "not_applicable"
                candidate_audit_status = "not_applicable"
            else:
                explicit_decision = (
                    "retain" if identifier in audited_retained else "reject"
                )
                candidate_audit_status = (
                    "agree" if nc_decision == explicit_decision else "disagree"
                )
            decisions.write(
                f"{identifier}\t{nc_decision}\t{explicit_decision}"
                f"\t{candidate_audit_status}\n"
            )

    with candidate_coordinates.open("rt", encoding="utf-8", newline="") as source, retained_coordinates_path.open(
        "wt", encoding="utf-8", newline="\n"
    ) as output:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("candidate coordinate table has no header")
        writer = csv.DictWriter(
            output, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        observed: set[str] = set()
        for row in reader:
            identifier = row["orf_id"]
            observed.add(identifier)
            if identifier in retained_set:
                writer.writerow(row)
        if observed != candidate_set:
            raise ValueError("candidate coordinate table identifiers differ from FASTA identifiers")

    Phylo.write(legacy_tree, str(rooted_tree_path), "newick")
    summary: dict[str, Any] = {
        "schema_version": "1.1",
        "status": "pass",
        "primary_method": "NC Tree_parser_2.R-compatible known-OR MRCA rerooting",
        "audit_method": "known-OR MRCA after explicit non-OR outgroup rooting",
        "audit_status": audit_status,
        "algorithms_agree": algorithms_agree,
        "counts": {
            "tree_leaves": len(leaf_names),
            "known_type1_or": len(type1),
            "known_type2_or": len(type2),
            "template_outgroups": len(outgroups),
            "input_candidates": len(candidate_order),
            "retained_candidates": len(retained_ids),
            "rejected_candidates": len(rejected_ids),
            "explicit_outgroup_retained_candidates": (
                len(audited_retained) if audited_retained is not None else None
            ),
            "legacy_only_candidates": (
                len(legacy_only) if legacy_only is not None else None
            ),
            "explicit_outgroup_only_candidates": (
                len(explicit_outgroup_only)
                if explicit_outgroup_only is not None
                else None
            ),
        },
        "outgroup_ids": outgroups,
        "legacy_only_ids": legacy_only,
        "explicit_outgroup_only_ids": explicit_outgroup_only,
    }
    with summary_path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary
