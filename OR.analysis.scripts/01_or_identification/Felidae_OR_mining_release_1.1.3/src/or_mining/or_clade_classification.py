"""NC-compatible evolutionary-clade assignment for a completed OR catalog."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any
import csv
import json
import re

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from or_mining.fasta import read_fasta


OR_CLADES = (
    "alpha", "beta", "gamma", "delta", "epsilon",
    "eta", "zeta", "kappa", "lambda", "theta",
)
CLASS_I_CLADES = {"alpha", "beta", "delta", "epsilon", "zeta"}
CLASS_II_CLADES = {"gamma"}
BLAST_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
)


def infer_reference_clade(identifier: str) -> str | None:
    """Infer a published OR clade from an NC reference sequence identifier."""

    lowered = identifier.lower()
    alternatives = "|".join(OR_CLADES)
    patterns = (
        rf"^({alternatives})---",
        rf"---({alternatives})-",
        rf"-(kappa|lambda|theta)(?:\d|[._-])",
    )
    for pattern in patterns:
        match = re.search(pattern, lowered)
        if match:
            return match.group(1)
    return None


def _read_ids(path: Path) -> set[str]:
    identifiers = {line.strip() for line in path.read_text(encoding="utf-8").splitlines()
                   if line.strip()}
    if not identifiers:
        raise ValueError(f"reference ID file is empty: {path}")
    return identifiers


def _fasta_ids(path: Path) -> list[str]:
    identifiers = [header.split()[0] for header, _ in read_fasta(path)]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"duplicate FASTA identifiers in {path}")
    return identifiers


def _template_labels(
    template_alignment: Path,
    type1_ids_file: Path,
    type2_ids_file: Path,
) -> tuple[dict[str, str], dict[str, str], set[str]]:
    template_ids = set(_fasta_ids(template_alignment))
    type1 = _read_ids(type1_ids_file)
    type2 = _read_ids(type2_ids_file)
    if type1 & type2:
        raise ValueError("Type I and Type II reference sets overlap")
    known = type1 | type2
    if not known <= template_ids:
        raise ValueError("known OR IDs are absent from the template alignment")

    clade_by_id: dict[str, str] = {}
    type_by_id: dict[str, str] = {}
    for identifier in sorted(known):
        clade = infer_reference_clade(identifier)
        if clade is None:
            raise ValueError(f"cannot infer OR clade from template ID {identifier!r}")
        clade_by_id[identifier] = clade
        type_by_id[identifier] = (
            "Known_OR_genes_TYPEI" if identifier in type1 else "Known_OR_genes_TYPEII"
        )

    nc_reference_set: dict[str, str] = {}
    for clade in sorted(set(clade_by_id.values())):
        observed = {type_by_id[identifier] for identifier, value in clade_by_id.items()
                    if value == clade}
        if len(observed) != 1:
            raise ValueError(f"template clade {clade} spans multiple OR types: {observed}")
        nc_reference_set[clade] = next(iter(observed))
    return clade_by_id, nc_reference_set, template_ids - known


def canonical_or_class(clade: str) -> str:
    """Return the classical vertebrate OR Class I/II designation.

    This is intentionally distinct from the NC tree-filter template files
    named ``Known_OR_genes_TYPEI/TYPEII``.  Those files encode two ancestral
    template sets and must not be interpreted as mammalian Class I/II labels.
    """

    if clade in CLASS_I_CLADES:
        return "class_I"
    if clade in CLASS_II_CLADES:
        return "class_II"
    if clade == "ambiguous":
        return "ambiguous"
    return "not_applicable"


def _normalize_support(value: float | None) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number / 100.0 if number > 1.0 else number


def _descendant_names(clade: Clade) -> set[str]:
    return {tip.name for tip in clade.get_terminals() if tip.name}


def _load_tree_assignments(
    tree_path: Path,
    functional_ids: set[str],
    clade_by_reference: dict[str, str],
) -> tuple[dict[str, dict[str, Any]], Tree, dict[str, Any]]:
    tree = Phylo.read(str(tree_path), "newick")
    leaf_names = [tip.name for tip in tree.get_terminals()]
    if any(name is None for name in leaf_names) or len(leaf_names) != len(set(leaf_names)):
        raise ValueError("tree has unnamed or duplicate leaves")
    leaf_set = set(leaf_names)
    required = functional_ids | set(clade_by_reference)
    if not required <= leaf_set:
        raise ValueError(f"tree is missing {len(required - leaf_set)} required leaves")

    # Clade membership is an unrooted property.  Testing only descendant MRCA
    # clades after rerooting can falsely split a valid group when the root lies
    # on its boundary (observed for gamma in the NC template).  Every tree edge
    # defines two partitions; a reference clade is monophyletic if either side
    # contains exactly its labelled references, while candidate leaves are
    # allowed within that partition.
    terminals = {tip.name: tip for tip in tree.get_terminals() if tip.name}
    all_leaves = set(terminals)
    reference_ids = set(clade_by_reference)
    nodes = list(tree.find_clades(order="postorder"))
    descendant_by_node_id: dict[int, set[str]] = {}
    for node in nodes:
        if node.is_terminal():
            descendant_by_node_id[id(node)] = {node.name} if node.name else set()
        else:
            descendants: set[str] = set()
            for child in node.clades:
                descendants.update(descendant_by_node_id[id(child)])
            descendant_by_node_id[id(node)] = descendants
    clade_references: dict[str, list[str]] = {}
    for identifier, clade in clade_by_reference.items():
        clade_references.setdefault(clade, []).append(identifier)

    monophyly: dict[str, bool] = {}
    members: dict[str, set[str]] = {}
    boundary_support: dict[str, float | None] = {}
    boundary_orientation: dict[str, str] = {}
    partition_sizes: dict[str, int] = {}
    for clade, references in clade_references.items():
        target = set(references)
        options: list[tuple[int, float, set[str], str, float | None]] = []
        for node in nodes:
            if node is tree.root:
                continue
            descendant = descendant_by_node_id[id(node)]
            for side, orientation in (
                (descendant, "descendant"),
                (all_leaves - descendant, "complement"),
            ):
                if side & reference_ids == target:
                    support = _normalize_support(node.confidence)
                    options.append((
                        len(side), -1.0 if support is None else support,
                        side, orientation, support,
                    ))
        monophyly[clade] = bool(options)
        if options:
            # The largest pure-reference partition includes lineage-specific
            # candidate expansions up to the first foreign reference boundary.
            _, _, selected, orientation, support = max(
                options, key=lambda item: (item[0], item[1])
            )
            members[clade] = selected & functional_ids
            boundary_support[clade] = support
            boundary_orientation[clade] = orientation
            partition_sizes[clade] = len(selected)
        else:
            members[clade] = set()
            boundary_support[clade] = None
            boundary_orientation[clade] = "not_resolved"
            partition_sizes[clade] = 0

    assignments: dict[str, dict[str, Any]] = {}
    adjacency: dict[int, list[tuple[int, float]]] = {id(node): [] for node in nodes}
    for parent in nodes:
        for child in parent.clades:
            length = float(child.branch_length or 0.0)
            adjacency[id(parent)].append((id(child), length))
            adjacency[id(child)].append((id(parent), length))
    candidate_by_node = {id(terminals[name]): name for name in functional_ids}
    distances_by_candidate: dict[str, list[tuple[float, str, str]]] = {
        name: [] for name in functional_ids
    }
    for reference_id, reference_clade in clade_by_reference.items():
        stack = [(id(terminals[reference_id]), -1, 0.0)]
        while stack:
            node_id, parent_id, distance = stack.pop()
            candidate_name = candidate_by_node.get(node_id)
            if candidate_name is not None:
                distances_by_candidate[candidate_name].append(
                    (distance, reference_id, reference_clade)
                )
            for neighbour_id, length in adjacency[node_id]:
                if neighbour_id != parent_id:
                    stack.append((neighbour_id, node_id, distance + length))
    for identifier in sorted(functional_ids):
        containing = [clade for clade, values in members.items() if identifier in values]
        ranked = sorted(distances_by_candidate[identifier])
        nearest_three = [item[2] for item in ranked[:3]]
        nearest_consensus = (
            nearest_three[0]
            if nearest_three and len(set(nearest_three)) == 1
            else "ambiguous"
        )
        if len(containing) == 1 and monophyly[containing[0]]:
            tree_clade = containing[0]
            method = "unrooted_reference_edge_partition"
            support = boundary_support[tree_clade]
        else:
            tree_clade = ranked[0][2]
            method = "nearest_reference_fallback"
            support = None
        assignments[identifier] = {
            "tree_clade": tree_clade,
            "tree_method": method,
            "tree_boundary_support": support,
            "nearest_reference_id": ranked[0][1],
            "nearest_reference_distance": ranked[0][0],
            "nearest_three_clades": ";".join(nearest_three),
            "nearest_three_consensus": nearest_consensus,
        }
    # Keep the original FastTree root only for reproducible serialization; all
    # classifications above are explicitly root-invariant.
    tree.rooted = False
    tree_audit = {
        "monophyly": monophyly,
        "boundary_support": boundary_support,
        "boundary_orientation": boundary_orientation,
        "partition_sizes": partition_sizes,
    }
    return assignments, tree, tree_audit


def _read_blast(path: Path) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            values = raw.rstrip("\n").split("\t")
            if len(values) != len(BLAST_FIELDS):
                raise ValueError(f"{path}:{line_number}: expected {len(BLAST_FIELDS)} columns")
            row = dict(zip(BLAST_FIELDS, values, strict=True))
            hits.setdefault(row["qseqid"], []).append(row)
    for values in hits.values():
        values.sort(key=lambda row: (-float(row["bitscore"]), float(row["evalue"]),
                                     -float(row["pident"]), row["sseqid"]))
    return hits


def _blast_evidence(
    hits: list[dict[str, str]],
    *,
    nucleotide_query: bool,
) -> dict[str, Any]:
    if not hits:
        return {
            "blast_clade": "unassigned", "blast_status": "no_hit",
            "blast_subject": "", "blast_pident": "", "blast_evalue": "",
            "blast_bitscore": "", "blast_query_coverage": "",
        }
    parsed = [(row, infer_reference_clade(row["sseqid"])) for row in hits]
    parsed = [(row, clade) for row, clade in parsed if clade]
    if not parsed:
        return {
            "blast_clade": "unassigned", "blast_status": "unlabelled_hits",
            "blast_subject": hits[0]["sseqid"], "blast_pident": hits[0]["pident"],
            "blast_evalue": hits[0]["evalue"], "blast_bitscore": hits[0]["bitscore"],
            "blast_query_coverage": "",
        }
    best, best_clade = parsed[0]
    best_score = float(best["bitscore"])
    near_best = {clade for row, clade in parsed
                 if float(row["bitscore"]) >= best_score * 0.95}
    status = "consistent_top_hits" if near_best == {best_clade} else "cross_clade_top_hits"
    coverage = (
        int(best["length"]) * (3 if nucleotide_query else 1) / int(best["qlen"])
    )
    return {
        "blast_clade": best_clade,
        "blast_status": status,
        "blast_subject": best["sseqid"],
        "blast_pident": best["pident"],
        "blast_evalue": best["evalue"],
        "blast_bitscore": best["bitscore"],
        "blast_query_coverage": f"{min(1.0, coverage):.6f}",
    }


def classify_or_clades(
    species: str,
    tree_path: Path,
    template_alignment: Path,
    type1_ids: Path,
    type2_ids: Path,
    catalog_path: Path,
    functional_protein: Path,
    functional_blastp: Path,
    pseudogene_blastx: Path,
    truncated_blastx: Path,
    ambiguous_blastx: Path,
    output_prefix: Path,
    minimum_tree_support: float = 0.80,
) -> dict[str, Any]:
    """Assign published evolutionary clades while keeping tree/BLAST audits separate."""

    clade_by_reference, nc_reference_set, outgroups = _template_labels(
        template_alignment, type1_ids, type2_ids
    )
    functional_ids = set(_fasta_ids(functional_protein))
    with catalog_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("complete OR catalog has no header")
        catalog_fields = list(reader.fieldnames)
        catalog_rows = list(reader)
    required_fields = {"catalog_id", "classification", "sseqid", "start_1", "end_1"}
    if not required_fields <= set(catalog_fields):
        raise ValueError(f"complete catalog lacks fields: {sorted(required_fields)}")
    catalog_ids = [row["catalog_id"] for row in catalog_rows]
    if len(catalog_ids) != len(set(catalog_ids)):
        raise ValueError("complete catalog contains duplicate IDs")
    expected_functional = {
        row["catalog_id"] for row in catalog_rows if row["classification"] == "functional"
    }
    if functional_ids != expected_functional:
        raise ValueError("functional protein FASTA and catalog IDs differ")

    tree_assignments, audited_tree, tree_audit = _load_tree_assignments(
        tree_path, functional_ids, clade_by_reference
    )
    monophyly = tree_audit["monophyly"]
    blast_by_class = {
        "functional": _read_blast(functional_blastp),
        "pseudogene": _read_blast(pseudogene_blastx),
        "truncated": _read_blast(truncated_blastx),
        "ambiguous": _read_blast(ambiguous_blastx),
    }

    decisions: list[dict[str, str]] = []
    final_counts: Counter[str] = Counter()
    confidence_counts: Counter[str] = Counter()
    for row in catalog_rows:
        identifier = row["catalog_id"]
        classification = row["classification"]
        blast = _blast_evidence(
            blast_by_class[classification].get(identifier, []),
            nucleotide_query=classification != "functional",
        )
        if classification == "functional":
            tree = tree_assignments[identifier]
            tree_clade = tree["tree_clade"]
            blast_clade = blast["blast_clade"]
            support = tree["tree_boundary_support"]
            topology_strong = (
                tree["tree_method"] == "unrooted_reference_edge_partition"
                and tree["nearest_three_consensus"] == tree_clade
                and (support is None or support >= minimum_tree_support)
            )
            if blast_clade in {tree_clade, "unassigned"}:
                final_clade = tree_clade
                audit_status = "tree_blast_agree" if blast_clade == tree_clade else "tree_only"
                confidence = "high" if topology_strong and blast_clade == tree_clade else "medium"
            elif topology_strong:
                final_clade = tree_clade
                audit_status = "blast_conflict_tree_retained"
                confidence = "medium"
            else:
                final_clade = "ambiguous"
                audit_status = "tree_blast_conflict_unresolved"
                confidence = "low"
            method = "phylogeny_primary_blastp_audit"
        else:
            tree = {
                "tree_clade": "not_applicable", "tree_method": "not_applicable",
                "tree_boundary_support": None, "nearest_reference_id": "",
                "nearest_reference_distance": "", "nearest_three_clades": "",
                "nearest_three_consensus": "not_applicable",
            }
            final_clade = blast["blast_clade"]
            method = "best_blastx_propagation"
            if final_clade == "unassigned":
                confidence = "unassigned"
                audit_status = blast["blast_status"]
            else:
                pident = float(blast["blast_pident"])
                coverage = float(blast["blast_query_coverage"])
                strong = (
                    blast["blast_status"] == "consistent_top_hits"
                    and float(blast["blast_evalue"]) <= 1e-20
                    and pident >= 40.0 and coverage >= 0.50
                )
                confidence = "high" if strong else "low"
                audit_status = blast["blast_status"]

        or_class = canonical_or_class(final_clade)
        decision = {
            "catalog_id": identifier,
            "classification": classification,
            "or_clade": final_clade,
            "or_class": or_class,
            "assignment_method": method,
            "assignment_confidence": confidence,
            "audit_status": audit_status,
            "tree_clade": str(tree["tree_clade"]),
            "tree_method": str(tree["tree_method"]),
            "tree_boundary_support": (
                "" if tree["tree_boundary_support"] is None
                else f"{float(tree['tree_boundary_support']):.6f}"
            ),
            "nearest_reference_id": str(tree["nearest_reference_id"]),
            "nearest_reference_distance": (
                "" if tree["nearest_reference_distance"] == ""
                else f"{float(tree['nearest_reference_distance']):.8f}"
            ),
            "nearest_three_clades": str(tree["nearest_three_clades"]),
            **{key: str(value) for key, value in blast.items()},
        }
        decisions.append(decision)
        final_counts[final_clade] += 1
        confidence_counts[confidence] += 1

    decision_fields = list(decisions[0]) if decisions else []
    assignments_path = Path(f"{output_prefix}.OR_clade_assignments.tsv")
    annotated_path = Path(f"{output_prefix}.OR_catalog_with_clades.tsv")
    counts_path = Path(f"{output_prefix}.OR_clade_counts.tsv")
    audited_tree_path = Path(f"{output_prefix}.OR_clades.audited.unrooted.treefile")
    summary_path = Path(f"{output_prefix}.OR_clade_classification.summary.json")
    with assignments_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=decision_fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(decisions)
    by_id = {row["catalog_id"]: row for row in decisions}
    added = ("or_clade", "or_class", "clade_assignment_method", "clade_assignment_confidence")
    with annotated_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=catalog_fields + list(added), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in catalog_rows:
            decision = by_id[row["catalog_id"]]
            writer.writerow({
                **row,
                "or_clade": decision["or_clade"],
                "or_class": decision["or_class"],
                "clade_assignment_method": decision["assignment_method"],
                "clade_assignment_confidence": decision["assignment_confidence"],
            })
    with counts_path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("classification\tor_clade\tor_class\tcount\n")
        grouped = Counter(
            (row["classification"], row["or_clade"], row["or_class"]) for row in decisions
        )
        for key, count in sorted(grouped.items()):
            handle.write("\t".join((*key, str(count))) + "\n")
    Phylo.write(audited_tree, str(audited_tree_path), "newick")

    observed_functional_clades = {
        row["or_clade"] for row in decisions
        if row["classification"] == "functional" and row["or_clade"] in monophyly
    }
    nonmonophyletic = sorted(clade for clade, state in monophyly.items() if not state)
    observed_nonmonophyletic = sorted(observed_functional_clades & set(nonmonophyletic))

    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "pass" if not nonmonophyletic else "pass_with_warnings",
        "species": species,
        "method": {
            "functional": "FastTree reference-clade phylogeny with BLASTP audit",
            "incomplete": "NC-compatible best BLASTX clade propagation",
        },
        "parameters": {"minimum_tree_support": minimum_tree_support},
        "reference": {
            "or_reference_sequences": len(clade_by_reference),
            "outgroup_sequences": len(outgroups),
            "class_I_clades": sorted(CLASS_I_CLADES),
            "class_II_clades": sorted(CLASS_II_CLADES),
            "not_applicable_to_class_I_II": ["eta", "kappa", "lambda", "theta"],
            "nc_filter_reference_set_by_clade": dict(sorted(nc_reference_set.items())),
            "reference_clade_monophyly": dict(sorted(monophyly.items())),
            "reference_clade_boundary_support": dict(
                sorted(tree_audit["boundary_support"].items())
            ),
            "reference_clade_boundary_orientation": dict(
                sorted(tree_audit["boundary_orientation"].items())
            ),
            "reference_clade_partition_sizes": dict(
                sorted(tree_audit["partition_sizes"].items())
            ),
            "nonmonophyletic_reference_clades": nonmonophyletic,
        },
        "counts": {
            "catalog_loci": len(catalog_rows),
            "functional_tree_assignments": len(functional_ids),
            "by_clade": dict(sorted(final_counts.items())),
            "by_confidence": dict(sorted(confidence_counts.items())),
        },
        "checks": {
            "all_catalog_loci_reported": len(decisions) == len(catalog_rows),
            "all_reference_clades_monophyletic": all(monophyly.values()),
            "all_observed_functional_clades_monophyletic": not observed_nonmonophyletic,
            "observed_nonmonophyletic_functional_clades": observed_nonmonophyletic,
            "functional_ids_match_catalog": functional_ids == expected_functional,
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary
