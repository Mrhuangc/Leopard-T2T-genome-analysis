#!/usr/bin/env python3
"""Reference-anchored OR subfamily assignment from 17 unrooted family trees."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


FAMILIES = [
    "OR1", "OR2", "OR3", "OR4", "OR5", "OR6", "OR7", "OR8", "OR9",
    "OR10", "OR11", "OR12", "OR13", "OR14", "OR51", "OR52", "OR56",
]


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    label: str = ""
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0

    def skip_space(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def token(self) -> str:
        self.skip_space()
        if self.i < len(self.text) and self.text[self.i] == "'":
            self.i += 1
            chars = []
            while self.i < len(self.text):
                if self.text[self.i] == "'":
                    if self.i + 1 < len(self.text) and self.text[self.i + 1] == "'":
                        chars.append("'")
                        self.i += 2
                        continue
                    self.i += 1
                    break
                chars.append(self.text[self.i])
                self.i += 1
            return "".join(chars)
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();":
            self.i += 1
        return self.text[start:self.i].strip()

    def branch_length(self) -> float:
        self.skip_space()
        if self.i >= len(self.text) or self.text[self.i] != ":":
            return 0.0
        self.i += 1
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();":
            self.i += 1
        raw = self.text[start:self.i].strip()
        return float(raw) if raw else 0.0

    def subtree(self) -> Node:
        self.skip_space()
        if self.text[self.i] == "(":
            self.i += 1
            children = [self.subtree()]
            while True:
                self.skip_space()
                if self.text[self.i] == ",":
                    self.i += 1
                    children.append(self.subtree())
                elif self.text[self.i] == ")":
                    self.i += 1
                    break
                else:
                    raise ValueError(f"unexpected Newick character at {self.i}: {self.text[self.i:self.i+20]}")
            label = self.token()
            node = Node(label=label, length=self.branch_length(), children=children)
            for child in children:
                child.parent = node
            return node
        name = self.token()
        if not name:
            raise ValueError(f"empty leaf name at Newick position {self.i}")
        return Node(name=name, length=self.branch_length())

    def parse(self) -> Node:
        root = self.subtree()
        self.skip_space()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self.skip_space()
        if self.i != len(self.text):
            raise ValueError(f"trailing Newick text: {self.text[self.i:self.i+50]}")
        return root


def iter_nodes(root: Node) -> Iterable[Node]:
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def support(label: str) -> tuple[float | None, float | None]:
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)\s*", label or "")
    if not match:
        return None, None
    return float(match.group(1)), float(match.group(2))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analyze_family(tree_path: Path, family: str, human_subfamily: dict[str, str]) -> list[dict[str, object]]:
    root = NewickParser(tree_path.read_text(encoding="utf-8")).parse()
    nodes = list(iter_nodes(root))
    leaves = [node for node in nodes if not node.children]
    by_name = {node.name: node for node in leaves}
    if len(by_name) != len(leaves):
        raise ValueError(f"duplicate tree leaf IDs in {family}")
    human_nodes = [node for node in leaves if node.name.startswith("HUMAN__")]
    query_nodes = [node for node in leaves if node.name.startswith("ORREP")]
    unknown_human = sorted(node.name for node in human_nodes if node.name not in human_subfamily)
    if unknown_human:
        raise ValueError(f"human tree IDs absent from metadata in {family}: {unknown_human[:10]}")
    if not human_nodes or not query_nodes:
        raise ValueError(f"family lacks human anchors or felid queries: {family}")

    adjacency: dict[Node, list[tuple[Node, float]]] = defaultdict(list)
    for node in nodes:
        for child in node.children:
            length = max(float(child.length or 0.0), 0.0)
            adjacency[node].append((child, length))
            adjacency[child].append((node, length))

    distances: dict[str, list[tuple[float, str]]] = defaultdict(list)
    query_set = set(query_nodes)
    for anchor in human_nodes:
        stack = [(anchor, None, 0.0)]
        while stack:
            node, parent, distance = stack.pop()
            if node in query_set:
                distances[node.name].append((distance, anchor.name))
            for neighbor, length in adjacency[node]:
                if neighbor is not parent:
                    stack.append((neighbor, node, distance + length))

    all_leaf_names = {node.name for node in leaves}
    descendant: dict[Node, set[str]] = {}

    def collect(node: Node) -> set[str]:
        if not node.children:
            result = {node.name}
        else:
            result = set().union(*(collect(child) for child in node.children))
        descendant[node] = result
        return result

    collect(root)
    pure: dict[str, dict[str, tuple[int, int, float | None, float | None]]] = defaultdict(dict)
    for node in nodes:
        if node is root:
            continue
        sh, uf = support(node.label)
        for side in (descendant[node], all_leaf_names - descendant[node]):
            anchors = [name for name in side if name.startswith("HUMAN__")]
            subfamilies = {human_subfamily[name] for name in anchors}
            if len(subfamilies) != 1:
                continue
            subfamily = next(iter(subfamilies))
            queries = [name for name in side if name.startswith("ORREP")]
            candidate = (len(side), len(anchors), sh, uf)
            for query in queries:
                previous = pure[query].get(subfamily)
                if previous is None or candidate[0] < previous[0]:
                    pure[query][subfamily] = candidate

    rows = []
    for query in sorted(node.name for node in query_nodes):
        refs = sorted(distances[query])
        if not refs:
            raise ValueError(f"no human-reference distance for {query} in {family}")
        by_subfamily: dict[str, tuple[float, str]] = {}
        for distance, ref in refs:
            subfamily = human_subfamily[ref]
            if subfamily not in by_subfamily or distance < by_subfamily[subfamily][0]:
                by_subfamily[subfamily] = (distance, ref)
        ranked = sorted((distance, subfamily, ref) for subfamily, (distance, ref) in by_subfamily.items())
        best_distance, best_subfamily, best_ref = ranked[0]
        second_distance = ranked[1][0] if len(ranked) > 1 else math.inf
        ratio = math.inf if best_distance <= 0 or math.isinf(second_distance) else second_distance / best_distance
        nearest = refs[: min(3, len(refs))]
        neighbor_counts = Counter(human_subfamily[ref] for _, ref in nearest)
        consensus_subfamily, consensus_count = sorted(
            neighbor_counts.items(), key=lambda item: (-item[1], item[0])
        )[0]
        pure_info = pure[query].get(best_subfamily)
        rows.append({
            "representative_id": query,
            "tree_family": family,
            "tree_nearest_reference_id": best_ref,
            "tree_nearest_subfamily": best_subfamily,
            "tree_nearest_distance": f"{best_distance:.10g}",
            "tree_second_subfamily_distance": "inf" if math.isinf(second_distance) else f"{second_distance:.10g}",
            "tree_distance_ratio": "inf" if math.isinf(ratio) else f"{ratio:.6f}",
            "nearest_3_consensus_subfamily": consensus_subfamily,
            "nearest_3_consensus_fraction": f"{consensus_count / len(nearest):.6f}",
            "reference_pure_split_present": "true" if pure_info else "false",
            "reference_pure_split_leaves": pure_info[0] if pure_info else "",
            "reference_pure_split_human_anchors": pure_info[1] if pure_info else "",
            "reference_pure_split_sh_alrt": "" if not pure_info or pure_info[2] is None else f"{pure_info[2]:g}",
            "reference_pure_split_ufboot": "" if not pure_info or pure_info[3] is None else f"{pure_info[3]:g}",
        })
    return rows


def as_float(value: object, default: float = math.nan) -> float:
    text = str(value or "").strip()
    if text.lower() == "inf":
        return math.inf
    try:
        return float(text)
    except ValueError:
        return default


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family-tree-root", required=True, type=Path)
    parser.add_argument("--candidate-assignments", required=True, type=Path)
    parser.add_argument("--human-metadata", required=True, type=Path)
    parser.add_argument("--identical-sequence-map", required=True, type=Path)
    parser.add_argument("--species-metadata", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--high-sh-alrt", type=float, default=80.0)
    parser.add_argument("--high-ufboot", type=float, default=95.0)
    parser.add_argument("--high-distance-ratio", type=float, default=1.10)
    parser.add_argument("--medium-distance-ratio", type=float, default=1.02)
    args = parser.parse_args()

    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    candidate_rows = read_tsv(args.candidate_assignments)
    human_rows = read_tsv(args.human_metadata)
    locus_rows = read_tsv(args.identical_sequence_map)
    species_rows = read_tsv(args.species_metadata)
    if len(candidate_rows) != 5338 or len({x["representative_id"] for x in candidate_rows}) != 5338:
        raise ValueError("candidate assignment table must contain 5,338 unique representatives")
    if len(human_rows) != 417 or len({x["reference_id"] for x in human_rows}) != 417:
        raise ValueError("human metadata must contain 417 unique references")
    if len(locus_rows) != 6419 or len({x["global_id"] for x in locus_rows}) != 6419:
        raise ValueError("identical-sequence map must contain 6,419 unique global loci")
    if len(species_rows) != 7 or len({x["analysis_unit"] for x in species_rows}) != 7:
        raise ValueError("species metadata must contain seven unique analysis units")

    human_subfamily = {row["reference_id"]: row["or_subfamily"] for row in human_rows}
    placement_rows = []
    for family in FAMILIES:
        tree_path = args.family_tree_root / "families" / family / "tree" / f"{family}.treefile"
        if not tree_path.is_file() or tree_path.stat().st_size == 0:
            raise ValueError(f"family tree missing or empty: {tree_path}")
        placement_rows.extend(analyze_family(tree_path, family, human_subfamily))
    placement_by_key = {(x["representative_id"], x["tree_family"]): x for x in placement_rows}
    if len(placement_rows) != 5957 or len(placement_by_key) != 5957:
        raise ValueError(f"expected 5,957 unique felid family-tree placements, observed {len(placement_rows)}")

    placement_fields = list(placement_rows[0])
    write_tsv(output / "all_candidate_tree_placements.tsv", placement_rows, placement_fields)

    final_rows = []
    for candidate in candidate_rows:
        representative = candidate["representative_id"]
        family = candidate.get("best_family", "").strip()
        placement = placement_by_key.get((representative, family)) if family else None
        status = candidate.get("candidate_status", "")
        final_subfamily = ""
        confidence = "unresolved"
        decision = "unresolved_no_qualified_reference_placement"
        if placement is not None and status not in {
            "not_assignable_with_human_reference", "no_qualified_reference_hit"
        }:
            tree_subfamily = str(placement["tree_nearest_subfamily"])
            blast_subfamily = candidate.get("best_subfamily", "").strip()
            agreement = tree_subfamily == blast_subfamily and bool(blast_subfamily)
            ratio = as_float(placement["tree_distance_ratio"])
            consensus = placement["nearest_3_consensus_subfamily"] == tree_subfamily
            consensus_fraction = as_float(placement["nearest_3_consensus_fraction"], 0.0)
            sh = as_float(placement["reference_pure_split_sh_alrt"])
            uf = as_float(placement["reference_pure_split_ufboot"])
            anchor_count = int(placement["reference_pure_split_human_anchors"] or 0)
            if not agreement:
                decision = "tree_blast_conflict"
                confidence = "unresolved"
            else:
                final_subfamily = tree_subfamily
                if (
                    ratio >= args.high_distance_ratio
                    and consensus and consensus_fraction == 1.0
                    and anchor_count >= 2
                    and not math.isnan(sh) and sh >= args.high_sh_alrt
                    and not math.isnan(uf) and uf >= args.high_ufboot
                ):
                    confidence = "high"
                    decision = "tree_blast_agree_strong_supported_split"
                elif ratio >= args.medium_distance_ratio and consensus and consensus_fraction >= 2 / 3:
                    confidence = "medium"
                    decision = "tree_blast_agree_neighborhood_supported"
                else:
                    confidence = "low"
                    decision = "tree_blast_agree_weak_separation_or_support"
        final_rows.append({
            **candidate,
            "final_family": family if placement is not None and decision != "unresolved_no_qualified_reference_placement" else "",
            "final_subfamily": final_subfamily,
            "tree_nearest_reference_id": placement["tree_nearest_reference_id"] if placement else "",
            "tree_nearest_subfamily": placement["tree_nearest_subfamily"] if placement else "",
            "tree_distance_ratio": placement["tree_distance_ratio"] if placement else "",
            "nearest_3_consensus_subfamily": placement["nearest_3_consensus_subfamily"] if placement else "",
            "nearest_3_consensus_fraction": placement["nearest_3_consensus_fraction"] if placement else "",
            "reference_pure_split_human_anchors": placement["reference_pure_split_human_anchors"] if placement else "",
            "reference_pure_split_sh_alrt": placement["reference_pure_split_sh_alrt"] if placement else "",
            "reference_pure_split_ufboot": placement["reference_pure_split_ufboot"] if placement else "",
            "assignment_confidence": confidence,
            "assignment_decision": decision,
        })

    final_fields = list(final_rows[0])
    write_tsv(output / "representative_OR_subfamily_assignments.tsv", final_rows, final_fields)
    final_by_rep = {row["representative_id"]: row for row in final_rows}
    if len(final_by_rep) != 5338:
        raise ValueError("final representative assignment table is not one row per representative")

    expanded_rows = []
    for locus in locus_rows:
        assignment = final_by_rep.get(locus["representative_id"])
        if assignment is None:
            raise ValueError(f"locus references unknown representative: {locus['representative_id']}")
        expanded_rows.append({
            **locus,
            "final_family": assignment["final_family"],
            "final_subfamily": assignment["final_subfamily"],
            "assignment_confidence": assignment["assignment_confidence"],
            "assignment_decision": assignment["assignment_decision"],
            "tree_nearest_reference_id": assignment["tree_nearest_reference_id"],
        })
    expanded_fields = list(expanded_rows[0])
    write_tsv(output / "locus_OR_subfamily_assignments.tsv", expanded_rows, expanded_fields)

    count_rows = []
    grouped = Counter(
        (row["analysis_unit"], row["final_family"] or "unresolved", row["final_subfamily"] or "unresolved", row["assignment_confidence"])
        for row in expanded_rows
    )
    for key, count in sorted(grouped.items()):
        count_rows.append({
            "analysis_unit": key[0], "or_family": key[1], "or_subfamily": key[2],
            "assignment_confidence": key[3], "functional_loci": count,
        })
    write_tsv(
        output / "species_OR_subfamily_counts.tsv", count_rows,
        ["analysis_unit", "or_family", "or_subfamily", "assignment_confidence", "functional_loci"],
    )

    confidence_counts = Counter(row["assignment_confidence"] for row in final_rows)
    decision_counts = Counter(row["assignment_decision"] for row in final_rows)
    checks = {
        "all_17_family_trees_parsed": True,
        "family_tree_placements_are_5957": len(placement_rows) == 5957,
        "representative_rows_are_5338": len(final_rows) == 5338,
        "representative_ids_are_unique": len(final_by_rep) == 5338,
        "locus_rows_are_6419": len(expanded_rows) == 6419,
        "locus_ids_are_unique": len({x["global_id"] for x in expanded_rows}) == 6419,
        "all_seven_analysis_units_present": {x["analysis_unit"] for x in expanded_rows}
        == {x["analysis_unit"] for x in species_rows},
        "exactly_ten_candidate_stage_unresolved": sum(
            row["candidate_status"] in {"not_assignable_with_human_reference", "no_qualified_reference_hit"}
            for row in candidate_rows
        ) == 10,
    }
    summary = {
        "schema_version": "1.0",
        "status": "pass" if all(checks.values()) else "fail",
        "method": "human-reference anchored within-family patristic placement with unrooted pure-split support audit",
        "thresholds": {
            "high_sh_alrt": args.high_sh_alrt,
            "high_ufboot": args.high_ufboot,
            "high_distance_ratio": args.high_distance_ratio,
            "medium_distance_ratio": args.medium_distance_ratio,
        },
        "counts": {
            "family_trees": 17,
            "human_reference_proteins": len(human_rows),
            "candidate_tree_placements": len(placement_rows),
            "unique_felid_representatives": len(final_rows),
            "functional_felid_loci": len(expanded_rows),
            "analysis_units": len(species_rows),
            "by_confidence": dict(sorted(confidence_counts.items())),
            "by_decision": dict(sorted(decision_counts.items())),
        },
        "checks": checks,
        "provenance": {
            "candidate_assignments_sha256": sha256(args.candidate_assignments),
            "human_metadata_sha256": sha256(args.human_metadata),
            "identical_sequence_map_sha256": sha256(args.identical_sequence_map),
            "species_metadata_sha256": sha256(args.species_metadata),
        },
        "outputs": {
            "all_candidate_tree_placements": str((output / "all_candidate_tree_placements.tsv").resolve()),
            "representative_assignments": str((output / "representative_OR_subfamily_assignments.tsv").resolve()),
            "locus_assignments": str((output / "locus_OR_subfamily_assignments.tsv").resolve()),
            "species_counts": str((output / "species_OR_subfamily_counts.tsv").resolve()),
        },
    }
    marker = output / "ASSIGNMENT_COMPLETE.json"
    marker.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(marker)
    if summary["status"] != "pass":
        raise ValueError("assignment completion checks failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

