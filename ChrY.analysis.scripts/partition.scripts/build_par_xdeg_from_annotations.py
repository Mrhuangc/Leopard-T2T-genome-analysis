#!/usr/bin/env python3
"""Build the reviewed four-track base_input.tsv from ChrY GFF3 and partition BED.

The GFF3/BED selection follows the source project's gene-order method:
assign each gene to the region with the greatest overlap, keep PAR/X-degenerate
genes, exclude provisional PAR ORFs, then join ChrX gene annotations. Rebuild
the historical 208-row table with FCA's original and PPA's boundary-revised
partition; a newer FCA partition changes the membership by design.
"""

from __future__ import annotations

import argparse
import csv
import tempfile
from pathlib import Path
from urllib.parse import unquote

from build_par_xdeg_base_input import build


Y_ORDER_COLUMNS = (
    "species", "region_class", "region_id", "order_in_class", "family",
    "gene_id", "y_start0", "y_end0", "y_strand", "annotation_class",
)


def attributes(raw: str) -> dict[str, str]:
    return {unquote(k): unquote(v) for part in raw.strip().strip(";").split(";")
            if "=" in part for k, v in [part.split("=", 1)]}


def read_genes(path: Path, default_class: str) -> list[dict[str, object]]:
    genes = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_no}: expected 9 GFF3 columns")
            if fields[2] != "gene":
                continue
            attr = attributes(fields[8])
            gene_id = attr.get("ID")
            if not gene_id:
                raise ValueError(f"{path}:{line_no}: gene lacks ID")
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid gene coordinates")
            family = (attr.get("gene_symbol") or attr.get("Name") or
                      attr.get("gene") or gene_id).strip()
            if family.upper() == "CCDC71LY":
                family = "CCDC71L"
            genes.append({
                "chrom": fields[0], "start": start, "end": end,
                "strand": fields[6], "gene_id": gene_id, "family": family,
                "annotation_class": attr.get("annotation_class", default_class),
            })
    return genes


def read_regions(path: Path) -> list[dict[str, object]]:
    """Accept old chrom/start/end/class/id or revised chrom/start/end/id/class BED."""
    regions = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0].lower() in {"chrom", "#chrom"}:
                continue
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_no}: expected at least 5 BED columns")
            if fields[3] in {"PAR", "X-degenerate", "Ampliconic", "Repeat-rich"} or fields[3].startswith("Repeat-rich/"):
                region_class, region_id = fields[3:5]
            else:
                region_id, region_class = fields[3:5]
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid BED interval")
            regions.append({"chrom": fields[0], "start": start, "end": end,
                            "region_class": region_class, "region_id": region_id})
    if not regions:
        raise ValueError(f"No regions in {path}")
    return regions


def y_order(species: str, high: Path, provisional: Path, partition: Path) -> list[dict[str, object]]:
    regions = read_regions(partition)
    genes = read_genes(high, "high_confidence") + read_genes(provisional, "provisional")
    seen = set()
    selected = []
    for gene in genes:
        key = gene["gene_id"]
        if key in seen:
            raise ValueError(f"Duplicate {species} gene ID across GFF3s: {key}")
        seen.add(key)
        overlaps = [max(0, min(int(gene["end"]), int(r["end"])) -
                        max(int(gene["start"]), int(r["start"])))
                    if gene["chrom"] == r["chrom"] else 0 for r in regions]
        best = max(range(len(regions)), key=lambda i: overlaps[i])
        region = regions[best]
        if overlaps[best] == 0 or region["region_class"] not in {"PAR", "X-degenerate"}:
            continue
        if region["region_class"] == "PAR" and gene["annotation_class"] == "provisional":
            continue
        selected.append((gene, region))
    selected.sort(key=lambda pair: (int(pair[0]["start"]), int(pair[0]["end"]), str(pair[0]["gene_id"])))
    counters = {"PAR": 0, "X-degenerate": 0}
    result = []
    for gene, region in selected:
        category = str(region["region_class"])
        counters[category] += 1
        result.append({
            "species": species, "region_class": category, "region_id": region["region_id"],
            "order_in_class": counters[category], "family": gene["family"],
            "gene_id": gene["gene_id"], "y_start0": gene["start"],
            "y_end0": gene["end"], "y_strand": gene["strand"],
            "annotation_class": gene["annotation_class"],
        })
    return result


def write_order(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=Y_ORDER_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    for species in ("fca", "ppa"):
        parser.add_argument(f"--{species}-y-high", required=True, type=Path)
        parser.add_argument(f"--{species}-y-provisional", required=True, type=Path)
        parser.add_argument(f"--{species}-partition", required=True, type=Path)
        parser.add_argument(f"--{species}-x-gff", type=Path)
    parser.add_argument("--x-catalog", type=Path, help="Frozen ChrX subset if full ChrX GFF3s are unavailable")
    parser.add_argument("--rab9-review", type=Path,
                        default=project / "results/02_fca_pseudogenes/RAB9AY.review.tsv")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--y-order-output", type=Path, help="Write selected Y gene and partition audit")
    args = parser.parse_args()
    if bool(args.fca_x_gff) != bool(args.ppa_x_gff):
        parser.error("Both ChrX GFF3s must be supplied together")
    if not args.fca_x_gff and not args.x_catalog:
        parser.error("Supply both ChrX GFF3s, or explicitly supply --x-catalog")
    rows = (y_order("FCA", args.fca_y_high, args.fca_y_provisional, args.fca_partition)
            + y_order("PPA", args.ppa_y_high, args.ppa_y_provisional, args.ppa_partition))
    if args.y_order_output:
        write_order(rows, args.y_order_output)
    with tempfile.TemporaryDirectory() as directory:
        order_file = Path(directory) / "PAR_XDG_gene_order.tsv"
        write_order(rows, order_file)
        result = build(order_file, args.fca_x_gff, args.ppa_x_gff,
                       args.x_catalog or project / "par_xdeg_4tracks/ChrX_gene_catalog.tsv",
                       args.rab9_review, args.output)
    print(f"Wrote {len(result)} four-track rows to {args.output}")


if __name__ == "__main__":
    main()
