#!/usr/bin/env python3
"""Independent file-level validation for curated FCA/PPA ChrY coding models."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from curate_felid_chry_annotation import (  # noqa: E402
    assign_region, finish_model, parse_attrs, parse_gff, read_fasta, read_regions,
)


def inspect_gff(path: Path, sequence_length: int) -> tuple[int, list[str]]:
    ids: set[str] = set()
    parents: list[tuple[str, str]] = []
    loci = 0
    errors: list[str] = []
    with path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                errors.append(f'{path.name}:{lineno}:not_9_columns')
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                errors.append(f'{path.name}:{lineno}:non_integer_coordinate')
                continue
            if start < 1 or end > sequence_length or start > end:
                errors.append(f'{path.name}:{lineno}:coordinate_out_of_bounds')
            attrs = parse_attrs(f[8])
            feature_id = attrs.get('ID')
            if feature_id:
                if feature_id in ids:
                    errors.append(f'{path.name}:{lineno}:duplicate_ID:{feature_id}')
                ids.add(feature_id)
            if attrs.get('Parent'):
                for parent in attrs['Parent'].split(','):
                    parents.append((parent, f'{path.name}:{lineno}'))
            if f[2] in {'gene','pseudogene'}:
                loci += 1
    for parent, where in parents:
        if parent not in ids:
            errors.append(f'{where}:missing_Parent:{parent}')
    return loci, errors


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--species', choices=('FCA','PPA'), required=True)
    p.add_argument('--genome', type=Path, required=True)
    p.add_argument('--regions', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    args = p.parse_args()
    species = args.species
    genome_records = read_fasta(args.genome)
    if set(genome_records) != {'ChrY'}:
        raise SystemExit('Genome must contain exactly one sequence named ChrY.')
    genome = genome_records['ChrY']
    regions = read_regions(args.regions)
    errors: list[str] = []
    rows: list[dict[str, object]] = []
    for cls in ('high_confidence','provisional'):
        gff = args.results / f'{species}.ChrY.protein_coding.{cls}.gff3'
        loci, file_errors = inspect_gff(gff, len(genome))
        errors.extend(file_errors)
        models = parse_gff(gff, f'QC_{cls}')
        cds_fa = read_fasta(args.results / f'{species}.ChrY.protein_coding.{cls}.cds.fa')
        pep_fa = read_fasta(args.results / f'{species}.ChrY.protein_coding.{cls}.pep.fa')
        if loci != len(models) or loci != len(cds_fa) or loci != len(pep_fa):
            errors.append(f'{cls}:count_mismatch:gff_loci={loci}:models={len(models)}:cds={len(cds_fa)}:pep={len(pep_fa)}')
        for model in models:
            finish_model(model, genome)
            assign_region(model, regions)
            gid = model.gene_attrs.get('ID', '')
            if cls == 'high_confidence' and model.orf_status != 'complete_intact':
                errors.append(f'{gid}:high_confidence_ORF={model.orf_status}')
            if gid not in cds_fa or cds_fa[gid] != model.cds_seq:
                errors.append(f'{gid}:CDS_FASTA_mismatch')
            if gid not in pep_fa or pep_fa[gid] != model.peptide.rstrip('*'):
                errors.append(f'{gid}:protein_FASTA_mismatch')
            declared = model.gene_attrs.get('sequence_class')
            if declared != model.region_class:
                errors.append(f'{gid}:sequence_class_mismatch:{declared}:{model.region_class}')
        rows.append({'file': gff.name, 'loci': loci, 'status': 'PASS' if not file_errors else 'FAIL'})
    pg_gff = args.results / f'{species}.ChrY.pseudogene.gff3'
    pg_loci, pg_errors = inspect_gff(pg_gff, len(genome))
    errors.extend(pg_errors)
    pg_fa = read_fasta(args.results / f'{species}.ChrY.pseudogene.cds_like.fa')
    if pg_loci != len(pg_fa):
        errors.append(f'pseudogene:count_mismatch:gff={pg_loci}:fasta={len(pg_fa)}')
    rows.append({'file': pg_gff.name, 'loci': pg_loci, 'status': 'PASS' if not pg_errors else 'FAIL'})
    with args.report.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['file','loci','status'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
        out.write(f'overall\t{sum(int(r["loci"]) for r in rows)}\t{"FAIL" if errors else "PASS"}\n')
        for error in errors:
            out.write(f'error\tNA\t{error}\n')
    if errors:
        for error in errors[:20]:
            print(error, file=sys.stderr)
        return 1
    print(f'File-level validation PASS: {sum(int(r["loci"]) for r in rows)} loci')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
