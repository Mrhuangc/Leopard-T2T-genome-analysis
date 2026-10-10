#!/usr/bin/env python3
"""Compare old and new coding loci within the same oriented ChrY assembly."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

FIELDS = ('status','family','old_id','old_class','old_start','old_end',
          'new_id','new_class','new_start','new_end','reciprocal_overlap')


@dataclass(frozen=True)
class Gene:
    gene_id: str
    family: str
    start: int
    end: int
    strand: str
    gene_class: str


def read_genes(path: Path, gene_class: str) -> list[Gene]:
    genes = []
    with path.open(encoding='utf-8') as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith('#'):
                continue
            fields = raw.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'gene':
                continue
            attrs = dict(part.split('=', 1) for part in fields[8].split(';') if '=' in part)
            family = unquote(attrs.get('gene_symbol') or attrs.get('Name') or '')
            if not family:
                raise ValueError(f'{path}: gene without symbol at {fields[3]}')
            genes.append(Gene(unquote(attrs['ID']), family, int(fields[3]), int(fields[4]),
                              fields[6], gene_class))
    return genes


def reciprocal_overlap(a: Gene, b: Gene) -> float:
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    return min(overlap/(a.end-a.start+1), overlap/(b.end-b.start+1))


def compare(old: list[Gene], new: list[Gene], min_reciprocal: float = .5) -> list[dict[str, str]]:
    candidates = []
    for i, a in enumerate(old):
        for j, b in enumerate(new):
            if a.family == b.family and a.strand == b.strand:
                score = reciprocal_overlap(a, b)
                if score >= min_reciprocal:
                    candidates.append((score, i, j))
    matched_old = set()
    matched_new = set()
    pairs = []
    for score, i, j in sorted(candidates, reverse=True):
        if i not in matched_old and j not in matched_new:
            pairs.append((i, j, score))
            matched_old.add(i)
            matched_new.add(j)
    output = []
    def row(a: Gene | None, b: Gene | None, status: str, score: float | None) -> dict[str, str]:
        return {
            'status': status, 'family': a.family if a else b.family,
            'old_id': a.gene_id if a else '', 'old_class': a.gene_class if a else '',
            'old_start': str(a.start) if a else '', 'old_end': str(a.end) if a else '',
            'new_id': b.gene_id if b else '', 'new_class': b.gene_class if b else '',
            'new_start': str(b.start) if b else '', 'new_end': str(b.end) if b else '',
            'reciprocal_overlap': f'{score:.4f}' if score is not None else '',
        }
    for i, j, score in pairs:
        a, b = old[i], new[j]
        output.append(row(a, b, 'retained' if a.gene_class == b.gene_class else 'reclassified', score))
    output.extend(row(a, None, 'baseline_only', None) for i, a in enumerate(old) if i not in matched_old)
    output.extend(row(None, b, 'new_only', None) for j, b in enumerate(new) if j not in matched_new)
    return sorted(output, key=lambda r: (r['family'], int(r['old_start'] or r['new_start']), r['status']))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ('old-high','old-provisional','new-high','new-provisional'):
        p.add_argument('--'+arg, required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--min-reciprocal', type=float, default=.5)
    args = p.parse_args()
    old = read_genes(args.old_high, 'high_confidence') + read_genes(args.old_provisional, 'provisional')
    new = read_genes(args.new_high, 'high_confidence') + read_genes(args.new_provisional, 'provisional')
    rows = compare(old, new, args.min_reciprocal)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    from collections import Counter
    print(dict(Counter(row['status'] for row in rows)))


if __name__ == '__main__':
    main()
