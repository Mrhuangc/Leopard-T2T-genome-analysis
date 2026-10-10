#!/usr/bin/env python3
"""Computational curation of felid ChrY protein-coding models from supplied evidence.

The workflow keeps gene evolutionary origin separate from the structural
sequence class.  It integrates the user's EVM annotation, named ab-initio/homology
callsets from a manifest, and independent felid ChrY protein evidence.
It emits conservative high-confidence, provisional, pseudogene, copy-status,
and audit files.  No transcript support or read-depth copy-number estimate is
asserted because those inputs are not part of this workflow.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import urllib.parse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
TARGET_SPECIES = 'FCA'
FELID_SUPPORT_SPECIES = ('PPA', 'PLE', 'PON', 'NEE')
TARGET_HOMOLOGY_SOURCE = 'GeMoMa_gp'

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

PAR_FAMILIES = {
    'AKAP17A','ANOS1','ARSD','ARSF','ARSH','ARSL','ASMT','ASMTL','CD99','CRLF2',
    'CSF2RA','DHRSX','GPR143','GTPBP6','GYG2','HADHB','IL3RA','MXRA5','NLGN4X',
    'OBP','P2RY8','PLCXD1','PNPLA4','PPP2R3B','PRKX','PUDP','SHOX','SLC25A6',
    'STAG2','STS','TBL1X','XG','ZBED1',
}
Y_GAMETOLOG_FAMILIES = {
    'AMELY','AP1S2Y','CLDN34Y','CUL4BY','DDX3Y','EEF1A1Y','EIF1AY','EIF2S3Y',
    'GPR160Y','HSFY','KDM5D','MBTPS2Y','OFD1Y','RAB9AY','RBMY','SHROOM2Y',
    'TMSB4Y','TXLNGY','UBE1Y','USP9Y','UTY','ZFY',
}
KNOWN_FAMILIES = PAR_FAMILIES | Y_GAMETOLOG_FAMILIES | {
    'CCDC71L','FELIDY_ORF1','SRY','TSPY',
}

ALIASES = {
    'AMELX':'AMELY', 'AP1S2':'AP1S2Y', 'CLDN34':'CLDN34Y',
    'CUL4A':'CUL4BY', 'CUL4B':'CUL4BY', 'DDX3X':'DDX3Y',
    'EEF1A1':'EEF1A1Y', 'EIF1AX':'EIF1AY', 'EIF2S3':'EIF2S3Y',
    'GPR160':'GPR160Y', 'HSFY1':'HSFY', 'KDM5C':'KDM5D',
    'KDM6A':'UTY', 'MBTPS2':'MBTPS2Y', 'OFD1':'OFD1Y',
    'RAB9A':'RAB9AY', 'RBMX':'RBMY', 'SHROOM2':'SHROOM2Y',
    'TMSB4':'TMSB4Y', 'TMSB4X':'TMSB4Y', 'TXLNG':'TXLNGY',
    'UBA1':'UBE1Y', 'UBA1Y':'UBE1Y', 'USP9X':'USP9Y', 'ZFX':'ZFY',
}


@dataclass
class Model:
    uid: str
    source: str
    seqid: str
    start: int
    end: int
    strand: str
    raw_tx_id: str
    raw_gene_id: str
    attrs: Dict[str, str]
    gene_attrs: Dict[str, str]
    cds: List[Tuple[int, int, int]] = field(default_factory=list)
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds_seq: str = ''
    peptide: str = ''
    orf_status: str = 'not_evaluated'
    hint: Optional[str] = None
    family: Optional[str] = None
    region_class: str = 'Unassigned'
    region_id: str = 'NA'
    preliminary_reject: str = ''
    matches: List[dict] = field(default_factory=list)
    felid_full_species: Set[str] = field(default_factory=set)
    felid_partial_species: Set[str] = field(default_factory=set)
    best_identity: float = 0.0
    best_qcov: float = 0.0
    best_tcov: float = 0.0
    best_reference: str = 'NA'
    best_reference_species: str = 'NA'

    @property
    def cds_len(self) -> int:
        return sum(e - s + 1 for s, e, _ in self.cds)

    @property
    def span(self) -> int:
        return self.end - self.start + 1


@dataclass
class ReferenceProtein:
    rid: str
    species: str
    raw_name: str
    symbol: Optional[str]
    sequence: str
    description: str


@dataclass
class Region:
    seqid: str
    start0: int
    end0: int
    sequence_class: str
    region_id: str


def parse_attrs(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in text.strip().strip(';').split(';'):
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            out[k] = urllib.parse.unquote(v.strip('"'))
        elif ' ' in item:
            k, v = item.split(' ', 1)
            out[k] = v.strip().strip('"')
    return out


def gff_escape(value: object) -> str:
    return urllib.parse.quote(str(value), safe='._,:|+-')


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    name: Optional[str] = None
    parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records[name] = ''.join(parts).upper()
    return records


def read_fasta_with_headers(path: Path) -> List[Tuple[str, str, str]]:
    records: List[Tuple[str, str, str]] = []
    name: Optional[str] = None
    header = ''
    parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, header, ''.join(parts).upper()))
                header = line[1:]
                name = header.split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records.append((name, header, ''.join(parts).upper()))
    return records


def reverse_complement(seq: str) -> str:
    return seq.translate(RC)[::-1]


def translate(seq: str) -> str:
    return ''.join(CODON_TABLE.get(seq[i:i+3], 'X')
                   for i in range(0, len(seq) - 2, 3))


def canonical_symbol(raw: Optional[str]) -> Optional[str]:
    if not raw:
        return None
    x = urllib.parse.unquote(raw).strip().split(',')[0]
    x = re.sub(r'^(?:gp|cat|human|bbear|dog|pig|sheep|cattle)_', '', x,
               flags=re.IGNORECASE)
    x = re.sub(r'_R\d+$', '', x)
    x = re.sub(r'\.t\d+$', '', x)
    x = re.sub(r'\.c\d+$', '', x)
    x = re.sub(r'-(?:\d+)$', '', x)
    x = x.upper()
    if re.match(r'^(?:LOC\d+|XM_|XP_|NM_|NP_|G\d+$|CHRY-G\d+$|FCA\.YG\.|PPA\.YG\.|AMEL\.[A-Z]+\.)', x):
        return None
    if x.startswith('TSPY'):
        return 'TSPY'
    if x.startswith('HSFY'):
        return 'HSFY'
    if x.startswith('SRY'):
        return 'SRY'
    if x.startswith('RBMY'):
        return 'RBMY'
    if x.startswith('CCDC71L'):
        return 'CCDC71L'
    if x in ALIASES:
        return ALIASES[x]
    return x


def model_hint(source: str, attrs: Dict[str, str], gene_attrs: Dict[str, str]) -> Optional[str]:
    candidates: List[Optional[str]] = []
    if source == 'EVM':
        candidates.extend([gene_attrs.get('Name'), gene_attrs.get('gene_symbol'),
                           gene_attrs.get('ID'), attrs.get('Name')])
    elif source.startswith('GeMoMa_'):
        candidates.extend([attrs.get('ref-gene'), attrs.get('Name')])
    else:
        candidates.extend([gene_attrs.get('Name'), attrs.get('Name')])
    for raw in candidates:
        value = canonical_symbol(raw)
        if value:
            return value
    return None


def parse_gff(path: Path, source: str) -> List[Model]:
    models: List[Model] = []
    gene_by_raw: Dict[str, Tuple[str, Dict[str, str]]] = {}
    active_tx: Dict[str, Model] = {}
    occurrence = collections.Counter()
    with path.open(errors='replace') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            seqid, _, ftype, start_s, end_s, _, strand, phase_s, attr_s = fields[:9]
            if not start_s.isdigit() or not end_s.isdigit():
                continue
            start, end = int(start_s), int(end_s)
            attrs = parse_attrs(attr_s)
            if ftype in {'gene', 'pseudogene'}:
                raw = attrs.get('ID', f'{source}.gene.{start}.{end}')
                occurrence[(source, 'gene', raw)] += 1
                unique = f'{source}:{raw}:{occurrence[(source,"gene",raw)]}:{start}-{end}'
                gene_by_raw[raw] = (unique, attrs)
            elif ftype in {'mRNA', 'transcript', 'pseudogenic_transcript'}:
                raw_tx = attrs.get('ID', f'{source}.tx.{start}.{end}')
                raw_gene = attrs.get('Parent', f'{source}.gene.{start}.{end}').split(',')[0]
                occurrence[(source, 'tx', raw_tx)] += 1
                uid = f'{source}:{raw_tx}:{occurrence[(source,"tx",raw_tx)]}:{start}-{end}:{strand}'
                gene_unique, gene_attrs = gene_by_raw.get(raw_gene, (raw_gene, {}))
                model = Model(uid=uid, source=source, seqid=seqid, start=start,
                              end=end, strand=strand, raw_tx_id=raw_tx,
                              raw_gene_id=gene_unique, attrs=attrs,
                              gene_attrs=dict(gene_attrs))
                model.hint = model_hint(source, attrs, gene_attrs)
                models.append(model)
                active_tx[raw_tx] = model
            elif ftype in {'CDS', 'exon', 'pseudogenic_exon'}:
                parents = attrs.get('Parent', '').split(',')
                for parent in parents:
                    model = active_tx.get(parent)
                    if model is None:
                        continue
                    if ftype == 'CDS':
                        phase = int(phase_s) if phase_s in {'0','1','2'} else 0
                        model.cds.append((start, end, phase))
                    else:
                        model.exons.append((start, end))
    return [m for m in models if m.cds]


def finish_model(model: Model, genome: str) -> None:
    ordered = sorted(model.cds, key=lambda x: x[0], reverse=model.strand == '-')
    chunks: List[str] = []
    for start, end, _ in ordered:
        seq = genome[start - 1:end]
        if model.strand == '-':
            seq = reverse_complement(seq)
        chunks.append(seq)
    seq = ''.join(chunks)
    first_phase = ordered[0][2] if ordered else 0
    if first_phase:
        seq = seq[first_phase:]
    pep = translate(seq)
    model.cds_seq = seq
    model.peptide = pep
    if not seq:
        model.orf_status = 'empty_cds'
    elif len(seq) % 3 != 0:
        model.orf_status = 'frameshift_length'
    elif 'X' in pep:
        model.orf_status = 'ambiguous_sequence'
    elif '*' in pep[:-1]:
        model.orf_status = 'internal_stop'
    elif not pep.endswith('*'):
        model.orf_status = 'missing_stop'
    elif not pep.startswith('M'):
        model.orf_status = 'complete_noncanonical_start'
    else:
        model.orf_status = 'complete_intact'


def read_regions(path: Path) -> List[Region]:
    regions: List[Region] = []
    header: List[str] | None = None
    with path.open() as fh:
        for line_number, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.lstrip('#').rstrip('\n').split('\t')
            if f[0].lower() == 'chrom':
                header = [x.lower() for x in f]
                continue
            if line.startswith('#'):
                continue
            if header is not None:
                if 'class' not in header or 'region_id' not in header:
                    raise ValueError(f'{path}: header needs class and region_id')
                region_class = f[header.index('class')]
                region_id = f[header.index('region_id')]
            else:
                if len(f) < 5:
                    raise ValueError(f'{path}:{line_number}: expected 5 columns')
                # Legacy FCA BED is chrom/start/end/class/region_id; most
                # revised beds are chrom/start/end/region_id/class.
                region_class, region_id = (f[3], f[4]) if f[3] in {
                    'PAR', 'Ampliconic', 'X-degenerate', 'Repeat-rich',
                    'Repeat-rich/ampliconic complex', 'Other MSY'
                } else (f[4], f[3])
            regions.append(Region(f[0], int(f[1]), int(f[2]), region_class, region_id))
    ordered = sorted(regions, key=lambda r: r.start0)
    if not ordered or ordered[0].start0 != 0:
        raise ValueError(f'{path}: partition must start at 0')
    for left, right in zip(ordered, ordered[1:]):
        if left.seqid != right.seqid or left.end0 != right.start0:
            raise ValueError(f'{path}: partition gap, overlap or mixed chromosome')
    return regions


def interval_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def assign_region(model: Model, regions: Sequence[Region]) -> None:
    best: Optional[Region] = None
    best_bp = -1
    intervals = [(s, e) for s, e, _ in model.cds]
    for region in regions:
        if region.seqid != model.seqid:
            continue
        bp = sum(interval_overlap(iv, (region.start0 + 1, region.end0)) for iv in intervals)
        if bp > best_bp:
            best_bp = bp
            best = region
    if best:
        model.region_class = best.sequence_class
        model.region_id = best.region_id


def model_overlap(a: Model, b: Model) -> float:
    if a.seqid != b.seqid or a.strand != b.strand:
        return 0.0
    total = 0
    for sa, ea, _ in a.cds:
        for sb, eb, _ in b.cds:
            total += interval_overlap((sa, ea), (sb, eb))
    denom = min(a.cds_len, b.cds_len)
    return total / denom if denom else 0.0


def parse_float(attrs: Dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(attrs.get(key, default))
    except (TypeError, ValueError):
        return default


def add_reference_fasta(refs: List[ReferenceProtein], path: Path, species: str) -> None:
    for idx, (name, header, seq) in enumerate(read_fasta_with_headers(path), 1):
        cleaned = seq.replace('*', '')
        if len(cleaned) < 20:
            continue
        refs.append(ReferenceProtein(
            rid=f'{species}.{idx:05d}', species=species, raw_name=name,
            symbol=canonical_symbol(name), sequence=cleaned, description=header))


def build_gff_references(refs: List[ReferenceProtein], gff: Path,
                         fasta: Path, species: str) -> None:
    genome_records = read_fasta(fasta)
    if len(genome_records) != 1:
        raise ValueError(f'Expected one {species} ChrY sequence in {fasta}')
    genome = next(iter(genome_records.values()))
    models = parse_gff(gff, f'{species}_reference')
    for idx, model in enumerate(models, 1):
        finish_model(model, genome)
        symbol = canonical_symbol(model.gene_attrs.get('gene_symbol') or
                                  model.gene_attrs.get('Name') or
                                  model.hint)
        seq = model.peptide.replace('*', '')
        if len(seq) < 20:
            continue
        refs.append(ReferenceProtein(
            rid=f'{species}.{idx:05d}', species=species,
            raw_name=model.gene_attrs.get('ID', model.raw_gene_id),
            symbol=symbol, sequence=seq,
            description=f'{model.seqid}:{model.start}-{model.end}({model.strand})'))


def protein_kmers(seq: str, k: int = 5) -> Set[str]:
    return {seq[i:i+k] for i in range(max(0, len(seq) - k + 1))
            if 'X' not in seq[i:i+k] and '*' not in seq[i:i+k]}


def build_alignment_pairs(models: Sequence[Model], refs: Sequence[ReferenceProtein],
                          pairs_path: Path) -> Dict[str, Tuple[int, int]]:
    kmer_to_refs: Dict[str, List[int]] = collections.defaultdict(list)
    ref_kmers: List[Set[str]] = []
    for ridx, ref in enumerate(refs):
        kmers = protein_kmers(ref.sequence)
        ref_kmers.append(kmers)
        for kmer in kmers:
            kmer_to_refs[kmer].append(ridx)
    common = {k for k, ids in kmer_to_refs.items() if len(ids) > 120}
    family_refs: Dict[str, List[int]] = collections.defaultdict(list)
    for ridx, ref in enumerate(refs):
        if ref.symbol:
            family_refs[ref.symbol].append(ridx)
    pair_map: Dict[str, Tuple[int, int]] = {}
    with pairs_path.open('w') as out:
        for midx, model in enumerate(models):
            seq = model.peptide.replace('*', '')
            if len(seq) < 15:
                continue
            counts: Dict[int, int] = collections.Counter()
            qkmers = protein_kmers(seq)
            for kmer in qkmers:
                if kmer in common:
                    continue
                for ridx in kmer_to_refs.get(kmer, []):
                    counts[ridx] += 1
            by_species: Dict[str, List[Tuple[float, int]]] = collections.defaultdict(list)
            for ridx, count in counts.items():
                denom = max(1, min(len(qkmers), len(ref_kmers[ridx])))
                by_species[refs[ridx].species].append((count / denom, ridx))
            selected: Set[int] = set()
            for species, values in by_species.items():
                values.sort(reverse=True)
                selected.update(ridx for _, ridx in values[:3])
            if model.hint:
                same = family_refs.get(model.hint, [])
                same = sorted(same, key=lambda i: abs(len(refs[i].sequence) - len(seq)))
                selected.update(same[:3])
            for ridx in sorted(selected):
                ref = refs[ridx]
                if len(seq) * len(ref.sequence) > 29000000:
                    continue
                pair_id = f'M{midx:05d}.R{ridx:05d}'
                pair_map[pair_id] = (midx, ridx)
                out.write(f'{pair_id}\t{seq}\t{ref.sequence}\n')
    return pair_map


def load_alignment_results(path: Path, pair_map: Dict[str, Tuple[int, int]],
                           models: Sequence[Model], refs: Sequence[ReferenceProtein]) -> None:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row['status'] != 'ok' or row['pair_id'] not in pair_map:
                continue
            midx, ridx = pair_map[row['pair_id']]
            ref = refs[ridx]
            match = {
                'ref_index': ridx, 'rid': ref.rid, 'species': ref.species,
                'raw_name': ref.raw_name, 'symbol': ref.symbol,
                'identity': float(row['identity']),
                'qcov': float(row['query_coverage']),
                'tcov': float(row['target_coverage']),
            }
            match['similarity'] = match['identity'] * min(match['qcov'], match['tcov'])
            models[midx].matches.append(match)


def resolve_family_and_support(model: Model) -> None:
    matches = sorted(model.matches, key=lambda x: x['similarity'], reverse=True)
    named = [m for m in matches if m['symbol'] and m['similarity'] >= 0.30]
    best_named = named[0] if named else None
    if model.hint in KNOWN_FAMILIES:
        same = [m for m in named if m['symbol'] == model.hint]
        best_same = same[0] if same else None
        # GeMoMa can force a paralogous reference through a neighbouring gene.
        # A substantially stronger named felid match therefore overrides the
        # input hint; otherwise the curated EVM/GeMoMa name is retained.
        if (best_named and best_named['symbol'] != model.hint and
                best_named['similarity'] >= 0.65 and
                (best_same is None or best_named['similarity'] >= best_same['similarity'] + 0.08)):
            model.family = best_named['symbol']
        elif best_same and best_same['similarity'] >= 0.35:
            model.family = model.hint
        elif best_named and best_named['similarity'] >= 0.55:
            model.family = best_named['symbol']
        else:
            model.family = model.hint
    elif best_named and best_named['similarity'] >= 0.40:
        model.family = best_named['symbol']
    elif model.hint:
        model.family = model.hint

    if matches:
        best = matches[0]
        model.best_identity = best['identity']
        model.best_qcov = best['qcov']
        model.best_tcov = best['tcov']
        model.best_reference = best['raw_name']
        model.best_reference_species = best['species']
    by_species: Dict[str, List[dict]] = collections.defaultdict(list)
    for match in matches:
        by_species[match['species']].append(match)
    for species in FELID_SUPPORT_SPECIES:
        if species not in by_species:
            continue
        compatible = [m for m in by_species[species]
                      if m['symbol'] is None or model.family is None or m['symbol'] == model.family]
        if not compatible:
            continue
        best = max(compatible, key=lambda x: x['similarity'])
        if best['identity'] >= 0.65 and best['qcov'] >= 0.80 and best['tcov'] >= 0.80:
            model.felid_full_species.add(species)
        elif best['identity'] >= 0.50 and best['qcov'] >= 0.40 and best['tcov'] >= 0.40:
            model.felid_partial_species.add(species)


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def preliminary_filter(model: Model) -> None:
    if model.cds_len < 60:
        model.preliminary_reject = 'CDS_shorter_than_60_bp'
        return
    # Family-specific architecture guardrails are derived conservatively from
    # the original PPA curated full-length loci.  They prevent repeat-rich sequence from
    # being stitched into implausibly long ampliconic gene models.
    max_span = {
        'TSPY': 10000, 'CCDC71L': 8000, 'HSFY': 6000, 'SRY': 3000,
        'GPR160Y': 5000, 'CLDN34Y': 3000, 'FELIDY_ORF1': 3000,
        'TMSB4Y': 3000, 'RAB9AY': 5000, 'EEF1A1Y': 5000,
    }.get(model.family)
    if max_span is not None and model.span > max_span:
        model.preliminary_reject = f'{model.family}_model_span_exceeds_{max_span}_bp'
        return
    ratio = model.span / max(1, model.cds_len)
    if model.source.startswith('GeMoMa_'):
        paa = parse_float(model.attrs, 'pAA')
        iaa = parse_float(model.attrs, 'iAA')
        if ratio > 1000 and paa < 0.80:
            model.preliminary_reject = 'pathological_homology_span'
        elif paa and paa < 0.45 and model.best_identity < 0.65:
            model.preliminary_reject = 'weak_GeMoMa_projection'
        elif iaa and iaa < 0.25 and model.best_identity < 0.65:
            model.preliminary_reject = 'weak_GeMoMa_identity'
    elif model.source != 'EVM' and ratio > 2000:
        model.preliminary_reject = 'pathological_ab_initio_span'


def cluster_models(models: Sequence[Model]) -> List[List[Model]]:
    eligible = [m for m in models if not m.preliminary_reject]
    uf = UnionFind(len(eligible))
    groups: Dict[Tuple[str, str, str], List[int]] = collections.defaultdict(list)
    for i, model in enumerate(eligible):
        family = model.family or 'UNRESOLVED'
        groups[(model.seqid, model.strand, family)].append(i)
    for _, indices in groups.items():
        indices.sort(key=lambda i: eligible[i].start)
        for pos, i in enumerate(indices):
            a = eligible[i]
            for j in indices[pos + 1:]:
                b = eligible[j]
                if b.start > a.end + 10000 and a.end < b.start:
                    break
                if model_overlap(a, b) >= 0.25:
                    uf.union(i, j)
    clusters: Dict[int, List[Model]] = collections.defaultdict(list)
    for i, model in enumerate(eligible):
        clusters[uf.find(i)].append(model)
    return list(clusters.values())


def representative_score(model: Model) -> float:
    status_score = {
        'complete_intact': 100, 'complete_noncanonical_start': 65,
        'missing_stop': 35, 'frameshift_length': 15,
        'internal_stop': 10, 'ambiguous_sequence': 5,
    }.get(model.orf_status, 0)
    source_score = 0
    if model.source == 'EVM':
        source_score = 45
    elif model.source == 'ANNEVO':
        source_score = 28
    elif model.source == 'Tiberius':
        source_score = 24
    elif model.source == 'Augustus':
        source_score = 18
    elif model.source.startswith('GeMoMa_'):
        source_score = 20 + 15 * parse_float(model.attrs, 'pAA') + 15 * parse_float(model.attrs, 'iAA')
    homology = 70 * model.best_identity * min(model.best_qcov, model.best_tcov)
    length_penalty = min(30.0, max(0.0, math.log10(max(1.0, model.span / max(1, model.cds_len))) - 2.0) * 15)
    return status_score + source_score + homology - length_penalty


def source_category(source: str) -> str:
    if source == 'EVM':
        return 'EVM'
    if source in {'ANNEVO','Tiberius','Augustus'}:
        return f'ab_initio:{source}'
    if source.startswith('GeMoMa_'):
        return f'homology:{source.split("_",1)[1]}'
    return source


def is_te_family(family: Optional[str]) -> bool:
    if not family:
        return False
    return bool(re.search(r'(ERV|RETRO|TRANSPOS|GAG|POL|LINE|LTR)', family, re.I))


def classify_cluster(cluster: List[Model]) -> Tuple[str, Model, dict]:
    representative = max(cluster, key=representative_score)
    intact = [m for m in cluster if m.orf_status == 'complete_intact']
    if intact:
        representative = max(intact, key=representative_score)
    sources = {source_category(m.source) for m in cluster}
    ab_initio = {m.source for m in cluster if m.source in {'ANNEVO','Tiberius','Augustus'}}
    homology_refs = {m.source.split('_',1)[1] for m in cluster if m.source.startswith('GeMoMa_')}
    has_evm = any(m.source == 'EVM' for m in cluster)
    family = representative.family or next((m.family for m in cluster if m.family), None)
    if family and representative.family != family:
        representative.family = family
    full_felid = set().union(*(m.felid_full_species for m in cluster))
    partial_felid = set().union(*(m.felid_partial_species for m in cluster))
    gp_y_support = any(
        m.source == TARGET_HOMOLOGY_SOURCE and m.hint in KNOWN_FAMILIES and
        parse_float(m.attrs, 'pAA') >= 0.70 and parse_float(m.attrs, 'iAA') >= 0.60
        for m in cluster)
    known = family in KNOWN_FAMILIES
    te_like = is_te_family(family)
    status = representative.orf_status
    classification = 'rejected'
    reason = 'insufficient_independent_support'
    if status == 'complete_intact':
        if te_like:
            if has_evm and len(sources) >= 3:
                classification, reason = 'provisional', 'intact_TE_like_ORF_with_integrated_support'
        elif len(full_felid) >= 2:
            classification, reason = 'high_confidence', 'full_length_support_from_at_least_two_felids'
        elif len(full_felid) >= 1 and (has_evm or len(sources) >= 2):
            classification, reason = 'high_confidence', 'full_length_felid_support_plus_integrated_or_independent_model'
        elif has_evm and known and len(sources) >= 2:
            classification, reason = 'high_confidence', 'intact_EVM_model_plus_independent_prediction'
        elif gp_y_support and ab_initio and homology_refs:
            classification, reason = 'high_confidence', 'intact_known_Y_model_with_gp_and_ab_initio_support'
        elif known and (has_evm or len(sources) >= 2 or partial_felid):
            classification, reason = 'provisional', 'intact_known_family_but_incomplete_independent_support'
        elif has_evm and len(sources) >= 2:
            classification, reason = 'provisional', 'intact_integrated_unresolved_ORF'
        elif full_felid:
            classification, reason = 'provisional', 'single_source_full_length_felid_projection'
    elif status == 'complete_noncanonical_start':
        if len(full_felid) >= 2 or (known and len(sources) >= 3):
            classification, reason = 'provisional', 'noncanonical_start_with_strong_homology_support'
    else:
        disrupted_models = sum(m.orf_status in {'frameshift_length','internal_stop'} for m in cluster)
        if known and (full_felid or partial_felid) and disrupted_models >= 2 and not intact:
            classification, reason = 'pseudogene', 'recurrent_ORF_disruption_in_homology_supported_locus'
        elif known and has_evm and (full_felid or partial_felid):
            classification, reason = 'provisional', 'likely_model_error_or_partial_gene_not_proven_pseudogene'
    meta = {
        'family': family or 'UNRESOLVED', 'sources': sorted(sources),
        'ab_initio_sources': sorted(ab_initio), 'homology_refs': sorted(homology_refs),
        'has_evm': has_evm, 'full_felid': sorted(full_felid),
        'partial_felid': sorted(partial_felid), 'gp_y_support': gp_y_support,
        'reason': reason, 'model_count': len(cluster),
    }
    return classification, representative, meta


def gene_origin(family: str, region_class: str) -> str:
    if family in PAR_FAMILIES and region_class == 'PAR':
        return 'PAR_shared'
    if family == 'SRY':
        return 'sex_determining_Y'
    if family == 'TSPY':
        return 'ancestral_Y_ampliconic'
    if family == 'CCDC71L':
        return 'autosomal_derived_ampliconic'
    if family == 'FELIDY_ORF1':
        return 'felid_Y_acquired_unresolved'
    if is_te_family(family):
        return 'transposable_element_derived'
    if family in Y_GAMETOLOG_FAMILIES or region_class in {'X-degenerate','Ampliconic'}:
        return 'X_derived_or_gametolog'
    if family in PAR_FAMILIES:
        return 'PAR_shared_or_relocated'
    return 'unresolved'


def format_fasta(seq: str, width: int = 60) -> str:
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))


def write_gff(path: Path, calls: List[dict], annotation_class: str,
              sequence_length: int) -> None:
    prefix = {'high_confidence':'HC','provisional':'PR','pseudogene':'PG'}[annotation_class]
    counters: Dict[str, int] = collections.Counter()
    with path.open('w') as out:
        out.write('##gff-version 3\n')
        out.write(f'##sequence-region ChrY 1 {sequence_length}\n')
        out.write(f'# {TARGET_SPECIES} ChrY protein-coding annotation, computationally curated with supplied evidence\n')
        out.write('# Evidence: source GFF manifest and named felid Y protein references; see run manifest.\n')
        out.write('# Coordinates are 1-based closed; structural sequence class is independent of gene origin.\n')
        if annotation_class == 'high_confidence':
            out.write('# Inclusion: intact ORF plus integrated support or full-length independent felid ChrY protein evidence.\n')
        elif annotation_class == 'provisional':
            out.write('# Inclusion: plausible coding locus lacking full-length, multi-source, or non-TE evidence; manual review required.\n')
        else:
            out.write('# Inclusion: homology-supported locus with recurrent ORF disruption and no intact alternative model; no functional CDS asserted.\n')
        for call in sorted(calls, key=lambda x: (x['model'].start, x['model'].end, x['family'])):
            model: Model = call['model']
            family = call['family']
            counters[family] += 1
            safe_family = re.sub(r'[^A-Za-z0-9_.-]+', '_', family)
            gene_id = f'{TARGET_SPECIES}.ChrY.{prefix}.{safe_family}.{counters[family]:03d}'
            call['gene_id'] = gene_id
            attrs = {
                'ID': gene_id, 'Name': family, 'gene_symbol': family,
                'gene_biotype': 'pseudogene' if annotation_class == 'pseudogene' else 'protein_coding',
                'annotation_class': annotation_class,
                'sequence_class': model.region_class, 'region_id': model.region_id,
                'gene_origin': gene_origin(family, model.region_class),
                'original_model': model.raw_gene_id, 'representative_source': model.source,
                'orf_status': model.orf_status,
                'evidence': call['reason'],
                'independent_model_sources': len(call['sources']),
                'model_sources': ','.join(call['sources']) or 'NA',
                'felid_full_length_sources': ','.join(call['full_felid']) or 'NA',
                'felid_partial_sources': ','.join(call['partial_felid']) or 'NA',
                'best_reference': model.best_reference,
                'best_reference_species': model.best_reference_species,
                'protein_identity': f'{model.best_identity:.4f}',
                'protein_query_coverage': f'{model.best_qcov:.4f}',
                'protein_target_coverage': f'{model.best_tcov:.4f}',
            }
            attr_text = ';'.join(f'{k}={gff_escape(v)}' for k, v in attrs.items())
            feature = 'pseudogene' if annotation_class == 'pseudogene' else 'gene'
            out.write(f'{model.seqid}\t{TARGET_SPECIES}_ChrY_curated\t{feature}\t{model.start}\t{model.end}\t.\t{model.strand}\t.\t{attr_text}\n')
            if annotation_class == 'pseudogene':
                tx_id = gene_id + '.pt1'
                out.write(f'{model.seqid}\t{TARGET_SPECIES}_ChrY_curated\tpseudogenic_transcript\t{model.start}\t{model.end}\t.\t{model.strand}\t.\tID={tx_id};Parent={gene_id};Name={gff_escape(family)}\n')
                for i, (start, end, _) in enumerate(sorted(model.cds), 1):
                    out.write(f'{model.seqid}\t{TARGET_SPECIES}_ChrY_curated\tpseudogenic_exon\t{start}\t{end}\t.\t{model.strand}\t.\tID={tx_id}.pexon{i};Parent={tx_id}\n')
            else:
                tx_id = gene_id + '.t1'
                out.write(f'{model.seqid}\t{TARGET_SPECIES}_ChrY_curated\tmRNA\t{model.start}\t{model.end}\t.\t{model.strand}\t.\tID={tx_id};Parent={gene_id};Name={gff_escape(family)};annotation_class={annotation_class}\n')
                exon_intervals = sorted(set(model.exons or [(s, e) for s, e, _ in model.cds]))
                for i, (start, end) in enumerate(exon_intervals, 1):
                    out.write(f'{model.seqid}\t{TARGET_SPECIES}_ChrY_curated\texon\t{start}\t{end}\t.\t{model.strand}\t.\tID={tx_id}.exon{i};Parent={tx_id}\n')
                for i, (start, end, phase) in enumerate(sorted(model.cds), 1):
                    out.write(f'{model.seqid}\t{TARGET_SPECIES}_ChrY_curated\tCDS\t{start}\t{end}\t.\t{model.strand}\t{phase}\tID={tx_id}.CDS{i};Parent={tx_id}\n')


def write_sequence_files(outdir: Path, calls_by_class: Dict[str, List[dict]]) -> None:
    for cls in ('high_confidence','provisional'):
        cds_path = outdir / f'{TARGET_SPECIES}.ChrY.protein_coding.{cls}.cds.fa'
        pep_path = outdir / f'{TARGET_SPECIES}.ChrY.protein_coding.{cls}.pep.fa'
        with cds_path.open('w') as cds_out, pep_path.open('w') as pep_out:
            for call in sorted(calls_by_class[cls], key=lambda x: x['model'].start):
                model: Model = call['model']
                gid = call['gene_id']
                desc = f'gene={call["family"]} ChrY:{model.start}-{model.end}({model.strand}) class={cls} orf={model.orf_status}'
                cds_out.write(f'>{gid} {desc}\n{format_fasta(model.cds_seq)}\n')
                pep_out.write(f'>{gid} {desc}\n{format_fasta(model.peptide.rstrip("*"))}\n')
    pg_cds = outdir / f'{TARGET_SPECIES}.ChrY.pseudogene.cds_like.fa'
    with pg_cds.open('w') as out:
        for call in sorted(calls_by_class['pseudogene'], key=lambda x: x['model'].start):
            model: Model = call['model']
            out.write(f'>{call["gene_id"]} gene={call["family"]} ChrY:{model.start}-{model.end}({model.strand}) orf={model.orf_status};not_functional_CDS\n')
            out.write(format_fasta(model.cds_seq) + '\n')


def write_audit(outdir: Path, all_calls: List[dict], rejected_models: Sequence[Model]) -> None:
    fields = [
        'annotation_class','family','gene_id','seqid','start','end','strand',
        'sequence_class','region_id','gene_origin','orf_status','representative_source',
        'original_model','cluster_model_count','model_sources','ab_initio_sources',
        'homology_references','felid_full_length_sources','felid_partial_sources',
        'best_reference_species','best_reference','protein_identity',
        'protein_query_coverage','protein_target_coverage','decision_reason',
    ]
    with (outdir / f'{TARGET_SPECIES}.ChrY.annotation_evidence.audit.tsv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for call in sorted(all_calls, key=lambda x: x['model'].start):
            m: Model = call['model']
            writer.writerow({
                'annotation_class': call['classification'], 'family': call['family'],
                'gene_id': call.get('gene_id','NA'), 'seqid': m.seqid,
                'start': m.start, 'end': m.end, 'strand': m.strand,
                'sequence_class': m.region_class, 'region_id': m.region_id,
                'gene_origin': gene_origin(call['family'], m.region_class),
                'orf_status': m.orf_status, 'representative_source': m.source,
                'original_model': m.raw_gene_id,
                'cluster_model_count': call['model_count'],
                'model_sources': ','.join(call['sources']) or 'NA',
                'ab_initio_sources': ','.join(call['ab_initio_sources']) or 'NA',
                'homology_references': ','.join(call['homology_refs']) or 'NA',
                'felid_full_length_sources': ','.join(call['full_felid']) or 'NA',
                'felid_partial_sources': ','.join(call['partial_felid']) or 'NA',
                'best_reference_species': m.best_reference_species,
                'best_reference': m.best_reference,
                'protein_identity': f'{m.best_identity:.4f}',
                'protein_query_coverage': f'{m.best_qcov:.4f}',
                'protein_target_coverage': f'{m.best_tcov:.4f}',
                'decision_reason': call['reason'],
            })
    with (outdir / f'{TARGET_SPECIES}.ChrY.rejected_models.audit.tsv').open('w', newline='') as fh:
        columns = ['uid','source','start','end','strand','hint','family','orf_status',
                   'preliminary_reject','best_reference_species','best_reference',
                   'protein_identity','protein_query_coverage','protein_target_coverage']
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for m in sorted(rejected_models, key=lambda x: (x.start, x.end, x.source)):
            writer.writerow({
                'uid':m.uid,'source':m.source,'start':m.start,'end':m.end,
                'strand':m.strand,'hint':m.hint or 'NA','family':m.family or 'NA',
                'orf_status':m.orf_status,'preliminary_reject':m.preliminary_reject or 'cluster_rejected',
                'best_reference_species':m.best_reference_species,
                'best_reference':m.best_reference,
                'protein_identity':f'{m.best_identity:.4f}',
                'protein_query_coverage':f'{m.best_qcov:.4f}',
                'protein_target_coverage':f'{m.best_tcov:.4f}',
            })


def write_copy_status(outdir: Path, calls_by_class: Dict[str, List[dict]]) -> None:
    by_family: Dict[str, Dict[str, List[dict]]] = collections.defaultdict(lambda: collections.defaultdict(list))
    for cls, calls in calls_by_class.items():
        for call in calls:
            by_family[call['family']][cls].append(call)
    fields = [
        'family','gene_origin','assembly_high_confidence_coding_cn',
        'assembly_provisional_cn','assembly_provisional_complete_intact_orf_cn',
        'assembly_pseudogene_cn','assembly_total_homologous_loci','depth_based_cn',
        'depth_based_cn_status','sequence_classes','high_confidence_loci_1based',
        'provisional_loci_1based','pseudogene_loci_1based','reference_sources',
        'copy_status','notes',
    ]
    with (outdir / f'{TARGET_SPECIES}.ChrY.copy_status.tsv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for family in sorted(by_family):
            d = by_family[family]
            all_calls = sum(d.values(), [])
            classes = sorted({c['model'].region_class for c in all_calls})
            refs = sorted({sp for c in all_calls for sp in c['full_felid'] + c['partial_felid']})
            def loci(cls: str) -> str:
                vals = []
                for c in sorted(d.get(cls, []), key=lambda x: x['model'].start):
                    m = c['model']
                    vals.append(f'{c["gene_id"]}@{m.seqid}:{m.start}-{m.end}({m.strand})')
                return ';'.join(vals) or 'NA'
            hc = len(d.get('high_confidence', []))
            pr = len(d.get('provisional', []))
            pg = len(d.get('pseudogene', []))
            writer.writerow({
                'family':family,
                'gene_origin':','.join(sorted({gene_origin(family,c['model'].region_class) for c in all_calls})),
                'assembly_high_confidence_coding_cn':hc,
                'assembly_provisional_cn':pr,
                'assembly_provisional_complete_intact_orf_cn':sum(c['model'].orf_status=='complete_intact' for c in d.get('provisional',[])),
                'assembly_pseudogene_cn':pg,
                'assembly_total_homologous_loci':hc+pr+pg,
                'depth_based_cn':'NA',
                'depth_based_cn_status':'not_estimated_no_raw_read_depth_input',
                'sequence_classes':','.join(classes),
                'high_confidence_loci_1based':loci('high_confidence'),
                'provisional_loci_1based':loci('provisional'),
                'pseudogene_loci_1based':loci('pseudogene'),
                'reference_sources':','.join(refs) or 'EVM_and_prediction_evidence_only',
                'copy_status':'multicopy' if hc+pr+pg > 1 else 'single_locus',
                'notes':'Assembly-locus count only; no read-depth CN estimate. Provisional and pseudogene loci are reported separately.',
            })


def write_original_disposition(outdir: Path, evm_models: Sequence[Model], all_calls: Sequence[dict]) -> None:
    with (outdir / f'{TARGET_SPECIES}.ChrY.original_EVM_model_disposition.tsv').open('w', newline='') as fh:
        fields = ['original_gene','original_transcript','start','end','original_orf_status',
                  'assigned_family','final_class','final_gene_id','action','reason']
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for evm in sorted(evm_models, key=lambda x: x.start):
            overlapping = []
            for call in all_calls:
                m = call['model']
                if m.seqid == evm.seqid and interval_overlap((m.start,m.end),(evm.start,evm.end)) > 0:
                    if model_overlap(m, evm) >= 0.20 or model_overlap(evm, m) >= 0.20:
                        overlapping.append(call)
            overlapping.sort(key=lambda c: (c['classification']!='high_confidence', -model_overlap(c['model'], evm)))
            best = overlapping[0] if overlapping else None
            if not best:
                final_class, gid, family, action, reason = 'rejected','NA',evm.family or evm.hint or 'UNRESOLVED','removed','no_final_supported_locus'
            else:
                final_class, gid, family = best['classification'], best.get('gene_id','NA'), best['family']
                if best['model'].source == 'EVM':
                    action = 'retained_and_standardized'
                else:
                    action = 'replaced_by_better_supported_model'
                reason = best['reason']
            writer.writerow({
                'original_gene':evm.gene_attrs.get('ID',evm.raw_gene_id),
                'original_transcript':evm.raw_tx_id,'start':evm.start,'end':evm.end,
                'original_orf_status':evm.orf_status,'assigned_family':family,
                'final_class':final_class,'final_gene_id':gid,'action':action,'reason':reason,
            })


def validate_outputs(outdir: Path, sequence_length: int, calls_by_class: Dict[str, List[dict]]) -> List[str]:
    messages: List[str] = []
    ids: Set[str] = set()
    errors: List[str] = []
    for cls, calls in calls_by_class.items():
        for call in calls:
            gid = call.get('gene_id')
            m: Model = call['model']
            if not gid or gid in ids:
                errors.append(f'duplicate_or_missing_gene_id:{gid}')
            ids.add(gid)
            if m.start < 1 or m.end > sequence_length or m.start > m.end:
                errors.append(f'coordinate_out_of_range:{gid}')
            if cls == 'high_confidence' and m.orf_status != 'complete_intact':
                errors.append(f'high_confidence_not_complete_intact:{gid}:{m.orf_status}')
            if cls != 'pseudogene' and not m.cds:
                errors.append(f'coding_gene_without_CDS:{gid}')
    messages.append(f'Validated {len(ids)} nonredundant final loci.')
    messages.append(f'Coordinate bounds: 1-{sequence_length}.')
    messages.append('High-confidence loci require complete_intact ORFs.')
    messages.append('Copy number is assembly-locus count only; raw-read depth was not supplied.')
    messages.append(f'RNA evidence status: not evaluated for {TARGET_SPECIES} ChrY; no expression claim is made.')
    if errors:
        messages.extend('ERROR: ' + e for e in errors)
        raise RuntimeError('\n'.join(errors))
    messages.append('QC status: PASS')
    return messages


def write_summary(outdir: Path, calls_by_class: Dict[str, List[dict]],
                  models: Sequence[Model], qc_messages: Sequence[str]) -> None:
    families = {c['family'] for calls in calls_by_class.values() for c in calls}
    with (outdir / f'{TARGET_SPECIES}.ChrY.annotation_summary.tsv').open('w') as out:
        out.write('metric\tvalue\n')
        out.write(f'input_candidate_transcripts\t{len(models)}\n')
        for cls in ('high_confidence','provisional','pseudogene'):
            out.write(f'{cls}_loci\t{len(calls_by_class[cls])}\n')
        out.write(f'final_gene_families\t{len(families)}\n')
        for cls, calls in calls_by_class.items():
            by_region = collections.Counter(c['model'].region_class for c in calls)
            for region, count in sorted(by_region.items()):
                out.write(f'{cls}_loci_in_{region}\t{count}\n')
    with (outdir / f'{TARGET_SPECIES}.ChrY.annotation_QC.txt').open('w') as out:
        out.write('\n'.join(qc_messages) + '\n')


def read_evidence_manifest(path: Path) -> List[Tuple[str, Path]]:
    """Read source/path pairs; paths relative to the TSV are allowed."""
    sources: List[Tuple[str, Path]] = []
    seen: Set[str] = {'EVM'}
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or not {'source','path'} <= set(reader.fieldnames):
            raise ValueError(f'{path}: expected source and path columns')
        for line_number, row in enumerate(reader, 2):
            source = (row.get('source') or '').strip()
            value = (row.get('path') or '').strip()
            if not source or not value or source in seen:
                raise ValueError(f'{path}:{line_number}: invalid/duplicate source or empty path')
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', source):
                raise ValueError(f'{path}:{line_number}: invalid source name {source!r}')
            source_path = (path.parent / value).resolve()
            if not source_path.is_file():
                raise FileNotFoundError(source_path)
            seen.add(source)
            sources.append((source, source_path))
    return sources


def parse_protein_refs(values: Sequence[str], target_species: str) -> List[Tuple[str, Path]]:
    refs: List[Tuple[str, Path]] = []
    seen: Set[str] = set()
    for item in values:
        if '=' not in item:
            raise ValueError(f'--protein-ref requires SPECIES=FASTA: {item}')
        species, file_name = item.split('=', 1)
        species = species.upper()
        path = Path(file_name).resolve()
        if not re.fullmatch(r'[A-Z]{3}', species) or species == target_species or species in seen:
            raise ValueError(f'Invalid/duplicate reference species: {species}')
        if not path.is_file():
            raise FileNotFoundError(path)
        seen.add(species)
        refs.append((species, path))
    return refs


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species', choices=('FCA', 'PPA'), required=True)
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--evm-gff', type=Path, required=True)
    parser.add_argument('--evidence-manifest', type=Path, required=True)
    parser.add_argument('--regions', type=Path, required=True)
    parser.add_argument('--reference-species', choices=('FCA', 'PPA'),
                        help='Optional existing GFF reference; supply all three --reference-* options together')
    parser.add_argument('--reference-gff', type=Path)
    parser.add_argument('--reference-genome', type=Path)
    parser.add_argument('--protein-ref', action='append', default=[], metavar='SPECIES=FASTA')
    parser.add_argument('--target-homology-source', default='GeMoMa_gp')
    parser.add_argument('--aligner', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    reference_supplied = [args.reference_species is not None,
                          args.reference_gff is not None, args.reference_genome is not None]
    if any(reference_supplied) and not all(reference_supplied):
        parser.error('--reference-species, --reference-gff and --reference-genome must be supplied together or omitted together')
    if args.reference_species == args.species:
        parser.error('--reference-species must differ from --species')
    if not args.aligner.is_file() or not os.access(args.aligner, os.X_OK):
        parser.error(f'aligner not executable: {args.aligner}')
    global TARGET_SPECIES, FELID_SUPPORT_SPECIES, TARGET_HOMOLOGY_SOURCE
    TARGET_SPECIES = args.species
    TARGET_HOMOLOGY_SOURCE = args.target_homology_source
    source_files = [('EVM', args.evm_gff.resolve())] + read_evidence_manifest(args.evidence_manifest)
    protein_files = parse_protein_refs(args.protein_ref, args.species)
    reference_species = [args.reference_species] if args.reference_species else []
    FELID_SUPPORT_SPECIES = tuple(dict.fromkeys(
        reference_species + [species for species, _ in protein_files if species != 'AME']
    ))
    args.outdir.mkdir(parents=True, exist_ok=True)

    genome_records = read_fasta(args.genome)
    if set(genome_records) != {'ChrY'}:
        raise ValueError(f'Expected a single sequence named ChrY in {args.genome}; found {list(genome_records)}')
    genome = genome_records['ChrY']
    sequence_length = len(genome)
    regions = read_regions(args.regions)
    if not regions or regions[0].start0 != 0 or regions[-1].end0 != sequence_length:
        raise ValueError(f'Primary partition does not span the full {TARGET_SPECIES} ChrY sequence.')

    models: List[Model] = []
    for source, path in source_files:
        if not path.exists():
            raise FileNotFoundError(path)
        parsed = parse_gff(path, source)
        for model in parsed:
            finish_model(model, genome)
            assign_region(model, regions)
        models.extend(parsed)

    refs: List[ReferenceProtein] = []
    if args.reference_gff is not None:
        build_gff_references(refs, args.reference_gff, args.reference_genome, args.reference_species)
    for species, path in protein_files:
        add_reference_fasta(refs, path, species)
    pairs_path = args.outdir / 'protein_alignment_pairs.tsv'
    alignment_path = args.outdir / 'protein_alignment_results.tsv'
    pair_map = build_alignment_pairs(models, refs, pairs_path)
    subprocess.run([str(args.aligner), str(pairs_path), str(alignment_path)], check=True)
    load_alignment_results(alignment_path, pair_map, models, refs)
    for model in models:
        resolve_family_and_support(model)
        preliminary_filter(model)

    clusters = cluster_models(models)
    calls_by_class: Dict[str, List[dict]] = {
        'high_confidence': [], 'provisional': [], 'pseudogene': []}
    rejected_cluster_models: List[Model] = []
    all_calls: List[dict] = []
    for cluster in clusters:
        classification, representative, meta = classify_cluster(cluster)
        family = representative.family or meta['family']
        call = dict(meta)
        call.update({'classification': classification, 'model': representative,
                     'family': family})
        all_calls.append(call)
        if classification in calls_by_class:
            calls_by_class[classification].append(call)
        else:
            rejected_cluster_models.extend(cluster)
    rejected_models = [m for m in models if m.preliminary_reject] + rejected_cluster_models

    write_gff(args.outdir / f'{TARGET_SPECIES}.ChrY.protein_coding.high_confidence.gff3',
              calls_by_class['high_confidence'], 'high_confidence', sequence_length)
    write_gff(args.outdir / f'{TARGET_SPECIES}.ChrY.protein_coding.provisional.gff3',
              calls_by_class['provisional'], 'provisional', sequence_length)
    write_gff(args.outdir / f'{TARGET_SPECIES}.ChrY.pseudogene.gff3',
              calls_by_class['pseudogene'], 'pseudogene', sequence_length)
    write_sequence_files(args.outdir, calls_by_class)
    write_copy_status(args.outdir, calls_by_class)
    write_audit(args.outdir, all_calls, rejected_models)
    evm_models = [m for m in models if m.source == 'EVM']
    write_original_disposition(args.outdir, evm_models, all_calls)
    qc = validate_outputs(args.outdir, sequence_length, calls_by_class)
    write_summary(args.outdir, calls_by_class, models, qc)
    input_files = [args.genome, args.regions, args.evidence_manifest] + [p for _, p in source_files] + [p for _, p in protein_files]
    if args.reference_gff is not None:
        input_files.extend([args.reference_gff, args.reference_genome])
    (args.outdir / f'{TARGET_SPECIES}.ChrY.run_manifest.json').write_text(
        json.dumps({
            'species': TARGET_SPECIES, 'reference_species': args.reference_species,
            'reference_gff': str(args.reference_gff.resolve()) if args.reference_gff else None,
            'reference_genome': str(args.reference_genome.resolve()) if args.reference_genome else None,
            'felid_support_species': FELID_SUPPORT_SPECIES,
            'target_homology_source': TARGET_HOMOLOGY_SOURCE,
            'genome_length': sequence_length,
            'genome_path': str(args.genome.resolve()),
            'regions_path': str(args.regions.resolve()),
            'sources': [(source, str(path.resolve())) for source, path in source_files],
            'protein_references': [(species, str(path)) for species, path in protein_files],
            'input_sha256': {str(p.resolve()): sha256(p.resolve()) for p in input_files},
            'counts': {key: len(value) for key, value in calls_by_class.items()},
        }, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'{TARGET_SPECIES} ChrY curation complete: HC={len(calls_by_class["high_confidence"])} '
          f'PR={len(calls_by_class["provisional"])} PG={len(calls_by_class["pseudogene"])}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
