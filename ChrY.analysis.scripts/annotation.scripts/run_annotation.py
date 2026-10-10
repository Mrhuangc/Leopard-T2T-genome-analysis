#!/usr/bin/env python3
"""Compile the protein aligner, curate one species, and independently validate it."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def required_path(config_dir: Path, value: str, label: str) -> Path:
    if not value:
        raise ValueError(f'{label} is required')
    path = (config_dir / value).resolve()
    if not path.is_file():
        raise FileNotFoundError(f'{label}: {path}')
    return path


def fasta_length(path: Path) -> int:
    headers = 0
    length = 0
    for line in path.open(encoding='utf-8'):
        if line.startswith('>'):
            headers += 1
            if line[1:].split()[0] != 'ChrY':
                raise ValueError(f'{path}: FASTA record must be ChrY')
        else:
            length += len(line.strip())
    if headers != 1:
        raise ValueError(f'{path}: expected one ChrY FASTA record, got {headers}')
    return length


def run(config_path: Path) -> None:
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text(encoding='utf-8'))
    base = config_path.parent
    species = config['species']
    if species not in ('PPA', 'FCA'):
        raise ValueError('species must be PPA or FCA')
    # An already available reference annotation may supplement the run, but
    # the target species must never require the other target's new output.
    reference_keys = ('reference_species', 'reference_gff', 'reference_genome')
    supplied = [bool(config.get(key)) for key in reference_keys]
    if any(supplied) and not all(supplied):
        raise ValueError('reference_species, reference_gff and reference_genome must be supplied together or omitted together')
    reference_args = []
    if all(supplied):
        reference = config['reference_species']
        if reference not in ('PPA', 'FCA') or reference == species:
            raise ValueError('reference_species must be the other felid species')
        reference_gff = required_path(base, config['reference_gff'], 'reference_gff')
        reference_genome = required_path(base, config['reference_genome'], 'reference_genome')
        reference_args = ['--reference-species', reference,
                          '--reference-gff', str(reference_gff),
                          '--reference-genome', str(reference_genome)]
    genome = required_path(base, config['genome'], 'genome')
    evm = required_path(base, config['evm_gff'], 'evm_gff')
    evidence = required_path(base, config['evidence_manifest'], 'evidence_manifest')
    regions = required_path(base, config['regions'], 'regions')
    expected = int(config['expected_length'])
    observed = fasta_length(genome)
    if observed != expected:
        raise ValueError(f'{species} ChrY length {observed}, expected {expected}')
    proteins = []
    for name, value in config.get('protein_references', {}).items():
        proteins.extend(['--protein-ref', f'{name.upper()}={required_path(base, value, name)}'])
    outdir = (base / config['outdir']).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent / 'scripts'
    aligner = outdir / 'protein_align_batch'
    subprocess.run(['cc', '-O3', '-std=c11', str(scripts / 'protein_align_batch.c'),
                    '-o', str(aligner)], check=True)
    subprocess.run([
        sys.executable, str(scripts / 'curate_felid_chry_annotation.py'),
        '--species', species, '--genome', str(genome), '--evm-gff', str(evm),
        '--evidence-manifest', str(evidence), '--regions', str(regions),
        *reference_args, '--target-homology-source',
        config.get('target_homology_source', 'GeMoMa_gp'), '--aligner', str(aligner),
        '--outdir', str(outdir), *proteins,
    ], check=True)
    subprocess.run([
        sys.executable, str(scripts / 'validate_felid_chry_annotation.py'),
        '--species', species, '--genome', str(genome), '--regions', str(regions),
        '--results', str(outdir), '--report', str(outdir / f'{species}.ChrY.file_validation.tsv'),
    ], check=True)
    print(f'{species}: annotation and independent validation complete: {outdir}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    args = parser.parse_args()
    run(args.config)
