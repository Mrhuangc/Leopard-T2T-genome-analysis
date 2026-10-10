"""Fail-fast validation of a new species run before cluster jobs are submitted."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import yaml

from or_mining.blast import _open_text, parse_blast_row, validate_merged_file
from or_mining.or_clade_classification import infer_reference_clade
from or_mining.species_code import build_species_code, validate_species_code

REQUIRED_PATHS = (
    "genome",
    "or_db",
    "gpcr_db",
    "template_alignment",
    "known_type1_ids",
    "known_type2_ids",
    "or_clade_db",
    "interproscan_dir",
)


class UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that refuses silently overwritten parameter keys."""


def _construct_unique_mapping(
    loader: UniqueKeyLoader, node: yaml.nodes.MappingNode, deep: bool = False
) -> dict[Any, Any]:
    loader.flatten_mapping(node)
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise ValueError(
                f"duplicate YAML key {key!r} at line {key_node.start_mark.line + 1}"
            )
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def load_params(path: Path) -> dict[str, Any]:
    data = yaml.load(path.read_text(encoding="utf-8"), Loader=UniqueKeyLoader)
    if not isinstance(data, dict):
        raise TypeError("parameter file must contain a YAML mapping")
    return data


def _fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("rt", encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            if raw.startswith(">"):
                identifier = raw[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"empty FASTA identifier in {path}")
                if identifier in ids:
                    raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
                ids.add(identifier)
    if not ids:
        raise ValueError(f"no FASTA records in {path}")
    return ids


def _fai_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed FASTA index row")
            if fields[0] in lengths:
                raise ValueError(f"duplicate contig {fields[0]!r} in {path}")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError(f"empty FASTA index: {path}")
    return lengths


def _id_file(path: Path) -> set[str]:
    values = {
        line.strip().split()[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    if not values:
        raise ValueError(f"empty identifier file: {path}")
    return values


def _add(checks: list[dict[str, Any]], name: str, status: str, detail: str) -> None:
    checks.append({"name": name, "status": status, "detail": detail})


def run_preflight(params_file: Path) -> dict[str, Any]:
    params_file = params_file.resolve()
    params = load_params(params_file)
    checks: list[dict[str, Any]] = []

    species = str(params.get("species") or "")
    try:
        validate_species_code(species)
        expected = build_species_code(
            str(params.get("taxon_genus") or ""),
            str(params.get("taxon_species") or ""),
            str(params.get("taxon_subspecies") or "") or None,
        )
        if species != expected:
            raise ValueError(f"species={species!r}, but Latin-name-derived code is {expected!r}")
        _add(checks, "species_code", "pass", f"{species} matches the Latin taxon fields")
    except ValueError as exc:
        _add(checks, "species_code", "fail", str(exc))

    assembly_id = str(params.get("assembly_id") or "").strip()
    _add(
        checks,
        "assembly_id",
        "pass" if assembly_id else "fail",
        assembly_id or "assembly_id is required and must be versioned",
    )

    initial_mode = str(params.get("initial_search_mode") or "reuse").strip().lower()
    _add(
        checks,
        "initial_search_mode",
        "pass" if initial_mode in {"reuse", "generate", "auto"} else "fail",
        initial_mode or "must be reuse, generate, or auto",
    )

    resolved: dict[str, Path] = {}
    for key in REQUIRED_PATHS:
        value = params.get(key)
        if not value:
            _add(checks, f"path:{key}", "fail", "parameter is missing")
            continue
        candidate = Path(str(value)).expanduser()
        resolved[key] = candidate
        expected_kind = "directory" if key == "interproscan_dir" else "file"
        exists = candidate.is_dir() if expected_kind == "directory" else candidate.is_file()
        _add(
            checks,
            f"path:{key}",
            "pass" if exists else "fail",
            f"{expected_kind}: {candidate}",
        )

    try:
        # Local import avoids a module cycle: initial_search reuses load_params.
        from or_mining.initial_search import resolve_merged_input

        merged_path, merged_source = resolve_merged_input(params)
        if merged_path.is_file():
            resolved["input_merged"] = merged_path
            _add(
                checks,
                "path:input_merged",
                "pass",
                f"{merged_source} file: {merged_path}",
            )
        elif merged_source == "generated" and initial_mode in {"generate", "auto"}:
            _add(
                checks,
                "path:input_merged",
                "warning",
                f"will be generated before stage 1: {merged_path}",
            )
        else:
            _add(checks, "path:input_merged", "fail", f"file: {merged_path}")
    except ValueError as exc:
        _add(checks, "path:input_merged", "fail", str(exc))

    genome = resolved.get("genome")
    fai_lengths: dict[str, int] = {}
    if genome and genome.is_file():
        fai = Path(str(genome) + ".fai")
        if fai.is_file():
            try:
                fai_lengths = _fai_lengths(fai)
                _add(checks, "genome_fai", "pass", f"{len(fai_lengths)} indexed sequences")
            except (ValueError, OSError) as exc:
                _add(checks, "genome_fai", "fail", str(exc))
        else:
            _add(checks, "genome_fai", "fail", f"missing {fai}; run samtools faidx")

    merged = resolved.get("input_merged")
    merged_subject_lengths: dict[str, int] = {}
    merged_query_ids: set[str] = set()
    if merged and merged.is_file():
        validation = validate_merged_file(merged)
        _add(
            checks,
            "merged_schema",
            validation["status"],
            f"{validation['counts']['valid_rows']} valid / "
            f"{validation['counts']['invalid_rows']} invalid rows",
        )
        if validation["status"] == "pass":
            with _open_text(merged) as handle:
                for raw in handle:
                    hit = parse_blast_row(raw)
                    merged_query_ids.add(hit.qseqid)
                    merged_subject_lengths[hit.sseqid] = hit.slen

    if fai_lengths and merged_subject_lengths:
        absent = sorted(set(merged_subject_lengths) - set(fai_lengths))
        conflicts = sorted(
            contig
            for contig, length in merged_subject_lengths.items()
            if contig in fai_lengths and fai_lengths[contig] != length
        )
        if absent or conflicts:
            _add(
                checks,
                "merged_genome_contract",
                "fail",
                f"missing_contigs={absent[:5]}; length_conflicts={conflicts[:5]}",
            )
        else:
            _add(
                checks,
                "merged_genome_contract",
                "pass",
                f"{len(merged_subject_lengths)} subjects match the genome index exactly",
            )

    fasta_sets: dict[str, set[str]] = {}
    for key in ("or_db", "gpcr_db", "template_alignment", "or_clade_db"):
        path = resolved.get(key)
        if path and path.is_file():
            try:
                fasta_sets[key] = _fasta_ids(path)
                _add(checks, f"fasta:{key}", "pass", f"{len(fasta_sets[key])} unique records")
            except (ValueError, OSError, UnicodeError) as exc:
                _add(checks, f"fasta:{key}", "fail", str(exc))

    if merged_query_ids and "or_db" in fasta_sets:
        missing_queries = sorted(merged_query_ids - fasta_sets["or_db"])
        _add(
            checks,
            "merged_query_reference_contract",
            "warning" if missing_queries else "pass",
            (
                f"{len(missing_queries)} merged queries are absent from or_db; "
                f"examples={missing_queries[:5]}"
                if missing_queries
                else f"all {len(merged_query_ids)} merged queries occur in or_db"
            ),
        )

    clade_reference_ids = fasta_sets.get("or_clade_db", set())
    if clade_reference_ids:
        observed_clades = sorted(
            {clade for identifier in clade_reference_ids if (clade := infer_reference_clade(identifier))}
        )
        expected_clades = {
            "alpha",
            "beta",
            "gamma",
            "delta",
            "epsilon",
            "eta",
            "zeta",
            "kappa",
            "lambda",
            "theta",
        }
        missing_clades = sorted(expected_clades - set(observed_clades))
        _add(
            checks,
            "or_clade_reference_labels",
            "fail" if missing_clades else "pass",
            f"observed={observed_clades}; missing={missing_clades}",
        )

    template_ids = fasta_sets.get("template_alignment", set())
    for key in ("known_type1_ids", "known_type2_ids"):
        path = resolved.get(key)
        if path and path.is_file() and template_ids:
            try:
                ids = _id_file(path)
                missing = sorted(ids - template_ids)
                _add(
                    checks,
                    f"template_membership:{key}",
                    "fail" if missing else "pass",
                    f"{len(ids)} IDs; missing_from_template={missing[:5]}",
                )
            except (ValueError, OSError) as exc:
                _add(checks, f"template_membership:{key}", "fail", str(exc))

    gpcr = resolved.get("gpcr_db")
    if gpcr and gpcr.is_file():
        text = gpcr.read_text(encoding="utf-8", errors="ignore")
        _add(
            checks,
            "gpcr_or_labels",
            "pass" if "OR-Receptor" in text else "fail",
            "database must contain the literal label OR-Receptor used by the decision rule",
        )

    ipr = resolved.get("interproscan_dir")
    if ipr and ipr.is_dir():
        executable = ipr / "interproscan.sh"
        _add(
            checks,
            "interproscan_executable",
            "pass" if executable.is_file() else "fail",
            str(executable),
        )

    numeric_rules = {
        "max_cpus": lambda x: x >= 1,
        "max_memory_gb": lambda x: x >= 4,
        "max_time_h": lambda x: x >= 1,
        "candidate_min_width": lambda x: x >= 1,
        "candidate_flank": lambda x: x >= 0,
        "min_orf_nt": lambda x: x >= 3 and x % 3 == 0,
        "maximum_intron_length": lambda x: x >= 50,
        "exonerate_top_queries_per_region": lambda x: x >= 1,
        "exonerate_pilot_pairs": lambda x: x >= 1,
        "exonerate_soft_timeout_minutes": lambda x: x >= 1,
        "strict_refinement_timeout_minutes": lambda x: x >= 1,
        "or_operational_cluster_gap_bp": lambda x: x >= 0,
        "initial_query_shards": lambda x: 1 <= x <= 16,
        "initial_tblastn_cpus": lambda x: x >= 1,
        "initial_tblastn_max_forks": lambda x: 1 <= x <= 16,
        "initial_tblastn_timeout_hours": lambda x: x >= 1,
    }
    initial_defaults = {
        "initial_query_shards": 2,
        "initial_tblastn_cpus": 8,
        "initial_tblastn_max_forks": 2,
        "initial_tblastn_timeout_hours": 24,
    }
    for key, rule in numeric_rules.items():
        try:
            value = int(params.get(key, initial_defaults.get(key)))
            valid = rule(value)
        except (KeyError, TypeError, ValueError):
            value = params.get(key)
            valid = False
        _add(checks, f"parameter:{key}", "pass" if valid else "fail", repr(value))

    try:
        initial_evalue = float(params.get("initial_tblastn_evalue", 1e-5))
        valid_initial_evalue = 0.0 < initial_evalue <= 1.0
    except (KeyError, TypeError, ValueError):
        initial_evalue = params.get("initial_tblastn_evalue")
        valid_initial_evalue = False
    _add(
        checks,
        "parameter:initial_tblastn_evalue",
        "pass" if valid_initial_evalue else "fail",
        repr(initial_evalue),
    )

    try:
        initial_cpus = int(params.get("initial_tblastn_cpus", 8))
        maximum_cpus = int(params["max_cpus"])
        resource_valid = initial_cpus <= maximum_cpus
    except (KeyError, TypeError, ValueError):
        initial_cpus = params.get("initial_tblastn_cpus")
        maximum_cpus = params.get("max_cpus")
        resource_valid = False
    _add(
        checks,
        "initial_search_cpu_contract",
        "pass" if resource_valid else "fail",
        f"initial_tblastn_cpus={initial_cpus}; max_cpus={maximum_cpus}",
    )

    outdir = Path(str(params.get("outdir") or ""))
    foreign = []
    if outdir.is_dir() and species:
        foreign = sorted(
            path.name
            for path in outdir.rglob("*.summary.json")
            if not path.name.startswith(species + ".")
        )
    _add(
        checks,
        "output_species_isolation",
        "fail" if foreign else "pass",
        f"foreign summaries={foreign[:5]}" if foreign else str(outdir),
    )
    _add(
        checks,
        "output_directory_name",
        "pass" if outdir.name == species else "fail",
        f"outdir basename must equal species code: {outdir.name!r} versus {species!r}",
    )

    publish_mode = str(params.get("publish_mode") or "link").strip()
    _add(
        checks,
        "publish_mode",
        "pass" if publish_mode in {"copy", "link"} else "fail",
        (
            f"{publish_mode}; link avoids duplicate payload I/O on one filesystem"
            if publish_mode in {"copy", "link"}
            else "must be 'link' (same filesystem) or 'copy' (cross-filesystem)"
        ),
    )

    for key in ("slurm_partition", "slurm_account", "slurm_qos"):
        value = str(params.get(key) or "").strip()
        _add(
            checks,
            f"cluster:{key}",
            "pass" if value and value != "REPLACE_ME" else "warning",
            value or "not set; acceptable only when the cluster does not require it",
        )

    counts = Counter(item["status"] for item in checks)
    return {
        "schema_version": "1.0",
        "status": "fail" if counts["fail"] else "pass_with_warnings" if counts["warning"] else "pass",
        "species": species,
        "assembly_id": assembly_id,
        "params_file": str(params_file),
        "counts": dict(sorted(counts.items())),
        "checks": checks,
    }


def write_preflight(report: dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
