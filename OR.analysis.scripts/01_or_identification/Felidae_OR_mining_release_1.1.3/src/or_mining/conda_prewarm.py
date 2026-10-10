"""Prepare every process Conda environment before Slurm jobs are submitted."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any


def discover_environment_files(root: Path) -> list[Path]:
    """Return one representative for every distinct production environment spec."""

    module_root = root.resolve() / "modules" / "local"
    if not module_root.is_dir():
        raise FileNotFoundError(f"module directory not found: {module_root}")

    unique: dict[str, Path] = {}
    for path in sorted(module_root.rglob("environment.yml")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        unique.setdefault(digest, path.resolve())
    if not unique:
        raise FileNotFoundError(f"no process environment.yml files found below {module_root}")
    return [unique[digest] for digest in sorted(unique)]


def render_prewarm_workflow(environment_files: list[Path]) -> str:
    """Render a tiny workflow whose only purpose is to materialize Conda caches."""

    lines = ["nextflow.enable.dsl=2", ""]
    process_names: list[str] = []
    for index, path in enumerate(environment_files, start=1):
        name = f"PREWARM_ENV_{index:03d}"
        process_names.append(name)
        conda_path = path.as_posix().replace("'", "\\'")
        lines.extend(
            [
                f"process {name} {{",
                f"    tag '{path.parent.name}'",
                f"    conda '{conda_path}'",
                "",
                "    output:",
                f"    path '{name}.ok'",
                "",
                "    script:",
                '    """',
                f"    printf 'status=ok\\n' > {name}.ok",
                '    """',
                "}",
                "",
            ]
        )

    lines.append("workflow {")
    lines.extend(f"    {name}()" for name in process_names)
    lines.extend(["}", ""])
    return "\n".join(lines)


def build_nextflow_command(
    workflow_path: Path,
    config_path: Path,
    cache_dir: Path,
    info_dir: Path,
) -> list[str]:
    """Build the CLI with global options before the ``run`` subcommand."""

    return [
        "nextflow",
        "-C",
        str(config_path),
        "-log",
        str(info_dir / "conda_prewarm.nextflow.log"),
        "run",
        str(workflow_path),
        "-w",
        str(cache_dir / ".prewarm-work"),
        "-resume",
    ]


def prewarm_conda_environments(
    root: Path,
    cache_dir: Path,
    *,
    generate_only: bool = False,
) -> dict[str, Any]:
    """Create all environments on the connected controller/login node."""

    root = root.resolve()
    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    info_dir = root / "pipeline_info"
    info_dir.mkdir(parents=True, exist_ok=True)

    environment_files = discover_environment_files(root)
    workflow_path = info_dir / "conda_prewarm.generated.nf"
    workflow_path.write_text(
        render_prewarm_workflow(environment_files), encoding="utf-8", newline="\n"
    )

    payload: dict[str, Any] = {
        "schema_version": "1.0",
        "status": "generated" if generate_only else "pending",
        "cache_dir": str(cache_dir),
        "environment_specs": [str(path) for path in environment_files],
        "environment_spec_count": len(environment_files),
    }
    summary_path = info_dir / "conda_prewarm.summary.json"

    if not generate_only:
        config_path = root / "conf" / "conda_prewarm.config"
        if not config_path.is_file():
            raise FileNotFoundError(f"prewarm config not found: {config_path}")
        environment = os.environ.copy()
        environment["NXF_CONDA_CACHEDIR"] = str(cache_dir)
        environment.pop("NXF_OFFLINE", None)
        environment.pop("CONDA_OFFLINE", None)
        environment.pop("MAMBA_OFFLINE", None)
        command = build_nextflow_command(workflow_path, config_path, cache_dir, info_dir)
        completed = subprocess.run(command, cwd=root, env=environment, check=False)
        if completed.returncode != 0:
            payload["status"] = "fail"
            payload["nextflow_exit_code"] = completed.returncode
            summary_path.write_text(
                json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
            )
            raise RuntimeError(
                "Conda prewarm failed on the login node; formal Slurm jobs were not submitted. "
                f"See {info_dir / 'conda_prewarm.nextflow.log'}"
            )
        payload["status"] = "pass"
        payload["nextflow_exit_code"] = completed.returncode

    summary_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return payload
