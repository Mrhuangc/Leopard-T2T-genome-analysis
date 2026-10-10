from pathlib import Path

from or_mining.conda_prewarm import (
    build_nextflow_command,
    discover_environment_files,
    render_prewarm_workflow,
)


def test_conda_prewarm_discovers_and_deduplicates_specs(tmp_path: Path) -> None:
    first = tmp_path / "modules" / "local" / "one" / "environment.yml"
    second = tmp_path / "modules" / "local" / "two" / "environment.yml"
    third = tmp_path / "modules" / "local" / "three" / "environment.yml"
    for path in (first, second, third):
        path.parent.mkdir(parents=True, exist_ok=True)
    first.write_text("dependencies:\n  - python=3.11\n", encoding="utf-8")
    second.write_text(first.read_text(encoding="utf-8"), encoding="utf-8")
    third.write_text("dependencies:\n  - blast=2.16.0\n", encoding="utf-8")

    specs = discover_environment_files(tmp_path)
    workflow = render_prewarm_workflow(specs)

    assert len(specs) == 2
    assert "process PREWARM_ENV_001" in workflow
    assert "process PREWARM_ENV_002" in workflow
    assert "PREWARM_ENV_003" not in workflow
    assert first.resolve().as_posix() in workflow or second.resolve().as_posix() in workflow


def test_conda_prewarm_places_global_config_before_run(tmp_path: Path) -> None:
    command = build_nextflow_command(
        tmp_path / "prewarm.nf",
        tmp_path / "prewarm.config",
        tmp_path / "cache",
        tmp_path / "info",
    )

    assert command.index("-C") < command.index("-log") < command.index("run")
    assert command.index("run") < command.index(str(tmp_path / "prewarm.nf"))
