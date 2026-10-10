"""Command-line entry point for the Felidae OR mining package."""

from __future__ import annotations

import json
from pathlib import Path

import typer

from or_mining import __version__
from or_mining.blast import validate_merged_file, write_validation_report
from or_mining.conda_prewarm import prewarm_conda_environments
from or_mining.initial_search import resolve_merged_input, write_resolved_params
from or_mining.preflight import load_params, run_preflight, write_preflight
from or_mining.species_code import build_species_code
from or_mining.stage_gate import check_stage

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Utilities used by the Felidae OR mining workflow.",
)


@app.callback()
def main() -> None:
    """Felidae OR mining command group."""


@app.command()
def version() -> None:
    """Print the installed package version."""

    typer.echo(__version__)


@app.command("prewarm-conda")
def prewarm_conda_command(
    root: Path = typer.Argument(..., exists=True, file_okay=False, resolve_path=True),
    cache_dir: Path = typer.Argument(..., file_okay=False, resolve_path=True),
    generate_only: bool = typer.Option(False, "--generate-only"),
) -> None:
    """Create every process environment before jobs enter offline compute nodes."""

    try:
        payload = prewarm_conda_environments(root, cache_dir, generate_only=generate_only)
    except (FileNotFoundError, RuntimeError) as exc:
        typer.echo(f"FAIL: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(
        f"status={payload['status']} environments={payload['environment_spec_count']} "
        f"cache={payload['cache_dir']}"
    )


@app.command("species-code")
def species_code_command(
    genus: str = typer.Argument(..., help="Latin genus, for example Canis."),
    species: str = typer.Argument(..., help="Latin species epithet, for example lupus."),
    subspecies: str | None = typer.Option(None, "--subspecies", "-u"),
) -> None:
    """Generate the required Gss or Gss_uu taxon code."""

    typer.echo(build_species_code(genus, species, subspecies))


@app.command("preflight")
def preflight_command(
    params_file: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True),
    report: Path | None = typer.Option(None, "--report"),
) -> None:
    """Validate paths, naming, references, merged hits and the genome index."""

    try:
        payload = run_preflight(params_file)
    except (OSError, TypeError, ValueError) as exc:
        typer.echo(f"FAIL: invalid parameter file: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    destination = report or Path(str(payload["species"]) + ".preflight.summary.json")
    write_preflight(payload, destination)
    typer.echo(
        f"status={payload['status']} pass={payload['counts'].get('pass', 0)} "
        f"warning={payload['counts'].get('warning', 0)} "
        f"fail={payload['counts'].get('fail', 0)} report={destination}"
    )
    if payload["status"] == "fail":
        raise typer.Exit(code=2)


@app.command("check-stage")
def check_stage_command(
    params_file: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True),
    stage: int = typer.Argument(..., min=0, max=19),
) -> None:
    """Require the primary summary for a completed production stage to pass."""

    try:
        path, payload = check_stage(params_file, stage)
    except (ValueError, FileNotFoundError, json.JSONDecodeError) as exc:
        typer.echo(f"FAIL stage={stage}: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"PASS stage={stage} status={payload['status']} summary={path}")


@app.command("resolve-input")
def resolve_input_command(
    params_file: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True),
    require_exists: bool = typer.Option(False, "--require-exists"),
) -> None:
    """Print the effective supplied or generated merged TBLASTN path."""

    try:
        params = load_params(params_file)
        path, source = resolve_merged_input(params)
    except (OSError, TypeError, ValueError) as exc:
        typer.echo(f"FAIL: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    if require_exists and not path.is_file():
        typer.echo(f"FAIL: resolved merged input is missing: {path}", err=True)
        raise typer.Exit(code=2)
    typer.echo(f"source={source}\tpath={path.resolve()}")


@app.command("write-resolved-params")
def write_resolved_params_command(
    params_file: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True),
    output: Path = typer.Option(..., "--output"),
) -> None:
    """Snapshot parameters with the effective merged input made explicit."""

    try:
        params = write_resolved_params(params_file, output)
    except (OSError, TypeError, ValueError) as exc:
        typer.echo(f"FAIL: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(
        f"status=pass source={params['resolved_initial_input_source']} "
        f"input_merged={params['input_merged']} output={output}"
    )


@app.command("validate-merged")
def validate_merged_command(
    input_merged: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Headerless 16-column merged TBLASTN file.",
    ),
    report: Path = typer.Option(
        Path("merged.validation.json"),
        "--report",
        help="Destination JSON audit report.",
    ),
    max_errors: int = typer.Option(20, min=1, help="Maximum errors included in the report."),
) -> None:
    """Validate the merged TBLASTN table before OR extraction."""

    result = validate_merged_file(input_merged, max_reported_errors=max_errors)
    write_validation_report(result, report)
    counts = result["counts"]
    typer.echo(
        f"status={result['status']} total={counts['total_rows']} "
        f"valid={counts['valid_rows']} invalid={counts['invalid_rows']} report={report}"
    )
    if result["status"] != "pass":
        raise typer.Exit(code=2)
