import tomllib
from pathlib import Path

from or_mining import __version__


def test_release_metadata_has_only_imported_runtime_dependencies() -> None:
    project = Path(__file__).resolve().parents[2]
    payload = tomllib.loads((project / "pyproject.toml").read_text(encoding="utf-8"))
    metadata = payload["project"]
    assert metadata["version"] == __version__ == "1.1.3"
    dependencies = metadata["dependencies"]
    assert dependencies == [
        "biopython>=1.85,<2",
        "pyyaml>=6,<7",
        "typer>=0.12,<1",
    ]
    forbidden = {"numpy", "polars", "pydantic", "pyfaidx", "intervaltree", "jsonschema"}
    assert not any(
        dependency.split("<", 1)[0].split(">", 1)[0].split("=", 1)[0].lower() in forbidden
        for dependency in dependencies
    )
