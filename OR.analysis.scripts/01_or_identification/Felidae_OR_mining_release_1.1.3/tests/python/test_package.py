"""Smoke tests for the installable Python package."""

from typer.testing import CliRunner

from or_mining import __version__
from or_mining.cli import app


def test_version_command() -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0
    assert result.stdout.strip() == __version__
