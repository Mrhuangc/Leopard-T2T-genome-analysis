"""Build a clean, checksummed release archive without run-time artifacts."""

from __future__ import annotations

import hashlib
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT.parent / f"{ROOT.name}.tar.gz"
ARCHIVE_CHECKSUM = ROOT.parent / f"{ROOT.name}.tar.gz.sha256"
MANIFEST = ROOT / "MANIFEST.sha256"
EXCLUDED_DIRS = {
    ".conda",
    ".git",
    ".nextflow",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "pipeline_info",
    "results",
    "work",
}
EXCLUDED_FILES = {"MANIFEST.sha256"}


def clean_caches() -> None:
    for path in sorted(ROOT.rglob("__pycache__"), reverse=True):
        shutil.rmtree(path)
    for path in ROOT.rglob("*.py[co]"):
        path.unlink()


def release_files(*, include_manifest: bool = False) -> list[Path]:
    files: list[Path] = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if any(part in EXCLUDED_DIRS for part in relative.parts):
            continue
        if not path.is_file():
            continue
        if not include_manifest and relative.as_posix() in EXCLUDED_FILES:
            continue
        files.append(path)
    return sorted(files, key=lambda item: item.relative_to(ROOT).as_posix())


def write_manifest() -> None:
    lines = []
    for path in release_files():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def normalized_tar_info(info: tarfile.TarInfo) -> tarfile.TarInfo:
    info.uid = 0
    info.gid = 0
    info.uname = "root"
    info.gname = "root"
    info.mode = 0o755 if info.name.endswith(".sh") else 0o644
    return info


def build_archive() -> None:
    with tarfile.open(ARCHIVE, "w:gz", compresslevel=6) as handle:
        for path in release_files(include_manifest=True):
            relative = path.relative_to(ROOT)
            archive_name = Path(ROOT.name) / relative
            handle.add(path, arcname=archive_name.as_posix(), filter=normalized_tar_info)


def main() -> None:
    clean_caches()
    write_manifest()
    build_archive()
    archive_digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    ARCHIVE_CHECKSUM.write_text(
        f"{archive_digest}  {ARCHIVE.name}\n", encoding="utf-8", newline="\n"
    )
    print(f"release_root={ROOT}")
    print(f"files={len(release_files(include_manifest=True))}")
    print(f"archive={ARCHIVE}")
    print(f"archive_sha256={archive_digest}")
    print(f"checksum={ARCHIVE_CHECKSUM}")


if __name__ == "__main__":
    main()
