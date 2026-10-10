"""Stable taxon codes used in catalog identifiers and output names."""

from __future__ import annotations

import re

SPECIES_CODE_RE = re.compile(r"^[A-Z][a-z]{2}(?:_[a-z]{2})?$")
LATIN_WORD_RE = re.compile(r"^[A-Za-z]+$")


def build_species_code(genus: str, species: str, subspecies: str | None = None) -> str:
    """Return ``Gss`` or ``Gss_uu`` from a Latin binomial/trinomial.

    Examples: ``Panthera pardus -> Ppa`` and
    ``Canis lupus familiaris -> Clu_fa``.
    """

    genus = genus.strip()
    species = species.strip()
    subspecies = subspecies.strip() if subspecies else None
    for label, value, minimum in (
        ("genus", genus, 1),
        ("species", species, 2),
        ("subspecies", subspecies, 2),
    ):
        if value is None:
            continue
        if len(value) < minimum or LATIN_WORD_RE.fullmatch(value) is None:
            raise ValueError(
                f"{label} must contain at least {minimum} unaccented Latin letters: {value!r}"
            )
    code = genus[0].upper() + species[:2].lower()
    if subspecies:
        code += "_" + subspecies[:2].lower()
    return code


def validate_species_code(code: str) -> None:
    """Raise when a code is not in the release naming convention."""

    if SPECIES_CODE_RE.fullmatch(code) is None:
        raise ValueError(
            "species must match Gss or Gss_uu (for example Ppa or Clu_fa); "
            "assembly labels such as T2T belong in assembly_id, not species"
        )
