import pytest

from or_mining.species_code import build_species_code, validate_species_code


def test_species_code_binomial_and_trinomial() -> None:
    assert build_species_code("Panthera", "pardus") == "Ppa"
    assert build_species_code("Canis", "lupus", "familiaris") == "Clu_fa"


@pytest.mark.parametrize("code", ["Ppa_T2T", "ppa", "Clup_fa", "Clu_FA", "Clu_fam"])
def test_invalid_species_codes(code: str) -> None:
    with pytest.raises(ValueError):
        validate_species_code(code)


def test_short_or_non_latin_taxon_rejected() -> None:
    with pytest.raises(ValueError):
        build_species_code("P", "x")
    with pytest.raises(ValueError):
        build_species_code("犬", "lupus")
