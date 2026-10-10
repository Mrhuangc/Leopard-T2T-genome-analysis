"""Tests for audited phylogenetic OR-clade filtering."""

from pathlib import Path

from or_mining.tree_filter import filter_or_clade


def test_tree_filter_reproduces_nc_rule_and_reports_outgroup_audit(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    template = tmp_path / "template.fa"
    type1 = tmp_path / "type1.txt"
    type2 = tmp_path / "type2.txt"
    proteins = tmp_path / "candidates.prot.fa"
    nucleotides = tmp_path / "candidates.nt.fa"
    coordinates = tmp_path / "candidates.tsv"
    prefix = tmp_path / "sp"

    tree.write_text(
        "((out1:1,out2:1):1,(((known1:1,cand1:1):1,(known2:1,cand2:1):1):1,"
        "cand_out:1):1);\n",
        encoding="utf-8",
    )
    template.write_text(
        ">out1\nAAAA\n>out2\nAAAA\n>known1\nAAAA\n>known2\nAAAA\n",
        encoding="utf-8",
    )
    type1.write_text("known1\n", encoding="utf-8")
    type2.write_text("known2\n", encoding="utf-8")
    proteins.write_text(
        ">cand1\nAAAA\n>cand2\nAAAA\n>cand_out\nAAAA\n", encoding="utf-8"
    )
    nucleotides.write_text(
        ">cand1\nGCTGCTGCTGCT\n>cand2\nGCTGCTGCTGCT\n"
        ">cand_out\nGCTGCTGCTGCT\n",
        encoding="utf-8",
    )
    coordinates.write_text(
        "orf_id\tsseqid\tstart_1\tend_1\tstrand\n"
        "cand1\tChr1\t1\t12\t+\n"
        "cand2\tChr1\t21\t32\t+\n"
        "cand_out\tChr1\t41\t52\t+\n",
        encoding="utf-8",
    )

    summary = filter_or_clade(
        tree,
        template,
        type1,
        type2,
        proteins,
        nucleotides,
        coordinates,
        prefix,
    )

    assert summary["algorithms_agree"] is False
    assert summary["counts"]["retained_candidates"] == 3
    assert summary["counts"]["rejected_candidates"] == 0
    assert summary["counts"]["explicit_outgroup_retained_candidates"] == 2
    assert summary["legacy_only_ids"] == ["cand_out"]
    assert summary["explicit_outgroup_only_ids"] == []
    assert Path(f"{prefix}.tree_filtered.ids.txt").read_text(encoding="utf-8") == (
        "cand1\ncand2\ncand_out\n"
    )


def test_tree_filter_allows_non_monophyletic_audit_outgroups(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    template = tmp_path / "template.fa"
    type1 = tmp_path / "type1.txt"
    type2 = tmp_path / "type2.txt"
    proteins = tmp_path / "candidates.prot.fa"
    nucleotides = tmp_path / "candidates.nt.fa"
    coordinates = tmp_path / "candidates.tsv"
    prefix = tmp_path / "sp"

    tree.write_text(
        "((out1:1,cand_out:1):1,((known1:1,cand1:1):1,(known2:1,out2:1):1):1);\n",
        encoding="utf-8",
    )
    template.write_text(
        ">out1\nAAAA\n>out2\nAAAA\n>known1\nAAAA\n>known2\nAAAA\n",
        encoding="utf-8",
    )
    type1.write_text("known1\n", encoding="utf-8")
    type2.write_text("known2\n", encoding="utf-8")
    proteins.write_text(">cand1\nAAAA\n>cand_out\nAAAA\n", encoding="utf-8")
    nucleotides.write_text(
        ">cand1\nGCTGCTGCTGCT\n>cand_out\nGCTGCTGCTGCT\n",
        encoding="utf-8",
    )
    coordinates.write_text(
        "orf_id\tsseqid\tstart_1\tend_1\tstrand\n"
        "cand1\tChr1\t1\t12\t+\n"
        "cand_out\tChr1\t21\t32\t+\n",
        encoding="utf-8",
    )

    summary = filter_or_clade(
        tree,
        template,
        type1,
        type2,
        proteins,
        nucleotides,
        coordinates,
        prefix,
    )

    assert summary["status"] == "pass"
    assert summary["audit_status"] == "not_applicable_non_monophyletic_outgroups"
    assert summary["algorithms_agree"] is None
    assert summary["counts"]["explicit_outgroup_retained_candidates"] is None
    assert Path(f"{prefix}.nc_legacy.rooted.treefile").is_file()
