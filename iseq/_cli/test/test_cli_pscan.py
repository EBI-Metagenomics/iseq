import os
from filecmp import cmp

from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath
from iseq.file import diff


def test_cli_pscan_GALNBKIG_00914_ont_01_plus_strand(tmp_path):
    # Related to https://iseq-py.s3.eu-west-2.amazonaws.com/Report+on+NMM+vs+tBLASTn.pdf
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    profile = example_filepath("PF00113.hmm")
    fasta = example_filepath("GALNBKIG_00914_ont_01_plus_strand.fasta")
    oamino = example_filepath("PF00113_GALNBKIG_00914_ont_01_plus_strand_oamino.fasta")
    ocodon = example_filepath("PF00113_GALNBKIG_00914_ont_01_plus_strand_ocodon.fasta")
    output = example_filepath("PF00113_GALNBKIG_00914_ont_01_plus_strand_output.gff")
    r = invoke(cli, ["pscan", str(profile), str(fasta)])
    assert r.exit_code == 0, r.output
    assert cmp(oamino, "oamino.fasta", shallow=False), diff(oamino, "oamino.fasta")
    assert cmp(ocodon, "ocodon.fasta", shallow=False), diff(ocodon, "ocodon.fasta")
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")


def test_cli_pscan_GALNBKIG_00914_ont_01_plus_strand_window(tmp_path):
    # Related to https://iseq-py.s3.eu-west-2.amazonaws.com/Report+on+NMM+vs+tBLASTn.pdf
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    profile = example_filepath("PF00113.hmm")
    fasta = example_filepath("GALNBKIG_00914_ont_01_plus_strand.fasta")
    oamino = example_filepath(
        "PF00113_GALNBKIG_00914_ont_01_plus_strand_oamino_window.fasta"
    )
    ocodon = example_filepath(
        "PF00113_GALNBKIG_00914_ont_01_plus_strand_ocodon_window.fasta"
    )
    output = example_filepath(
        "PF00113_GALNBKIG_00914_ont_01_plus_strand_output_window.gff"
    )
    r = invoke(cli, ["pscan", str(profile), str(fasta), "--window", "-1"])
    assert r.exit_code == 0, r.output
    assert cmp(oamino, "oamino.fasta", shallow=False), diff(oamino, "oamino.fasta")
    assert cmp(ocodon, "ocodon.fasta", shallow=False), diff(ocodon, "ocodon.fasta")
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")


def test_cli_pscan_nofile_output(tmp_path, GALNBKIG_cut):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    PF03373 = example_filepath("PF03373.hmm")
    r = invoke(cli, ["pscan", str(PF03373), str(fasta)])
    assert r.exit_code == 0, r.output


def test_cli_pscan_gff_output(tmp_path, GALNBKIG_cut):
    os.chdir(tmp_path)
    PF03373 = example_filepath("PF03373.hmm")
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    output = GALNBKIG_cut["gff"]
    codon = GALNBKIG_cut["codon.fasta"]
    amino = GALNBKIG_cut["amino.fasta"]
    r = invoke(
        cli,
        [
            "pscan",
            str(PF03373),
            str(fasta),
            "--output",
            "output.gff",
            "--ocodon",
            "codon.fasta",
            "--oamino",
            "amino.fasta",
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
    assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")


def test_cli_pscan_window0(tmp_path, large_rna):
    os.chdir(tmp_path)
    PF03373 = example_filepath("PF03373.hmm")
    invoke = CliRunner().invoke
    fasta = large_rna["fasta"]
    output = large_rna["output0"]
    codon = large_rna["codon0"]
    amino = large_rna["amino0"]
    r = invoke(
        cli,
        [
            "pscan",
            str(PF03373),
            str(fasta),
            "--output",
            "output.gff",
            "--ocodon",
            "codon.fasta",
            "--oamino",
            "amino.fasta",
            "--window",
            "0",
            "--no-e-value",
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
    assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")


def test_cli_pscan_window48(tmp_path, large_rna):
    os.chdir(tmp_path)
    PF03373 = example_filepath("PF03373.hmm")
    invoke = CliRunner().invoke
    fasta = large_rna["fasta"]
    output = large_rna["output48"]
    codon = large_rna["codon48"]
    amino = large_rna["amino48"]
    r = invoke(
        cli,
        [
            "pscan",
            str(PF03373),
            str(fasta),
            "--output",
            "output.gff",
            "--ocodon",
            "codon.fasta",
            "--oamino",
            "amino.fasta",
            "--window",
            "48",
            "--no-e-value",
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
    assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")


def test_cli_pscan_large_dataset_window():
    profile = example_filepath("PF00113.hmm")
    target = example_filepath("A0ALD9_dna_huge.fasta")
    output = example_filepath("PF00113_A0ALD9_dna_huge_output1776.gff")
    codon = example_filepath("PF00113_A0ALD9_dna_huge_codon1776.fasta")
    amino = example_filepath("PF00113_A0ALD9_dna_huge_amino1776.fasta")

    invoke = CliRunner().invoke
    r = invoke(
        cli,
        [
            "pscan",
            str(profile),
            str(target),
            "--output",
            "output.gff",
            "--ocodon",
            "codon.fasta",
            "--oamino",
            "amino.fasta",
            "--window",
            "1776",
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
    assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")


def test_cli_pscan_window0_evalue(tmp_path):
    os.chdir(tmp_path)
    database1 = example_filepath("database1.hmm")
    invoke = CliRunner().invoke
    fasta = example_filepath("PF03373_rna_most_likely.fasta")
    output = example_filepath("PF03373_rna_output.gff")
    oamino = example_filepath("PF03373_rna_oamino.fasta")
    r = invoke(
        cli,
        [
            "pscan",
            str(database1),
            str(fasta),
            "--output",
            "output.gff",
            "--oamino",
            "amino.fasta",
            "--window",
            "0",
            "--e-value",
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(oamino, "amino.fasta", shallow=False), diff(oamino, "amino.fasta")
