import os

from assertpy import assert_that, contents_of
from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath


def test_cli_bscan_GALNBKIG_pfam10(tmp_path):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    profile = example_filepath("Pfam-A.33.1_10.hmm")
    fasta = example_filepath("GALNBKIG_00914_ont_01_plus_strand.fasta")

    r = invoke(
        cli,
        [
            "pscan",
            str(profile),
            str(fasta),
            "--output",
            "pscan_output.gff",
            "--ocodon",
            "pscan_ocodon.fasta",
            "--oamino",
            "pscan_oamino.fasta",
            "--quiet",
        ],
    )
    assert r.exit_code == 0, r.output

    r = invoke(cli, ["press", str(profile)])
    assert r.exit_code == 0, r.output

    r = invoke(
        cli,
        [
            "bscan",
            str(profile),
            str(fasta),
            "--output",
            "bscan_output.gff",
            "--ocodon",
            "bscan_ocodon.fasta",
            "--oamino",
            "bscan_oamino.fasta",
            "--quiet",
        ],
    )
    assert r.exit_code == 0, r.output

    assert_that(contents_of("bscan_oamino.fasta")).is_equal_to(
        contents_of("pscan_oamino.fasta")
    )
    assert_that(contents_of("bscan_ocodon.fasta")).is_equal_to(
        contents_of("pscan_ocodon.fasta")
    )
    assert_that(contents_of("bscan_output.gff")).is_equal_to(
        contents_of("pscan_output.gff")
    )


def test_cli_bscan_GALNBKIG_pfam10_ncpus2(tmp_path):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    profile = example_filepath("Pfam-A.33.1_10.hmm")
    fasta = example_filepath("GALNBKIG_00914_ont_01_plus_strand.fasta")

    r = invoke(
        cli,
        [
            "pscan",
            str(profile),
            str(fasta),
            "--output",
            "pscan_output.gff",
            "--ocodon",
            "pscan_ocodon.fasta",
            "--oamino",
            "pscan_oamino.fasta",
            "--quiet",
        ],
    )
    assert r.exit_code == 0, r.output

    r = invoke(cli, ["press", str(profile)])
    assert r.exit_code == 0, r.output

    r = invoke(
        cli,
        [
            "bscan",
            str(profile),
            str(fasta),
            "--output",
            "bscan_output.gff",
            "--ocodon",
            "bscan_ocodon.fasta",
            "--oamino",
            "bscan_oamino.fasta",
            "--quiet",
            "--ncpus",
            2,
        ],
    )
    assert r.exit_code == 0, r.output

    assert_that(contents_of("bscan_oamino.fasta")).is_equal_to(
        contents_of("pscan_oamino.fasta")
    )
    assert_that(contents_of("bscan_ocodon.fasta")).is_equal_to(
        contents_of("pscan_ocodon.fasta")
    )
    assert_that(contents_of("bscan_output.gff")).is_equal_to(
        contents_of("pscan_output.gff")
    )
