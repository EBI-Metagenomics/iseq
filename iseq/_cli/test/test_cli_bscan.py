import os
from filecmp import cmp

from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath
from iseq.file import diff


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

    r = invoke(
        cli,
        [
            "pscan",
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

    assert cmp("bscan_oamino.fasta", "pscan_oamino.fasta", shallow=False), diff(
        "bscan_oamino.fasta", "pscan_oamino.fasta"
    )
    assert cmp("bscan_ocodon.fasta", "pscan_ocodon.fasta", shallow=False), diff(
        "bscan_ocodon.fasta", "pscan_ocodon.fasta"
    )
    assert cmp("bscan_output.gff", "pscan_output.gff", shallow=False), diff(
        "bscan_output.gff", "pscan_output.gff"
    )
