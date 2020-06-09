import os
from filecmp import cmp

from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath

from .misc import diff


def test_cli_scan_nofile_output(tmp_path, GALNBKIG_cut):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    PF03373 = example_filepath("PF03373.hmm")
    r = invoke(cli, ["scan", str(PF03373), str(fasta)])
    assert r.exit_code == 0, r.output


def test_cli_scan_gff_output(tmp_path, GALNBKIG_cut):
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
            "scan",
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
