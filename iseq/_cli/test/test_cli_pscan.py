import os
import shutil
from filecmp import cmp

import pytest
from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath

from .misc import diff


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
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
    assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")


@pytest.mark.skipif(
    shutil.which("hmmsearch") is None, reason="requires HMMER3 software"
)
def test_cli_score(tmp_path):
    os.chdir(tmp_path)
    output1 = example_filepath("output1.gff")
    shutil.copyfile(output1, tmp_path / "output.gff")

    database1 = example_filepath("database1.hmm")
    amino1 = example_filepath("amino1.fasta")
    output1_evalue = example_filepath("output1_evalue.gff")

    invoke = CliRunner().invoke
    r = invoke(cli, ["score", str(database1), str(amino1), "output.gff"])
    assert r.exit_code == 0, r.output
    assert cmp(output1_evalue, "output.gff", shallow=False), diff(
        output1_evalue, "output.gff"
    )
