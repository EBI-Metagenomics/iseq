import os
import shutil
from filecmp import cmp

import pytest
from click.testing import CliRunner

from iseq import cli, file_example
from .misc import diff


@pytest.fixture
def GALNBKIG_cut():
    return {
        "fasta": file_example("GALNBKIG_cut.fasta"),
        "gff": file_example("PF03373_GALNBKIG_cut.gff"),
        "amino.fasta": file_example("PF03373_GALNBKIG_cut.amino.fasta"),
        "codon.fasta": file_example("PF03373_GALNBKIG_cut.codon.fasta"),
    }


@pytest.fixture
def large_rna():
    return {
        "fasta": file_example("large_rna_seq.fasta"),
        "amino0": file_example("large_rna_seq_amino0.fasta"),
        "codon0": file_example("large_rna_seq_codon0.fasta"),
        "output0": file_example("large_rna_seq_output0.gff"),
        "amino48": file_example("large_rna_seq_amino48.fasta"),
        "codon48": file_example("large_rna_seq_codon48.fasta"),
        "output48": file_example("large_rna_seq_output48.gff"),
    }


def test_cli_scan_nofile_output(tmp_path, GALNBKIG_cut):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    PF03373 = file_example("PF03373.hmm")
    r = invoke(cli, ["scan", str(PF03373), str(fasta)])
    assert r.exit_code == 0, r.output


def test_cli_scan_gff_output(tmp_path, GALNBKIG_cut):
    os.chdir(tmp_path)
    PF03373 = file_example("PF03373.hmm")
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


def test_cli_scan_window0(tmp_path, large_rna):
    os.chdir(tmp_path)
    PF03373 = file_example("PF03373.hmm")
    invoke = CliRunner().invoke
    fasta = large_rna["fasta"]
    output = large_rna["output0"]
    codon = large_rna["codon0"]
    amino = large_rna["amino0"]
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
            "--window",
            "0",
        ],
    )
    assert r.exit_code == 0, r.output
    assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
    assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
    assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")


def test_cli_scan_window48(tmp_path, large_rna):
    os.chdir(tmp_path)
    PF03373 = file_example("PF03373.hmm")
    invoke = CliRunner().invoke
    fasta = large_rna["fasta"]
    output = large_rna["output48"]
    codon = large_rna["codon48"]
    amino = large_rna["amino48"]
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
    output1 = file_example("output1.gff")
    shutil.copyfile(output1, tmp_path / "output.gff")

    database1 = file_example("database1.hmm")
    amino1 = file_example("amino1.fasta")
    output1_evalue = file_example("output1_evalue.gff")

    invoke = CliRunner().invoke
    r = invoke(cli, ["score", str(database1), str(amino1), "output.gff"])
    assert r.exit_code == 0, r.output
    assert cmp(output1_evalue, "output.gff", shallow=False), diff(
        output1_evalue, "output.gff"
    )
