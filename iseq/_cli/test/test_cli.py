import shutil
from filecmp import cmp

from click.testing import CliRunner

from iseq import cli
from iseq._misc import diff


def test_cli_scan_nofile_output(tmpdir, PF03373, GALNBKIG_cut):
    tmpdir.chdir()
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    r = invoke(cli, ["scan", str(PF03373), str(fasta)])
    assert r.exit_code == 0, r.output


def test_cli_scan_gff_output(tmpdir, PF03373, GALNBKIG_cut):
    tmpdir.chdir()
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


def test_cli_scan_window0(tmpdir, PF03373, large_rna):
    tmpdir.chdir()
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


def test_cli_scan_window48(tmpdir, PF03373, large_rna):
    tmpdir.chdir()
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


def test_cli_score(tmpdir, database1, amino1, output1, output1_evalue):
    tmpdir.chdir()
    shutil.copyfile(output1, tmpdir / "output.gff")
    invoke = CliRunner().invoke
    r = invoke(cli, ["score", str(database1), str(amino1), "output.gff"])
    assert r.exit_code == 0, r.output
    assert cmp(output1_evalue, "output.gff", shallow=False), diff(
        output1_evalue, "output.gff"
    )
