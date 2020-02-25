import filecmp
import shutil

from click.testing import CliRunner
from numpy.testing import assert_equal

from iseq import cli


def test_cli_scan_nofile_output(tmpdir, PF03373, GALNBKIG_cut):
    tmpdir.chdir()
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    r = invoke(cli, ["scan", str(PF03373), str(fasta)])
    assert_equal(r.exit_code, 0)


def test_cli_scan_gff_output(tmpdir, PF03373, GALNBKIG_cut):
    tmpdir.chdir()
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    gff = GALNBKIG_cut["gff"]
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
    assert_equal(r.exit_code, 0)
    assert_equal(filecmp.cmp(gff, "output.gff", shallow=False), True)
    assert_equal(filecmp.cmp(codon, "codon.fasta", shallow=False), True)
    assert_equal(filecmp.cmp(amino, "amino.fasta", shallow=False), True)


def test_cli_score(tmpdir, database1, amino1, output1, output1_evalue):
    tmpdir.chdir()
    shutil.copyfile(output1, tmpdir / "output.gff")
    invoke = CliRunner().invoke
    r = invoke(cli, ["score", str(database1), str(amino1), "output.gff"])
    assert_equal(r.exit_code, 0)
    assert_equal(filecmp.cmp(output1_evalue, "output.gff", shallow=False), True)
