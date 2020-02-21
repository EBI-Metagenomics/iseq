import filecmp

from click.testing import CliRunner
from numpy.testing import assert_equal

from iseq import cli


def test_cli_search_nofile_output(tmpdir, PF03373, GALNBKIG_cut):
    tmpdir.chdir()
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    r = invoke(cli, ["search", str(PF03373), str(fasta)])
    assert_equal(r.exit_code, 0)


def test_cli_search_gff_output(tmpdir, PF03373, GALNBKIG_cut):
    tmpdir.chdir()
    invoke = CliRunner().invoke
    fasta = GALNBKIG_cut["fasta"]
    gff = GALNBKIG_cut["gff"]
    codon = GALNBKIG_cut["codon.fasta"]
    amino = GALNBKIG_cut["amino.fasta"]
    r = invoke(
        cli,
        [
            "search",
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
