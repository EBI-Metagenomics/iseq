import os

from assertpy import assert_that, contents_of
from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath


def test_cli_pscan2_pfam24(tmp_path):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke
    profile = example_filepath("Pfam-A_24.hmm")
    fasta = example_filepath("AE014075.1_subset_nucl.fasta")
    oamino = example_filepath("AE014075.1_subset_oamino.fasta")
    ocodon = example_filepath("AE014075.1_subset_ocodon.fasta")
    output = example_filepath("AE014075.1_subset_output.gff")
    r = invoke(
        cli, ["pscan2", str(profile), str(fasta), "--max-e-value", "1e-10", "--quiet"]
    )
    assert r.exit_code == 0, r.output

    assert_that(contents_of("oamino.fasta")).is_equal_to(contents_of(oamino))
    assert_that(contents_of("ocodon.fasta")).is_equal_to(contents_of(ocodon))
    assert_that(contents_of("output.gff")).is_equal_to(contents_of(output))
