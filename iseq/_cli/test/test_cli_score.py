import os
import shutil
from filecmp import cmp

import pytest
from click.testing import CliRunner

from iseq import cli
from iseq.example import example_filepath

from .misc import diff


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
