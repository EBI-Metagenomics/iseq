import pytest
import filecmp
import shutil

from click.testing import CliRunner
from numpy.testing import assert_equal

from iseq import cli
from iseq._misc import download, tempdir, brotli_decompress


@pytest.mark.slow
def test_cli_large_dataset_window():
    with tempdir():
        base = "https://rest.s3for.me/iseq"
        target = download(f"{base}/A0ALD9_dna_huge.fasta.br")
        profile = download(f"{base}/PF00113.hmm.br")
        output = download(f"{base}/PF00113_A0ALD9_dna_huge_output1776.gff")
        codon = download(f"{base}/PF00113_A0ALD9_dna_huge_codon1776.fasta")
        amino = download(f"{base}/PF00113_A0ALD9_dna_huge_amino1776.fasta")
    pass
