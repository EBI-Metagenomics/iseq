from filecmp import cmp

import pytest
from click.testing import CliRunner

from iseq import cli
from iseq._misc import brotli_decompress, diff, download, tempdir


@pytest.mark.slow
def test_cli_large_dataset_window():
    with tempdir():
        base = "https://rest.s3for.me/iseq"
        profile_zip = download(f"{base}/PF00113.hmm.br")
        target_zip = download(f"{base}/A0ALD9_dna_huge.fasta.br")
        output = download(f"{base}/PF00113_A0ALD9_dna_huge_output1776.gff")
        codon = download(f"{base}/PF00113_A0ALD9_dna_huge_codon1776.fasta")
        amino = download(f"{base}/PF00113_A0ALD9_dna_huge_amino1776.fasta")

        profile = brotli_decompress(profile_zip)
        target = brotli_decompress(target_zip)

        invoke = CliRunner().invoke
        r = invoke(
            cli,
            [
                "scan",
                str(profile),
                str(target),
                "--output",
                "output.gff",
                "--ocodon",
                "codon.fasta",
                "--oamino",
                "amino.fasta",
                "--window",
                "1776",
            ],
        )
        assert r.exit_code == 0, r.output
        assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
        assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
        assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")
