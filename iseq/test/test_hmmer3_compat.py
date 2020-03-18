from filecmp import cmp

import pytest
from click.testing import CliRunner

from iseq import cli
from iseq._misc import brotli_decompress, diff, download, tmp_cwd


# @pytest.mark.slow
def test_hmmer3_viterbi_scores_compat():
    with tmp_cwd():
        base = "https://rest.s3for.me/iseq"
        profiles_zip = download(f"{base}/Pfam-A.hmm.br")
        target_zip = download(f"{base}/A0ALD9.fasta.br")
        viterbi_scores_zip = download(f"{base}/Pfam-A_hmmer3.3_viterbi_scores.txt.br")

        # output = download(f"{base}/PF00113_A0ALD9_dna_huge_output1776.gff")
        # codon = download(f"{base}/PF00113_A0ALD9_dna_huge_codon1776.fasta")
        # amino = download(f"{base}/PF00113_A0ALD9_dna_huge_amino1776.fasta")

        profiles = brotli_decompress(profiles_zip)
        target = brotli_decompress(target_zip)
        viterbi_scores = brotli_decompress(viterbi_scores_zip)
        pass

        # invoke = CliRunner().invoke
        # r = invoke(
        #     cli,
        #     [
        #         "scan",
        #         str(profile),
        #         str(target),
        #         "--output",
        #         "output.gff",
        #         "--ocodon",
        #         "codon.fasta",
        #         "--oamino",
        #         "amino.fasta",
        #         "--window",
        #         "1776",
        #     ],
        # )
        # assert r.exit_code == 0, r.output
        # assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
        # assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
        # assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")
