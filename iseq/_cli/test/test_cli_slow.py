from filecmp import cmp

import pytest
from click.testing import CliRunner

from iseq import cli
from iseq._misc import diff, download, tmp_cwd
from iseq._file import brotli_decompress


# @pytest.mark.slow
# def test_cli_large_dataset_window():
#     with tmp_cwd():
#         profile = get_filepath("PF00113.hmm")
#         target = get_filepath("A0ALD9_dna_huge.fasta")
#         output = get_filepath("PF00113_A0ALD9_dna_huge_output1776.gff")
#         codon = get_filepath("PF00113_A0ALD9_dna_huge_codon1776.fasta")
#         amino = get_filepath("PF00113_A0ALD9_dna_huge_amino1776.fasta")

#         invoke = CliRunner().invoke
#         r = invoke(
#             cli,
#             [
#                 "scan",
#                 str(profile),
#                 str(target),
#                 "--output",
#                 "output.gff",
#                 "--ocodon",
#                 "codon.fasta",
#                 "--oamino",
#                 "amino.fasta",
#                 "--window",
#                 "1776",
#             ],
#         )
#         assert r.exit_code == 0, r.output
#         assert cmp(output, "output.gff", shallow=False), diff(output, "output.gff")
#         assert cmp(codon, "codon.fasta", shallow=False), diff(codon, "codon.fasta")
#         assert cmp(amino, "amino.fasta", shallow=False), diff(amino, "amino.fasta")
