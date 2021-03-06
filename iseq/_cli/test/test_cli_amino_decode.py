import os
from typing import Dict

import pytest
from assertpy import assert_that, contents_of
from click.testing import CliRunner

from iseq import cli

_amino = """>Homoserine_dh-consensus
PIISTLKESLTGDRITRIEGILNGTLNYILTEMEEEGASFSEALKEAQELGYAEADPTDD
VEGLDAARKLAILARLAFGLEVELEDVEVEGIEKLTAEDIEEAKEEGKVLKLVASAVEAR
VKPELVPKSHPLASVKGSDNAVAVETERVGELVVQGPGAGAEPTASAVLADLL
>AA_kinase-consensus
KRVVVKLGGSSLTDKEEASLRRLAEQIAALKESGNKLVVVHGGGSFTDGLLALKSGLSSG
ELAAGLRSTLEEAGEVATRDALASLGERLVAALLAAGLPAVGLSAAALDATEAGRDEGSD
GNVESVDAEAIEELLEAGVVPVLTGFIGLDEEGELGRGSSDTIAALLAEALGADKLIILT
DVDGVYDADPKKVPDARLLPEISVDEAEESASELATGGMKVKHPAALAAARRGGIPVVIT
N
>23ISL-consensus
QGLDNANRSLVRATKAESSDIRKEVTNGIAKGLKLDSLETAAESKNCSSAQKGGSLAWAT
NSQPQPLRESKLEPLEDSPRKALKTPVLQKTSSTITLQAVKVQPEPRAPVSGALSPSGEE
RKRPAASAPATLPTRQSGLGSQEVVSKVATRKIPMESQREST
"""

_nucl_transl: Dict[int, str] = {}

_nucl_transl[
    1
] = """>Homoserine_dh-consensus
CCTATCATTTCGACGCTCAAGGAGTCGCTGACAGGTGACCGTATTACTCGAATCGAAGGG
ATATTAAACGGCACCCTGAATTACATTCTCACTGAGATGGAGGAAGAGGGGGCTTCATTC
TCTGAGGCGCTGAAGGAGGCACAGGAATTGGGCTACGCGGAAGCGGATCCTACGGACGAT
GTGGAAGGGCTAGATGCTGCTAGAAAGCTGGCAATTCTAGCCAGATTGGCATTTGGGTTA
GAGGTCGAGTTGGAGGACGTAGAGGTGGAAGGAATTGAAAAGCTGACTGCCGAAGATATT
GAAGAAGCGAAGGAAGAGGGTAAAGTTTTAAAACTAGTGGCAAGCGCCGTCGAAGCCAGG
GTCAAGCCTGAGCTGGTACCTAAGTCACATCCATTAGCCTCGGTAAAAGGCTCTGACAAC
GCCGTGGCTGTAGAAACGGAACGGGTAGGCGAACTCGTAGTGCAGGGACCAGGGGCTGGC
GCAGAGCCAACCGCATCCGCTGTACTCGCTGACCTTCTC
>AA_kinase-consensus
AAACGTGTAGTTGTAAAGCTTGGGGGTAGTTCTCTGACAGATAAGGAAGAGGCATCACTC
AGGCGTTTAGCTGAGCAGATTGCAGCATTAAAAGAGAGTGGCAATAAACTAGTGGTCGTG
CATGGAGGCGGCAGCTTCACTGATGGTCTGCTGGCATTGAAAAGTGGCCTGAGCTCGGGC
GAATTAGCTGCGGGGTTGAGGAGCACGTTAGAAGAGGCCGGAGAAGTAGCGACGAGGGAC
GCCCTAGCTAGCTTAGGGGAACGGCTTGTTGCAGCGCTGCTGGCGGCGGGTCTCCCTGCT
GTAGGACTCAGCGCCGCTGCGTTAGATGCGACGGAGGCGGGCCGGGATGAAGGCAGCGAC
GGGAACGTCGAGTCCGTGGACGCAGAAGCAATTGAGGAGTTGCTTGAGGCCGGGGTGGTC
CCCGTCCTAACAGGATTTATCGGCTTAGACGAAGAAGGGGAACTGGGAAGGGGATCTTCT
GACACCATCGCTGCGTTACTCGCTGAAGCTTTAGGCGCGGACAAACTCATAATACTGACC
GACGTAGACGGCGTTTACGATGCCGACCCTAAAAAGGTCCCAGACGCGAGGCTCTTGCCA
GAGATAAGTGTGGACGAGGCCGAGGAAAGCGCCTCCGAATTAGCGACCGGTGGGATGAAG
GTCAAACATCCAGCGGCTCTTGCTGCAGCTAGACGGGGGGGTATTCCGGTCGTGATAACG
AAT
>23ISL-consensus
CAGGGTCTGGATAACGCTAATCGTTCGCTAGTTCGCGCTACAAAAGCAGAAAGTTCAGAT
ATACGGAAAGAGGTGACTAACGGCATCGCTAAAGGGCTGAAGCTAGACAGTCTGGAAACA
GCTGCAGAGTCGAAGAACTGCTCAAGCGCACAGAAAGGCGGATCGCTAGCTTGGGCAACC
AACTCCCAACCACAGCCTCTCCGTGAAAGTAAGCTTGAGCCATTGGAAGACTCCCCACGT
AAGGCTTTAAAAACACCTGTGTTGCAAAAGACATCCAGTACCATAACTTTACAAGCAGTC
AAGGTTCAACCTGAACCCCGCGCTCCCGTCTCCGGGGCGCTGTCCCCGAGCGGGGAGGAA
CGCAAGCGCCCAGCTGCGTCTGCTCCCGCTACCTTACCGACACGACAGAGTGGTCTAGGT
TCTCAGGAAGTCGTTTCGAAGGTGGCGACTCGCAAAATTCCAATGGAGTCACAACGCGAG
TCGACT
"""

_nucl_transl[
    9
] = """>Homoserine_dh-consensus
CCTATCATTTCGACGCTCAAGGAGTCCCTCACCGGAGATCGGATAACTCGTATTGAAGGC
ATATTAAACGGCACCCTGAATTACATTCTCACTGAGATGGAGGAAGAGGGGGCTTCATTC
TCTGAGGCGCTGAAGGAGGCGCAAGAGTTAGGCTACGCCGAGGCTGACCCAACTGACGAC
GTAGAGGGACTAGATGCTGCTCGTAAGCTGGCCATATTAGCTCGCCTAGCCTTTGGATTA
GAGGTCGAGTTGGAGGACGTAGAGGTGGAAGGAATTGAAAAGCTCACCGCTGAGGATATT
GAGGAAGCAAAGGAAGAGGGGAAGGTTCTCAAGCTAGTTGCTTCTGCAGTGGAAGCGCGA
GTGAAGCCGGAGCTGGTCCCTAAGTCCCACCCCCTCGCTAGGGTCAAGGGATCTGACAAA
GCTGTAGCTGTCGAGACAGAACGCGTTGGGGAGTTGGTGGTTCAAGGACCGGGAGCGGGG
GCGGAACCCACAGCATCGGCAGTGCTCGCAGATCTCCTA
>AA_kinase-consensus
AAGCGCGTAGTGGTAAAGTTGGGAGGCAGTTCACTCACTGACAAGGAAGAGGCTTCTTTA
CGACGGCTTGCGGAACAAATTGCCGCATTAAAGGAGTCAGGGAAAAAGCTTGTGGTCGTT
CATGGTGGGGGGTCTTTTACAGATGGTCTCCTAGCCTTAAAGAGTGGGCTGTCGAGTGGA
GAGTTGGCCGCGGGTTTACGTAGCACCTTGGAAGAAGCCGGAGAGGTCGCGACCCGTGAT
GCTCTCGCGTCCCTGGGGGAGCGGTTAGTTGCGGCCCTTTTAGCAGCGGGGCTGCCCGCC
GTGGGACTATCTGCCGCTGCGCTTGACGCGACAGAAGCAGGGCGAGATGAGGGCAGGGAC
GGTAATGTTGAAAGAGTGGACGCCGAAGCGATTGAAGAGCTCTTGGAGGCCGGGGTAGTT
CCCGTCCTCACGGGGTTCATCGGCCTCGATGAAGAAGGTGAGTTGGGACGGGGCAGGAGG
GACACCATCGCTGCACTTTTAGCCGAGGCTCTAGGTGCGGATAAGCTGATAATCTTAACA
GATGTCGACGGCGTTTACGATGCGGATCCTAAGAAGGTTCCTGACGCGCGGCTACTCCCA
GAAATCTCCGTCGATGAGGCCGAAGAGTCAGCCAGCGAACTAGCCACCGGAGGGATGAAG
GTGAAGCACCCGGCGGCGTTGGCAGCGGCACGTCGGGGGGGCATCCCCGTAGTCATCACC
AAA
>23ISL-consensus
CAAGGGTTGGATAACGCCAAACGTAGATTAGTACGTGCAACTAAGGCTGAGTCGTCGGAT
ATTCGGAAGGAGGTGACAAATGGGATTGCGAAGGGGCTGAAGCTTGACTCTCTAGAGACT
GCGGCTGAATCCAAGAATTGTAGAAGAGCTCAAAAGGGAGGTAGACTCGCATGGGCGACT
AACAGCCAGCCTCAACCGCTGCGGGAAAGGAAGCTAGAGCCACTTGAAGATTCCCCGCGG
AAGGCGTTGAAGACACCCGTACTCCAAAAGACCTCATCGACGATTACTCTTCAGGCGGTC
AAGGTACAACCCGAACCACGGGCTCCAGTTTCGGGAGCCCTTTCCCCTAGCGGCGAAGAA
CGGAAGCGTCCCGCTGCGTCTGCTCCAGCTACGTTGCCTACGCGACAGAGGGGTTTGGGA
AGTCAAGAAGTAGTCAGGAAGGTTGCTACTCGTAAGATCCCGATGGAAAGACAGCGTGAG
AGCACC
"""

_nucl_transl[
    11
] = """>Homoserine_dh-consensus
CCTATCATTTCGACGCTCAAGGAGTCGCTGACAGGTGACCGTATTACTCGAATCGAAGGG
ATATTAAACGGCACCCTGAATTACATTCTCACTGAGATGGAGGAAGAGGGGGCTTCATTC
TCTGAGGCGCTGAAGGAGGCACAGGAATTGGGCTACGCGGAAGCGGATCCTACGGACGAT
GTGGAAGGGCTAGATGCTGCTAGAAAGCTGGCAATTCTAGCCAGATTGGCATTTGGGTTA
GAGGTCGAGTTGGAGGACGTAGAGGTGGAAGGAATTGAAAAGCTGACTGCCGAAGATATT
GAAGAAGCGAAGGAAGAGGGTAAAGTTTTAAAACTAGTGGCAAGCGCCGTCGAAGCCAGG
GTCAAGCCTGAGCTGGTACCTAAGTCACATCCATTAGCCTCGGTAAAAGGCTCTGACAAC
GCCGTGGCTGTAGAAACGGAACGGGTAGGCGAACTCGTAGTGCAGGGACCAGGGGCTGGC
GCAGAGCCAACCGCATCCGCTGTACTCGCTGACCTTCTC
>AA_kinase-consensus
AAACGTGTAGTTGTAAAGCTTGGGGGTAGTTCTCTGACAGATAAGGAAGAGGCATCACTC
AGGCGTTTAGCTGAGCAGATTGCAGCATTAAAAGAGAGTGGCAATAAACTAGTGGTCGTG
CATGGAGGCGGCAGCTTCACTGATGGTCTGCTGGCATTGAAAAGTGGCCTGAGCTCGGGC
GAATTAGCTGCGGGGTTGAGGAGCACGTTAGAAGAGGCCGGAGAAGTAGCGACGAGGGAC
GCCCTAGCTAGCTTAGGGGAACGGCTTGTTGCAGCGCTGCTGGCGGCGGGTCTCCCTGCT
GTAGGACTCAGCGCCGCTGCGTTAGATGCGACGGAGGCGGGCCGGGATGAAGGCAGCGAC
GGGAACGTCGAGTCCGTGGACGCAGAAGCAATTGAGGAGTTGCTTGAGGCCGGGGTGGTC
CCCGTCCTAACAGGATTTATCGGCTTAGACGAAGAAGGGGAACTGGGAAGGGGATCTTCT
GACACCATCGCTGCGTTACTCGCTGAAGCTTTAGGCGCGGACAAACTCATAATACTGACC
GACGTAGACGGCGTTTACGATGCCGACCCTAAAAAGGTCCCAGACGCGAGGCTCTTGCCA
GAGATAAGTGTGGACGAGGCCGAGGAAAGCGCCTCCGAATTAGCGACCGGTGGGATGAAG
GTCAAACATCCAGCGGCTCTTGCTGCAGCTAGACGGGGGGGTATTCCGGTCGTGATAACG
AAT
>23ISL-consensus
CAGGGTCTGGATAACGCTAATCGTTCGCTAGTTCGCGCTACAAAAGCAGAAAGTTCAGAT
ATACGGAAAGAGGTGACTAACGGCATCGCTAAAGGGCTGAAGCTAGACAGTCTGGAAACA
GCTGCAGAGTCGAAGAACTGCTCAAGCGCACAGAAAGGCGGATCGCTAGCTTGGGCAACC
AACTCCCAACCACAGCCTCTCCGTGAAAGTAAGCTTGAGCCATTGGAAGACTCCCCACGT
AAGGCTTTAAAAACACCTGTGTTGCAAAAGACATCCAGTACCATAACTTTACAAGCAGTC
AAGGTTCAACCTGAACCCCGCGCTCCCGTCTCCGGGGCGCTGTCCCCGAGCGGGGAGGAA
CGCAAGCGCCCAGCTGCGTCTGCTCCCGCTACCTTACCGACACGACAGAGTGGTCTAGGT
TCTCAGGAAGTCGTTTCGAAGGTGGCGACTCGCAAAATTCCAATGGAGTCACAACGCGAG
TCGACT
"""


@pytest.mark.parametrize("trans_tbl", [1, 9, 11])
def test_cli_amino_decode_transl9(tmp_path, trans_tbl):
    os.chdir(tmp_path)
    invoke = CliRunner().invoke

    with open("amino.fa", "w") as file:
        file.write(_amino)

    with open("desired.fa", "w") as file:
        file.write(_nucl_transl[trans_tbl])

    r = invoke(
        cli,
        [
            "amino-decode",
            "amino.fa",
            "--output",
            "nucl.fa",
            "--transl-table",
            trans_tbl,
        ],
    )
    assert r.exit_code == 0, r.output
    assert_that(contents_of("nucl.fa")).is_equal_to(contents_of("desired.fa"))
