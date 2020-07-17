from imm.testing import assert_allclose
from iseq.example import example_filepath
from iseq.profmark import ProfMark


def test_profmark():
    hmmer = example_filepath("Pfam-A_24.hmm")

    domtblout = example_filepath("AE014075.1_domtblout.txt")
    cds_amino = example_filepath("AE014075.1_cds_amino.fasta")

    cds_nucl = example_filepath("AE014075.1_cds_nucl.fasta")
    output = example_filepath("AE014075.1_output.gff")

    pm = ProfMark(hmmer, cds_amino, domtblout, cds_nucl, output)

    tpr = [
        0.0,
        0.07142857142857142,
        0.14285714285714285,
        0.21428571428571427,
        0.2857142857142857,
        0.35714285714285715,
        0.42857142857142855,
        0.5,
        0.5714285714285714,
        0.5714285714285714,
        0.6428571428571429,
        0.7142857142857143,
        0.7142857142857143,
        0.7857142857142857,
        0.8571428571428571,
        0.9285714285714286,
        1.0,
    ]
    breakpoint()
    assert_allclose(pm.confusion_matrix.tpr[:len(tpr)], tpr)
    fpr = [0.0, 0.0, 0.004424778761061954, 0.11504424778761058, 0.995575221238938, 1.0]
    assert_allclose(pm.confusion_matrix.fpr[[0, 5, 10, 40, -2, -1]], fpr)
