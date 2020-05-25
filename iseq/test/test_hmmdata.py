from math import log

from numpy.testing import assert_allclose

from hmmer_reader import open_hmmer
from iseq import HMMData, file_example
from nmm import CanonicalAminoAlphabet, DNAAlphabet


def test_hmmdata_amino():
    filepath = file_example("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmm_data = HMMData(reader.read_model())

    assert isinstance(hmm_data.alphabet, CanonicalAminoAlphabet)
    assert_allclose(hmm_data.null_lprobs[3], -2.7056061901315998)
    assert_allclose(hmm_data.match_lprobs(6)[3], -3.36585)
    assert_allclose(hmm_data.insert_lprobs(1)[0], -2.68618)


def test_hmmdata_dna():
    filepath = file_example("ecori.hmm")
    with open_hmmer(filepath) as reader:
        hmm_data = HMMData(reader.read_model())

    assert isinstance(hmm_data.alphabet, DNAAlphabet)
    assert_allclose(hmm_data.null_lprobs[3], log(1 / 4))
    assert_allclose(hmm_data.match_lprobs(2)[2], -2.42319)
    assert_allclose(hmm_data.insert_lprobs(1)[0], -1.38629)
