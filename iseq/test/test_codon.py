# from numpy.testing import assert_allclose, assert_equal

# from hmmer_reader import open_hmmer
# from iseq.codon import create_profile
# from nmm import GeneticCode
# from nmm.alphabet import IUPACAminoAlphabet, RNAAlphabet
# from nmm.sequence import Sequence


# def test_codon_profile1(PF03373):
#     with open_hmmer(PF03373) as reader:
#         hmmer = create_profile(reader.read_profile(), RNAAlphabet())

#     rna_abc = hmmer.alphabet
#     most_likely_rna_seq = b"CCU GGU AAA GAA GAU AAU AAC AAA".replace(b" ", b"")
#     most_likely_seq = Sequence.create(most_likely_rna_seq, rna_abc)
#     r = hmmer.search(most_likely_seq).results[0]

#     assert_allclose(r.loglikelihood, 125.83363182422178)
#     frags = r.fragments
#     assert_equal(len(frags), 1)
#     frag = frags[0]
#     assert_equal(frag.homologous, True)
#     assert_equal(bytes(frag.sequence), bytes(most_likely_seq))
#     desired = "('CCU', '<M1,3>'),('GGU', '<M2,3>'),('AAA', '<M3,3>')"
#     assert_equal(str(frag)[:53], desired)
