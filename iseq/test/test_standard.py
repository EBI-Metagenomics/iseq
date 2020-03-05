from numpy.testing import assert_allclose, assert_equal

from hmmer_reader import open_hmmer
from fasta_reader import open_fasta
from iseq.standard import create_standard_profile, create_hmmer3_profile
from nmm.sequence import Sequence


def test_standard_profile_unihit_homologous_1(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    most_likely_seq = Sequence.create(b"PGKEDNNK", alphabet)
    r = hmmer.search(most_likely_seq).results[0]

    assert_allclose(r.loglikelihood, 11.867796719423442)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(most_likely_seq))

    hmmer.multiple_hits = False
    r = hmmer.search(most_likely_seq).results[0]
    assert_allclose(r.loglikelihood, 11.94063404337571)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(most_likely_seq))


def test_hmmer3_profile_problematic1(problematic1):
    with open_hmmer(problematic1["hmm"]) as reader:
        prof = create_hmmer3_profile(reader.read_profile())

    with open_fasta(problematic1["fasta"]) as reader:
        item = reader.read_items()[0]

    sequence = Sequence.create(item.sequence.encode(), prof.alphabet)
    prof._set_special_transitions(len(sequence))
    prof._alt_model.update_special_transitions(hmmer3=True)
    results = prof.alt_model.viterbi(sequence)

    assert len(results) == 1
    assert_allclose(results[0].loglikelihood, -2.103729125681)


def test_standard_profile_unihit_homologous_2(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PGKENNK", alphabet)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 3.299501501364073)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(seq))
    assert_equal(str(frag)[:31], "('P', '<M1,1>'),('G', '<M2,1>')")


def test_standard_profile_unihit_homologous_3(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PGKEPNNK", alphabet)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 6.883636719423446)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(seq))


def test_standard_profile_nonhomo_and_homologous(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"KKKPGKEDNNK", alphabet)
    assert_equal(hmmer.multiple_hits, True)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 10.707618955640605)
    frags = r.fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"KKK")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")

    hmmer.multiple_hits = False
    assert_equal(hmmer.multiple_hits, False)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 10.96037578075283)
    frags = r.fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"KKK")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")


def test_standard_profile_multihit_homologous1(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PPPPGKEDNNKDDDPGKEDNNKEEEE", alphabet)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 20.329227532144742)
    frags = r.fragments
    assert_equal(len(frags), 5)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"PPP")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")
    assert_equal(frags[3].homologous, True)
    assert_equal(bytes(frags[3].sequence), b"PGKEDNNK")
    assert_equal(frags[4].homologous, False)
    assert_equal(bytes(frags[4].sequence), b"EEEE")

    items = list(iter(frags[0]))

    assert_equal(bytes(items[0].sequence), b"")
    assert_equal(str(items[0].step), "<S,0>")
    assert_equal(bytes(items[1].sequence), b"P")
    assert_equal(str(items[1].step), "<N,1>")
    assert_equal(bytes(items[4].sequence), b"")
    assert_equal(str(items[4].step), "<B,0>")

    items = list(frags[1])

    assert_equal(bytes(items[0].sequence), b"P")
    assert_equal(str(items[0].step), "<M1,1>")
    assert_equal(bytes(items[1].sequence), b"G")
    assert_equal(str(items[1].step), "<M2,1>")
    assert_equal(bytes(items[2].sequence), b"K")
    assert_equal(str(items[2].step), "<M3,1>")
    assert_equal(bytes(items[3].sequence), b"E")
    assert_equal(str(items[3].step), "<M4,1>")
    assert_equal(bytes(items[4].sequence), b"D")
    assert_equal(str(items[4].step), "<M5,1>")
    assert_equal(bytes(items[5].sequence), b"N")
    assert_equal(str(items[5].step), "<M6,1>")
    assert_equal(bytes(items[6].sequence), b"N")
    assert_equal(str(items[6].step), "<M7,1>")
    assert_equal(bytes(items[7].sequence), b"K")
    assert_equal(str(items[7].step), "<M8,1>")

    hmmer.multiple_hits = False
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 8.666478660222928)
    frags = r.fragments
    assert_equal(len(frags), 3)
    assert_equal(frags[0].homologous, False)
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")
    assert_equal(frags[2].homologous, False)


def test_standard_profile_window(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PPPPGKEDNNKDDDPGKEDNNKEEEE", alphabet)
    r = hmmer.search(seq, 15)

    assert_allclose(r.results[0].loglikelihood, 8.868483903457452)
    frags = r.results[0].fragments
    assert_equal(len(frags), 3)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"PPP")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")
    assert_equal(frags[2].homologous, False)
    assert_equal(bytes(frags[2].sequence), b"DDDP")

    assert_allclose(r.results[1].loglikelihood, 12.123181172648234)
    frags = r.results[1].fragments
    assert_equal(len(frags), 3)
    assert_equal(frags[0].homologous, True)
    assert_equal(bytes(frags[0].sequence), b"DNNK")
    assert_equal(frags[1].homologous, False)
    assert_equal(bytes(frags[1].sequence), b"DDD")
    assert_equal(frags[2].homologous, True)
    assert_equal(bytes(frags[2].sequence), b"PGKEDNNK")

    assert_allclose(r.results[2].loglikelihood, 9.082860795403896)
    frags = r.results[2].fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, True)
    assert_equal(bytes(frags[0].sequence), b"PGKEDNNK")
    assert_equal(frags[1].homologous, False)
    assert_equal(bytes(frags[1].sequence), b"EEEE")
