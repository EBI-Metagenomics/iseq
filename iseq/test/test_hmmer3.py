import os

from fasta_reader import read_fasta
from hmmer_reader import open_hmmer
from imm import Sequence
from imm.testing import assert_allclose
from numpy.testing import assert_equal

from iseq.example import example_filepath
from iseq.hmmer3 import create_profile
from iseq.hmmer_model import HMMERModel
from iseq.model import EntryDistr


def test_hmmer3_profile_unihit_homologous_1():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmmdata = HMMERModel(reader.read_model())

    hmmer = create_profile(hmmdata, entry_distr=EntryDistr.UNIFORM)

    alphabet = hmmer.alphabet
    most_likely_seq = Sequence.create(b"PGKEDNNK", alphabet)
    r = hmmer.search(most_likely_seq).results[0]

    assert_allclose(r.loglikelihood, 13.39978964458627)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(most_likely_seq))

    hmmer.multiple_hits = False
    r = hmmer.search(most_likely_seq).results[0]
    assert_allclose(r.loglikelihood, 13.472626968538535)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(most_likely_seq))


def test_hmmer3_profile_problematic1():
    with open_hmmer(example_filepath("problematic1.hmm")) as reader:
        hmmdata = HMMERModel(reader.read_model())

    prof = create_profile(hmmdata, True)

    with read_fasta(example_filepath("problematic1.fasta")) as reader:
        item = reader.read_items()[0]

    sequence = Sequence.create(item.sequence.encode(), prof.alphabet)
    r = prof.search(sequence)
    assert len(r.results) == 1
    assert_allclose(r.results[0].alt_viterbi_score, -5.103729125680717)


def test_hmmer3_profile_small_viterbi_score(tmp_path):
    os.chdir(tmp_path)
    profile = example_filepath("PF15449.6.hmm")
    fasta = example_filepath("A0ALD9.fasta")

    with open_hmmer(profile) as reader:
        hmmdata = HMMERModel(reader.read_model())

    prof = create_profile(hmmdata, True)

    with read_fasta(fasta) as reader:
        item = reader.read_items()[0]

    sequence = Sequence.create(item.sequence.encode(), prof.alphabet)
    r = prof.search(sequence)
    assert len(r.results) == 1
    assert_allclose(r.results[0].alt_viterbi_score, -18.424065160005625)


def test_hmmer3_profile_large_viterbi_score(tmp_path):
    os.chdir(tmp_path)
    profile = example_filepath("PF07476.11.hmm")
    fasta = example_filepath("A0ALD9.fasta")

    with open_hmmer(profile) as reader:
        hmmdata = HMMERModel(reader.read_model())

    prof = create_profile(hmmdata, True)

    with read_fasta(fasta) as reader:
        item = reader.read_items()[0]

    sequence = Sequence.create(item.sequence.encode(), prof.alphabet)
    r = prof.search(sequence)
    assert len(r.results) == 1
    assert_allclose(r.results[0].alt_viterbi_score, 3.480341268180834)


def test_hmmer3_profile_unihit_homologous_2():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmmdata = HMMERModel(reader.read_model())

    hmmer = create_profile(hmmdata, entry_distr=EntryDistr.UNIFORM)

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PGKENNK", alphabet)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 4.679028352468922)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(seq))
    assert_equal(str(frag)[:31], "('P', '<M1,1>'),('G', '<M2,1>')")


def test_hmmer3_profile_unihit_homologous_3():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmmdata = HMMERModel(reader.read_model())

    hmmer = create_profile(hmmdata, entry_distr=EntryDistr.UNIFORM)

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PGKEPNNK", alphabet)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 8.55629338255035)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(bytes(frag.sequence), bytes(seq))


def test_hmmer3_profile_nonhomo_and_homologous():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmmdata = HMMERModel(reader.read_model())

    hmmer = create_profile(hmmdata, entry_distr=EntryDistr.UNIFORM)

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"KKKPGKEDNNK", alphabet)
    assert_equal(hmmer.multiple_hits, True)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 12.23961188080344)
    frags = r.fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"KKK")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")

    hmmer.multiple_hits = False
    assert_equal(hmmer.multiple_hits, False)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 12.492368705915666)
    frags = r.fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"KKK")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")


def test_hmmer3_profile_multihit_homologous1():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmmdata = HMMERModel(reader.read_model())

    hmmer = create_profile(hmmdata, entry_distr=EntryDistr.OCCUPANCY)

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PPPPGKEDNNKDDDPGKEDNNKEEEE", alphabet)
    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 23.396127957291384)
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
    assert_allclose(r.loglikelihood, 10.19992887279622)
    frags = r.fragments
    assert_equal(len(frags), 3)
    assert_equal(frags[0].homologous, False)
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")
    assert_equal(frags[2].homologous, False)


def test_hmmer3_profile_window():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        hmmdata = HMMERModel(reader.read_model())

    hmmer = create_profile(hmmdata, entry_distr=EntryDistr.UNIFORM, window_length=15)

    alphabet = hmmer.alphabet
    seq = Sequence.create(b"PPPPGKEDNNKDDDPGKEDNNKEEEE", alphabet)
    r = hmmer.search(seq)

    assert_allclose(r.results[0].loglikelihood, 10.400476828620292)
    frags = r.results[0].fragments
    assert_equal(len(frags), 3)
    assert_equal(frags[0].homologous, False)
    assert_equal(bytes(frags[0].sequence), b"PPP")
    assert_equal(frags[1].homologous, True)
    assert_equal(bytes(frags[1].sequence), b"PGKEDNNK")
    assert_equal(frags[2].homologous, False)
    assert_equal(bytes(frags[2].sequence), b"DDDP")

    assert_allclose(r.results[1].loglikelihood, 14.513128891939548)
    frags = r.results[1].fragments
    assert_equal(len(frags), 3)
    assert_equal(frags[0].homologous, True)
    assert_equal(bytes(frags[0].sequence), b"DNNK")
    assert_equal(frags[1].homologous, False)
    assert_equal(bytes(frags[1].sequence), b"DDD")
    assert_equal(frags[2].homologous, True)
    assert_equal(bytes(frags[2].sequence), b"PGKEDNNK")

    assert_allclose(r.results[2].loglikelihood, 10.614853720566728)
    frags = r.results[2].fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, True)
    assert_equal(bytes(frags[0].sequence), b"PGKEDNNK")
    assert_equal(frags[1].homologous, False)
    assert_equal(bytes(frags[1].sequence), b"EEEE")
