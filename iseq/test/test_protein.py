from hmmer_reader import open_hmmer
from imm import Sequence
from imm.testing import assert_allclose
from nmm import IUPACAminoAlphabet, RNAAlphabet

from iseq.codon_table import CodonTable
from iseq.example import example_filepath
from iseq.hmmer_model import HMMERModel
from iseq.protein import create_profile


def test_protein_profile_frame1():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet())

    rna_abc = hmmer.alphabet
    most_likely_rna_seq = b"CCU GGU AAA GAA GAU AAU AAC AAA".replace(b" ", b"")
    most_likely_seq = Sequence.create(most_likely_rna_seq, rna_abc)
    r = hmmer.search(most_likely_seq).results[0]

    assert_allclose(r.loglikelihood, 125.83363182422178)
    frags = r.fragments
    assert len(frags) == 1
    frag = frags[0]
    assert frag.homologous
    assert bytes(frag.sequence) == bytes(most_likely_seq)
    desired = "('CCU', '<M1,3>'),('GGU', '<M2,3>'),('AAA', '<M3,3>')"
    assert str(frag)[:53] == desired


def test_protein_profile_frame2():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet(), epsilon=0.1)

    rna_abc = hmmer.alphabet
    rna_seq = b"AAA AAA AAA CCU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")
    seq = Sequence.create(rna_seq, rna_abc)

    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 168.23071232889802)
    frags = r.fragments
    assert len(frags) == 2
    assert not frags[0].homologous
    assert bytes(frags[0].sequence) == b"AAAAAAAAA"
    assert frags[1].homologous
    assert bytes(frags[1].sequence) == b"CCUGGUAAAGAAGAUAAUAACAAA"


def test_protein_profile_frame3():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet(), epsilon=0.0)

    rna_abc = hmmer.alphabet
    rna_seq = b"CCU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")
    seq = Sequence.create(rna_seq, rna_abc)

    r = hmmer.search(seq).results[0]
    frags = r.fragments
    assert len(frags) == 1
    assert frags[0].homologous
    assert bytes(frags[0].sequence) == b"CCUGGUAAAGAAGAUAAUAACAAA"


def test_protein_profile_frame4():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet(), epsilon=0.0)

    rna_abc = hmmer.alphabet
    rna_seq = b"CCUU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")
    seq = Sequence.create(rna_seq, rna_abc)

    r = hmmer.search(seq).results[0]
    frags = r.fragments
    assert len(frags) == 0


def test_protein_profile_frame5():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet(), epsilon=0.00001)

    rna_abc = hmmer.alphabet
    rna_seq = b"CCUU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")
    seq = Sequence.create(rna_seq, rna_abc)

    r = hmmer.search(seq).results[0]
    frags = r.fragments
    assert len(frags) == 1
    assert frags[0].homologous
    assert bytes(frags[0].sequence) == b"CCUUGGUAAAGAAGAUAAUAACAAA"


def test_protein_profile_frame6():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet(), epsilon=0.00001)

    rna_abc = hmmer.alphabet
    rna_seq = b"CCUU GGU AAA GAA GAU AAU AAC AAA GAA GAA CCU GGU AAA GAA GAU AAU AAC AAA GAA GAA GA"
    rna_seq = rna_seq.replace(b" ", b"")
    seq = Sequence.create(rna_seq, rna_abc)

    r = hmmer.search(seq).results[0]
    frags = r.fragments
    assert len(frags) == 4
    assert frags[0].homologous
    assert bytes(frags[0].sequence) == b"CCUUGGUAAAGAAGAUAAUAACAAA"
    assert not frags[1].homologous
    assert bytes(frags[1].sequence) == b"GAAGAA"
    assert frags[2].homologous
    assert bytes(frags[2].sequence) == b"CCUGGUAAAGAAGAUAAUAACAAA"
    assert not frags[3].homologous
    assert bytes(frags[3].sequence) == b"GAAGAAGA"

    hmmer.multiple_hits = False
    r = hmmer.search(seq).results[0]
    frags = r.fragments
    assert_allclose(r.loglikelihood, 1445.0314253859958)
    assert len(frags) == 3
    assert not frags[0].homologous
    assert bytes(frags[0].sequence) == b"CCUUGGUAAAGAAGAUAAUAACAAAGAAGAA"
    assert frags[1].homologous
    assert bytes(frags[1].sequence) == b"CCUGGUAAAGAAGAUAAUAACAAA"
    assert not frags[2].homologous
    assert bytes(frags[2].sequence) == b"GAAGAAGA"


def test_protein_profile_codons():
    filepath = example_filepath("PF03373.hmm")
    with open_hmmer(filepath) as reader:
        plain_model = reader.read_model()
        hmmer = create_profile(HMMERModel(plain_model), RNAAlphabet(), epsilon=0.1)

    rna_abc = hmmer.alphabet
    amino_abc = IUPACAminoAlphabet()
    gcode = CodonTable(rna_abc, amino_abc)

    rna_seq = b"AAGA AAA AAA CCU GGU AAA GAA GAU AAU AAC AAA G"
    rna_seq = rna_seq.replace(b" ", b"")
    seq = Sequence.create(rna_seq, rna_abc)

    r = hmmer.search(seq).results[0]
    assert_allclose(r.loglikelihood, 175.35113397356454)
    frags = r.fragments
    cfrags = [f.decode() for f in frags]
    aafrags = [f.decode(gcode) for f in cfrags]

    assert len(frags) == 2
    assert len(cfrags) == 2
    assert len(aafrags) == 2

    assert not frags[0].homologous
    assert not cfrags[0].homologous
    assert not aafrags[0].homologous

    assert bytes(frags[0].sequence) == b"AAGAAAAAAA"
    assert bytes(cfrags[0].sequence) == b"AAGAAAAAA"
    assert bytes(aafrags[0].sequence) == b"KKK"

    items = list(iter(frags[0]))
    citems = list(iter(cfrags[0]))
    aaitems = list(iter(aafrags[0]))

    assert bytes(items[0].sequence) == b""
    assert str(items[0].step) == "<S,0>"
    assert bytes(citems[0].sequence) == b""
    assert str(citems[0].step) == "<S,0>"
    assert bytes(aaitems[0].sequence) == b""
    assert str(aaitems[0].step) == "<S,0>"

    assert bytes(items[1].sequence) == b"AAG"
    assert str(items[1].step) == "<N,3>"
    assert bytes(citems[1].sequence) == b"AAG"
    assert str(citems[1].step) == "<N,3>"
    assert bytes(aaitems[1].sequence) == b"K"
    assert str(aaitems[1].step) == "<N,1>"

    assert bytes(items[2].sequence) == b"AAAA"
    assert str(items[2].step) == "<N,4>"
    assert bytes(citems[2].sequence) == b"AAA"
    assert str(citems[2].step) == "<N,3>"
    assert bytes(aaitems[2].sequence) == b"K"
    assert str(aaitems[2].step) == "<N,1>"

    assert bytes(items[3].sequence) == b"AAA"
    assert str(items[3].step) == "<N,3>"
    assert bytes(citems[3].sequence) == b"AAA"
    assert str(citems[3].step) == "<N,3>"
    assert bytes(aaitems[3].sequence) == b"K"
    assert str(aaitems[3].step) == "<N,1>"

    assert bytes(items[4].sequence) == b""
    assert str(items[4].step) == "<B,0>"
    assert bytes(citems[4].sequence) == b""
    assert str(citems[4].step) == "<B,0>"
    assert bytes(aaitems[4].sequence) == b""
    assert str(aaitems[4].step) == "<B,0>"

    assert frags[1].homologous
    assert cfrags[1].homologous
    assert aafrags[1].homologous
    assert bytes(frags[1].sequence) == b"CCUGGUAAAGAAGAUAAUAACAAAG"
    assert bytes(cfrags[1].sequence) == b"CCUGGUAAAGAAGAUAAUAACAAG"

    assert bytes(aafrags[1].sequence) == b"PGKEDNNK"

    items = list(iter(frags[1]))
    citems = list(iter(cfrags[1]))
    aaitems = list(iter(aafrags[1]))

    assert bytes(items[0].sequence) == b"CCU"
    assert str(items[0].step) == "<M1,3>"
    assert bytes(citems[0].sequence) == b"CCU"
    assert str(citems[0].step) == "<M1,3>"
    assert bytes(aaitems[0].sequence) == b"P"
    assert str(aaitems[0].step) == "<M1,1>"

    assert bytes(items[7].sequence) == b"AAAG"
    assert str(items[7].step) == "<M8,4>"
    assert bytes(citems[7].sequence) == b"AAG"
    assert str(citems[7].step) == "<M8,3>"
    assert bytes(aaitems[7].sequence) == b"K"
    assert str(aaitems[7].step) == "<M8,1>"
