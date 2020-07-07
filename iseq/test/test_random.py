from nmm import DNAAlphabet, IUPACAminoAlphabet, RNAAlphabet

from iseq.codon_table import CodonTable
from iseq.codon_usage import CodonUsage
from iseq.gencode import GeneticCode
from iseq.random import RandomState


def test_random_codon_dna():

    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode(id=1))
    usage = CodonUsage(table, "9606")

    random = RandomState(0)

    codons = random.codon(usage)
    assert len(codons) == 1
    assert codons[0].symbols == b"GAC"

    codons = random.codon(usage, 3)
    assert len(codons) == 3
    assert codons[0].symbols == b"GGC"
    assert codons[1].symbols == b"GAG"
    assert codons[2].symbols == b"GAC"

    codons = random.codon(usage, 3)
    assert len(codons) == 3
    assert codons[0].symbols == b"CGG"
    assert codons[1].symbols == b"GCC"
    assert codons[2].symbols == b"CTA"

    random = RandomState(0)

    amino_acids = random.amino_acid(usage)
    assert len(amino_acids) == 1
    assert amino_acids[0] == b"D"

    amino_acids = random.amino_acid(usage, 3)
    assert len(amino_acids) == 3
    assert amino_acids[0] == b"G"
    assert amino_acids[1] == b"E"
    assert amino_acids[2] == b"D"

    amino_acids = random.amino_acid(usage, 3)
    assert len(amino_acids) == 3
    assert amino_acids[0] == b"R"
    assert amino_acids[1] == b"A"
    assert amino_acids[2] == b"L"

    amino_acids = random.amino_acid(usage, 500)
    assert amino_acids[442] == b"*"

    amino_acids = random.amino_acid(usage, 1000, include_stop=False)
    assert b"*" not in amino_acids
    assert len(amino_acids) == 1000


def test_random_codon_rna():

    base_abc = RNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode(id=1))
    usage = CodonUsage(table, "9606")

    random = RandomState(0)

    codons = random.codon(usage)
    assert len(codons) == 1
    assert codons[0].symbols == b"GAC"

    codons = random.codon(usage, 3)
    assert len(codons) == 3
    assert codons[0].symbols == b"GGC"
    assert codons[1].symbols == b"GAG"
    assert codons[2].symbols == b"GAC"

    codons = random.codon(usage, 3)
    assert len(codons) == 3
    assert codons[0].symbols == b"CGG"
    assert codons[1].symbols == b"GCC"
    assert codons[2].symbols == b"CUA"

    random = RandomState(0)

    amino_acids = random.amino_acid(usage)
    assert len(amino_acids) == 1
    assert amino_acids[0] == b"D"

    amino_acids = random.amino_acid(usage, 3)
    assert len(amino_acids) == 3
    assert amino_acids[0] == b"G"
    assert amino_acids[1] == b"E"
    assert amino_acids[2] == b"D"

    amino_acids = random.amino_acid(usage, 3)
    assert len(amino_acids) == 3
    assert amino_acids[0] == b"R"
    assert amino_acids[1] == b"A"
    assert amino_acids[2] == b"L"

    amino_acids = random.amino_acid(usage, 500)
    assert amino_acids[442] == b"*"

    amino_acids = random.amino_acid(usage, 1000, include_stop=False)
    assert b"*" not in amino_acids
    assert len(amino_acids) == 1000


def test_random_codon_dna_id():

    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode(id=3))
    usage = CodonUsage(table, "8084")

    random = RandomState(0)

    codons = random.codon(usage)
    assert len(codons) == 1
    assert codons[0].symbols == b"GAC"

    codons = random.codon(usage, 3)
    assert len(codons) == 3
    assert codons[0].symbols == b"GGC"
    assert codons[1].symbols == b"GAT"
    assert codons[2].symbols == b"GAC"

    codons = random.codon(usage, 3)
    assert len(codons) == 3
    assert codons[0].symbols == b"CGG"
    assert codons[1].symbols == b"GCC"
    assert codons[2].symbols == b"CTA"

    random = RandomState(0)

    amino_acids = random.amino_acid(usage)
    assert len(amino_acids) == 1
    assert amino_acids[0] == b"D"

    amino_acids = random.amino_acid(usage, 3)
    assert len(amino_acids) == 3
    assert amino_acids[0] == b"G"
    assert amino_acids[1] == b"D"
    assert amino_acids[2] == b"D"

    amino_acids = random.amino_acid(usage, 3)
    assert len(amino_acids) == 3
    assert amino_acids[0] == b"R"
    assert amino_acids[1] == b"A"
    assert amino_acids[2] == b"T"

    amino_acids = random.amino_acid(usage, 5000)
    assert b"*" in amino_acids

    amino_acids = random.amino_acid(usage, 5000, include_stop=False)
    assert b"*" not in amino_acids
