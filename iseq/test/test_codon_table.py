import pytest

from nmm import CanonicalAminoAlphabet, Codon, DNAAlphabet, RNAAlphabet
from iseq.codon_table import CodonTable
from iseq.gencode import get_genetic_code


def test_genetic_code_dna():
    base_abc = DNAAlphabet()
    amino_abc = CanonicalAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCT", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    assert len(table.codons(b"W")) == 1
    assert Codon.create(b"TGG", base_abc) in table.codons(b"W")

    assert table.amino_acid(Codon.create(b"ATG", base_abc)) == b"M"
    assert len(table.amino_acids()) == 20
    assert b"R" in table.amino_acids()


def test_genetic_code_dna_id33():
    base_abc = DNAAlphabet()
    amino_abc = CanonicalAminoAlphabet()

    gencode = get_genetic_code(id=33)
    table = CodonTable(base_abc, amino_abc, gencode)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCT", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    assert len(table.codons(b"W")) == 2
    assert Codon.create(b"TGG", base_abc) in table.codons(b"W")
    assert Codon.create(b"TGA", base_abc) in table.codons(b"W")

    assert table.amino_acid(Codon.create(b"ATG", base_abc)) == b"M"
    assert len(table.amino_acids()) == 20
    assert b"R" in table.amino_acids()


def test_genetic_code_rna():
    base_abc = RNAAlphabet()
    amino_abc = CanonicalAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCU", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    assert table.amino_acid(Codon.create(b"AUG", base_abc)) == b"M"
    assert len(table.amino_acids()) == 20
    assert b"R" in table.amino_acids()
    # from nmm._gencode import get_codon_probs

    # with pytest.raises(ValueError):
    #     get_codon_probs("homo sapiens")

    # table = get_codon_probs("Homo sapiens")


def test_genetic_code_rna_prob():
    base_abc = RNAAlphabet()
    amino_abc = CanonicalAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCU", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    codon = Codon.create(b"CCG", base_abc)
    assert abs(1 / 4 - table.codons_prob(b"P")[codon]) < 1e-7

    assert table.amino_acid(Codon.create(b"AUG", base_abc)) == b"M"
    assert len(table.amino_acids()) == 20
    assert b"R" in table.amino_acids()
