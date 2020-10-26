from nmm import Codon, DNAAlphabet, IUPACAminoAlphabet, RNAAlphabet

from iseq.codon_table import CodonTable
from iseq.gencode import GeneticCode


def test_codon_table_dna():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCT", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    assert len(table.codons(b"W")) == 1
    assert Codon.create(b"TGG", base_abc) in table.codons(b"W")

    assert table.amino_acid(Codon.create(b"ATG", base_abc)) == b"M"
    assert len(table.amino_acids) == 20
    assert b"R" in table.amino_acids

    assert len(table.stop_codons) == 3
    assert Codon.create(b"TAA", base_abc) in table.stop_codons
    assert Codon.create(b"TAG", base_abc) in table.stop_codons
    assert Codon.create(b"TGA", base_abc) in table.stop_codons

    assert len(set(table.codons())) == 64


def test_codon_table_dna_id33():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    gencode = GeneticCode(id=33)
    table = CodonTable(base_abc, amino_abc, gencode)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCT", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    assert len(table.codons(b"W")) == 2
    assert Codon.create(b"TGG", base_abc) in table.codons(b"W")
    assert Codon.create(b"TGA", base_abc) in table.codons(b"W")

    assert len(table.codons(b"T")) == 4
    assert Codon.create(b"ACT", base_abc) in table.codons(b"T")
    assert Codon.create(b"ACC", base_abc) in table.codons(b"T")
    assert Codon.create(b"ACA", base_abc) in table.codons(b"T")
    assert Codon.create(b"ACG", base_abc) in table.codons(b"T")

    assert table.amino_acid(Codon.create(b"ATG", base_abc)) == b"M"
    assert len(table.amino_acids) == 20
    assert b"R" in table.amino_acids

    assert len(table.start_codons) == 4
    assert Codon.create(b"TTG", base_abc) in table.start_codons
    assert Codon.create(b"CTG", base_abc) in table.start_codons
    assert Codon.create(b"ATG", base_abc) in table.start_codons
    assert Codon.create(b"GTG", base_abc) in table.start_codons

    assert len(table.stop_codons) == 1
    assert Codon.create(b"TAG", base_abc) in table.stop_codons

    assert len(set(table.codons())) == 64


def test_codon_table_rna():
    base_abc = RNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)

    assert len(table.codons(b"P")) == 4

    assert Codon.create(b"CCU", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in table.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in table.codons(b"P")

    assert table.amino_acid(Codon.create(b"AUG", base_abc)) == b"M"
    assert len(table.amino_acids) == 20
    assert b"R" in table.amino_acids

    assert len(table.start_codons) == 3
    assert Codon.create(b"UUG", base_abc) in table.start_codons
    assert Codon.create(b"CUG", base_abc) in table.start_codons
    assert Codon.create(b"AUG", base_abc) in table.start_codons

    assert len(table.stop_codons) == 3
    assert Codon.create(b"UAA", base_abc) in table.stop_codons
    assert Codon.create(b"UAG", base_abc) in table.stop_codons
    assert Codon.create(b"UGA", base_abc) in table.stop_codons

    assert len(set(table.codons())) == 64
