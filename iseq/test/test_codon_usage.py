from nmm import Codon, DNAAlphabet, IUPACAminoAlphabet, RNAAlphabet

from iseq.codon_table import CodonTable
from iseq.codon_usage import CodonUsage
from iseq.gencode import GeneticCode

EPS = 1e-9


def test_codon_usage_dna_100():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)
    usage = CodonUsage(table, "100")

    codon = Codon.create(b"CCT", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0018050541516245488) < EPS
    codon = Codon.create(b"CCC", base_abc)
    assert abs(usage.codon_prob(codon) - 0.016245487364620937) < EPS
    codon = Codon.create(b"CCA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS
    codon = Codon.create(b"CCG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.03429602888086643) < EPS

    assert abs(usage.amino_acid_prob(b"P") - 0.05324909747292419) < EPS
    assert abs(usage.amino_acid_prob(b"C") - 0.009927797833935019) < EPS
    assert abs(usage.amino_acid_prob(b"M") - 0.017148014440433214) < EPS

    codon = Codon.create(b"TAG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS

    codon = Codon.create(b"TAA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0) < EPS

    codon = Codon.create(b"TGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.002707581227436823) < EPS

    assert abs(usage.stop_prob - 0.0036101083032490976) < EPS
    assert abs(usage.amino_acid_prob(b"W") - 0.009927797833935019) < EPS

    assert abs(usage.start_prob - 0.058664259927797835) < EPS

    p = sum(usage.amino_acid_prob(aa) for aa in table.amino_acids)
    p += usage.stop_prob
    assert abs(1 - p) < EPS


def test_codon_usage_rna_100():
    base_abc = RNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc)
    usage = CodonUsage(table, "100")

    codon = Codon.create(b"CCU", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0018050541516245488) < EPS
    codon = Codon.create(b"CCC", base_abc)
    assert abs(usage.codon_prob(codon) - 0.016245487364620937) < EPS
    codon = Codon.create(b"CCA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS
    codon = Codon.create(b"CCG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.03429602888086643) < EPS

    assert abs(usage.amino_acid_prob(b"P") - 0.05324909747292419) < EPS
    assert abs(usage.amino_acid_prob(b"C") - 0.009927797833935019) < EPS
    assert abs(usage.amino_acid_prob(b"M") - 0.017148014440433214) < EPS

    codon = Codon.create(b"UAG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS

    codon = Codon.create(b"UAA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0) < EPS

    codon = Codon.create(b"UGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.002707581227436823) < EPS

    assert abs(usage.stop_prob - 0.0036101083032490976) < EPS
    assert abs(usage.amino_acid_prob(b"W") - 0.009927797833935019) < EPS

    assert abs(usage.start_prob - 0.058664259927797835) < EPS

    p = sum(usage.amino_acid_prob(aa) for aa in table.amino_acids)
    p += usage.stop_prob
    assert abs(1 - p) < EPS


def test_codon_usage_dna_100_gencode2():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode("Vertebrate Mitochondrial"))
    usage = CodonUsage(table, "100")

    codon = Codon.create(b"CCT", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0018050541516245488) < EPS
    codon = Codon.create(b"CCC", base_abc)
    assert abs(usage.codon_prob(codon) - 0.016245487364620937) < EPS
    codon = Codon.create(b"CCA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS
    codon = Codon.create(b"CCG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.03429602888086643) < EPS

    assert abs(usage.amino_acid_prob(b"P") - 0.05324909747292419) < EPS
    assert abs(usage.amino_acid_prob(b"C") - 0.009927797833935019) < EPS
    assert abs(usage.amino_acid_prob(b"M") - 0.018050541516245487) < EPS

    codon = Codon.create(b"TAG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS

    codon = Codon.create(b"TAA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0) < EPS

    codon = Codon.create(b"TGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.002707581227436823) < EPS

    codon = Codon.create(b"AGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0) < EPS

    codon = Codon.create(b"AGG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.002707581227436823) < EPS

    assert abs(usage.stop_prob - 0.0036101083032490976) < EPS
    assert abs(usage.amino_acid_prob(b"W") - 0.012635379061371842) < EPS

    assert abs(usage.start_prob - 0.07671480144404333) < EPS

    p = sum(usage.amino_acid_prob(aa) for aa in table.amino_acids)
    p += usage.stop_prob
    assert abs(1 - p) < EPS


def test_codon_usage_dna_100_gencode16():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode(id=16))
    usage = CodonUsage(table, "100")

    codon = Codon.create(b"CCT", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0018050541516245488) < EPS
    codon = Codon.create(b"CCC", base_abc)
    assert abs(usage.codon_prob(codon) - 0.016245487364620937) < EPS
    codon = Codon.create(b"CCA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS
    codon = Codon.create(b"CCG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.03429602888086643) < EPS

    assert abs(usage.amino_acid_prob(b"P") - 0.05324909747292419) < EPS
    assert abs(usage.amino_acid_prob(b"C") - 0.009927797833935019) < EPS
    assert abs(usage.amino_acid_prob(b"M") - 0.017148014440433214) < EPS

    codon = Codon.create(b"TAG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009025270758122744) < EPS

    codon = Codon.create(b"TAA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0) < EPS

    codon = Codon.create(b"TGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.002707581227436823) < EPS

    codon = Codon.create(b"AGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0) < EPS

    codon = Codon.create(b"AGG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.002707581227436823) < EPS

    assert abs(usage.stop_prob - 0.002707581227436823) < EPS
    assert abs(usage.amino_acid_prob(b"W") - 0.009927797833935019) < EPS

    assert abs(usage.start_prob - 0.017148014440433214) < EPS

    p = sum(usage.amino_acid_prob(aa) for aa in table.amino_acids)
    p += usage.stop_prob
    assert abs(1 - p) < EPS


def test_codon_usage_dna_homo_sapiens_gencode1():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode(id=1))
    usage = CodonUsage(table, "9606")

    codon = Codon.create(b"CCT", base_abc)
    assert abs(usage.codon_prob(codon) - 0.01754027818499081) < EPS
    codon = Codon.create(b"CCC", base_abc)
    assert abs(usage.codon_prob(codon) - 0.019787725235942962) < EPS
    codon = Codon.create(b"CCA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.016920666769267137) < EPS
    codon = Codon.create(b"CCG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.006924547978778131) < EPS

    assert abs(usage.amino_acid_prob(b"P") - 0.06117321816897904) < EPS
    assert abs(usage.amino_acid_prob(b"C") - 0.023199190843316342) < EPS
    assert abs(usage.amino_acid_prob(b"M") - 0.022035123101626947) < EPS
    assert abs(usage.amino_acid_prob(b"Q") - 0.046575596207835496) < EPS

    codon = Codon.create(b"TAG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0007896448877742195) < EPS

    codon = Codon.create(b"TAA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009907142640376355) < EPS

    codon = Codon.create(b"TGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0015551644015129192) < EPS

    codon = Codon.create(b"AGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.012165533413495484) < EPS

    codon = Codon.create(b"AGG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.011963406553966494) < EPS

    assert abs(usage.stop_prob - 0.003335523553324774) < EPS
    assert abs(usage.amino_acid_prob(b"W") - 0.013171691851737305) < EPS

    assert abs(usage.start_prob - 0.07460160793527573) < EPS

    p = sum(usage.amino_acid_prob(aa) for aa in table.amino_acids)
    p += usage.stop_prob
    assert abs(1 - p) < EPS


def test_codon_usage_dna_homo_sapiens_gencode15():
    base_abc = DNAAlphabet()
    amino_abc = IUPACAminoAlphabet()

    table = CodonTable(base_abc, amino_abc, GeneticCode(id=15))
    usage = CodonUsage(table, "9606")

    codon = Codon.create(b"CCT", base_abc)
    assert abs(usage.codon_prob(codon) - 0.01754027818499081) < EPS
    codon = Codon.create(b"CCC", base_abc)
    assert abs(usage.codon_prob(codon) - 0.019787725235942962) < EPS
    codon = Codon.create(b"CCA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.016920666769267137) < EPS
    codon = Codon.create(b"CCG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.006924547978778131) < EPS

    assert abs(usage.amino_acid_prob(b"P") - 0.06117321816897904) < EPS
    assert abs(usage.amino_acid_prob(b"C") - 0.023199190843316342) < EPS
    assert abs(usage.amino_acid_prob(b"M") - 0.022035123101626947) < EPS
    assert abs(usage.amino_acid_prob(b"Q") - 0.04736524109560972) < EPS

    codon = Codon.create(b"TAG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0007896448877742195) < EPS

    codon = Codon.create(b"TAA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0009907142640376355) < EPS

    codon = Codon.create(b"TGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.0015551644015129192) < EPS

    codon = Codon.create(b"AGA", base_abc)
    assert abs(usage.codon_prob(codon) - 0.012165533413495484) < EPS

    codon = Codon.create(b"AGG", base_abc)
    assert abs(usage.codon_prob(codon) - 0.011963406553966494) < EPS

    assert abs(usage.stop_prob - 0.0025458786655505545) < EPS
    assert abs(usage.amino_acid_prob(b"W") - 0.013171691851737305) < EPS

    assert abs(usage.start_prob - 0.022035123101626947) < EPS

    p = sum(usage.amino_acid_prob(aa) for aa in table.amino_acids)
    p += usage.stop_prob
    assert abs(1 - p) < EPS
