import pytest


@pytest.fixture
def GALNBKIG_cut():
    from iseq.example import example_filepath

    return {
        "fasta": example_filepath("GALNBKIG_cut.fasta"),
        "gff": example_filepath("PF03373_GALNBKIG_cut.gff"),
        "amino.fasta": example_filepath("PF03373_GALNBKIG_cut.amino.fasta"),
        "codon.fasta": example_filepath("PF03373_GALNBKIG_cut.codon.fasta"),
    }


@pytest.fixture
def large_rna():
    from iseq.example import example_filepath

    return {
        "fasta": example_filepath("large_rna_seq.fasta"),
        "amino0": example_filepath("large_rna_seq_amino0.fasta"),
        "codon0": example_filepath("large_rna_seq_codon0.fasta"),
        "output0": example_filepath("large_rna_seq_output0.gff"),
        "amino48": example_filepath("large_rna_seq_amino48.fasta"),
        "codon48": example_filepath("large_rna_seq_codon48.fasta"),
        "output48": example_filepath("large_rna_seq_output48.gff"),
    }
