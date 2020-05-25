from typing import Optional, Union

from fasta_reader import FASTAParser
from nmm import CanonicalAminoAlphabet, DNAAlphabet, RNAAlphabet

Alphabets = Union[DNAAlphabet, RNAAlphabet, CanonicalAminoAlphabet]

__all__ = ["Alphabets", "infer_alphabet"]


def infer_alphabet(sequence: bytes) -> Optional[Alphabets]:
    """
    Infer alphabet from a sequence of symbols.

    Parameters
    ----------
    sequence
        Sequence of symbols.
    """
    dna = DNAAlphabet()
    rna = RNAAlphabet()
    amino = CanonicalAminoAlphabet()

    abc = set(sequence)

    if len(abc - set(dna.symbols)) == 0:
        return dna

    if len(abc - set(rna.symbols)) == 0:
        return rna

    if len(abc - set(amino.symbols)) == 0:
        return amino

    return None


def infer_fasta_alphabet(parser: FASTAParser) -> Optional[Alphabets]:
    """
    Infer alphabet from fasta file.

    Parameters
    ----------
    parser
        FASTA parser.
    """

    for item in parser:
        alphabet = infer_alphabet(item.sequence.encode())
        if alphabet is not None:
            return alphabet

    return None
