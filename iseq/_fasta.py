from typing import Optional, Union

from fasta_reader import FASTAParser
from nmm import CanonicalAminoAlphabet, DNAAlphabet, RNAAlphabet

from ._alphabet import infer_alphabet

Alphabets = Union[DNAAlphabet, RNAAlphabet, CanonicalAminoAlphabet]


def infer_fasta_alphabet(parser: FASTAParser) -> Optional[Alphabets]:

    for item in parser:
        alphabet = infer_alphabet(item.sequence.encode())
        if alphabet is not None:
            return alphabet

    return None
