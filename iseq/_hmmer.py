from typing import Optional, Union

from hmmer_reader import HMMERParser
from nmm import CanonicalAminoAlphabet, DNAAlphabet, RNAAlphabet

Alphabets = Union[DNAAlphabet, RNAAlphabet, CanonicalAminoAlphabet]


def infer_hmmer_alphabet(parser: HMMERParser) -> Optional[Alphabets]:

    for prof in parser:
        alph = dict(prof.metadata)["ALPH"]
        if alph == "amino":
            return CanonicalAminoAlphabet()
        if alph == "dna":
            return DNAAlphabet()
        if alph == "rna":
            return RNAAlphabet()

    return None
