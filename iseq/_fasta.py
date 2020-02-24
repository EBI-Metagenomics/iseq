from typing import Union

from fasta_reader import FASTAParser
from nmm.alphabet import DNAAlphabet, RNAAlphabet

from ._alphabet import infer_alphabet


def infer_fasta_alphabet(parser: FASTAParser) -> Union[DNAAlphabet, RNAAlphabet, None]:
    for item in parser:
        alphabet = infer_alphabet(item.sequence.encode())
        if alphabet is not None:
            return alphabet

    return None
