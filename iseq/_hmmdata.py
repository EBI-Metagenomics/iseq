from math import log
from nmm.prob import lprob_zero
from typing import Union, List
from hmmer_reader import HMMERModel
from nmm.alphabet import RNAAlphabet, DNAAlphabet, CanonicalAminoAlphabet

from ._alphabet import infer_alphabet

HMMERAlphabet = Union[RNAAlphabet, DNAAlphabet, CanonicalAminoAlphabet]


class HMMData:
    def __init__(self, hmmer: HMMERModel):
        alphabet = infer_alphabet(hmmer.alphabet.encode())
        if alphabet is None:
            raise ValueError("Could not infer alphabet from HMMER model.")
        self._alphabet = alphabet

        if isinstance(self._alphabet, CanonicalAminoAlphabet):
            self._null_lprobs = _null_amino_lprobs(alphabet)
        else:
            k = alphabet.length
            self._null_lprobs = [log(1 / k)] * k

    @property
    def alphabet(self) -> HMMERAlphabet:
        return self._alphabet

    @property
    def null_lprobs(self) -> List[float]:
        return self._null_lprobs


def _null_amino_lprobs(alphabet: CanonicalAminoAlphabet):
    """
    Copy/paste from HMMER3 amino acid frequences infered form Swiss-Prot 50.8,
    (Oct 2006), counting over 85956127 (86.0M) residues.
    """
    lprobs = {
        "A": log(0.0787945),
        "C": log(0.0151600),
        "D": log(0.0535222),
        "E": log(0.0668298),
        "F": log(0.0397062),
        "G": log(0.0695071),
        "H": log(0.0229198),
        "I": log(0.0590092),
        "K": log(0.0594422),
        "L": log(0.0963728),
        "M": log(0.0237718),
        "N": log(0.0414386),
        "P": log(0.0482904),
        "Q": log(0.0395639),
        "R": log(0.0540978),
        "S": log(0.0683364),
        "T": log(0.0540687),
        "V": log(0.0673417),
        "W": log(0.0114135),
        "Y": log(0.0304133),
    }
    abc = list(alphabet.symbols.decode())
    return [lprobs.get(sym, lprob_zero()) for sym in abc]
