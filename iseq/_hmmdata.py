from math import log
from typing import List, Mapping

from hmmer_reader import HMMERModel
from imm import lprob_zero
from nmm import CanonicalAminoAlphabet

from ._alphabet import infer_alphabet
from ._model import Transitions
from ._typing import HMMERAlphabet

__all__ = ["HMMData"]


class HMMData:
    def __init__(self, hmmer: HMMERModel):
        self._original_symbols: str = hmmer.alphabet
        alphabet = infer_alphabet(self._original_symbols.encode())

        if alphabet is None:
            raise ValueError("Could not infer alphabet from HMMER model.")
        self._alphabet = alphabet

        if isinstance(self._alphabet, CanonicalAminoAlphabet):
            self._null_lprobs = _null_amino_lprobs(self._original_symbols)
        else:
            k = alphabet.length
            self._null_lprobs = [log(1 / k)] * k

        self._match_lprobs = [self._sort(hmmer.match(i)) for i in range(1, hmmer.M + 1)]
        self._insert_lprobs = [
            self._sort(hmmer.insert(i)) for i in range(1, hmmer.M + 1)
        ]
        self._model_length = hmmer.M

        self._transitions: List[Transitions] = []
        for m in range(0, hmmer.M + 1):
            t = Transitions(**hmmer.trans(m))
            self._transitions.append(t)

    @property
    def transitions(self) -> List[Transitions]:
        return self._transitions

    @property
    def model_length(self) -> int:
        return self._model_length

    @property
    def alphabet(self) -> HMMERAlphabet:
        return self._alphabet

    @property
    def null_lprobs(self) -> List[float]:
        return self._null_lprobs

    def match_lprobs(self, m: int) -> List[float]:
        return self._match_lprobs[m - 1]

    def insert_lprobs(self, m: int) -> List[float]:
        return self._insert_lprobs[m - 1]

    def _sort(self, lprobs: Mapping[str, float]) -> List[float]:
        symbols = self._alphabet.symbols.decode()
        return [lprobs.get(sym, lprob_zero()) for sym in symbols]


def _null_amino_lprobs(symbols: str):
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
    return [lprobs.get(sym, lprob_zero()) for sym in list(symbols)]
