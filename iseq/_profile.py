from abc import ABC, abstractmethod
from math import log, exp
from typing import Generic, Dict

from nmm.alphabet import Alphabet
from nmm.prob import lprob_zero
from nmm.sequence import Sequence

from ._model import AltModel, NullModel, SpecialTransitions
from ._result import SearchResults
from ._typing import TAlphabet, TState

__all__ = ["Profile"]


class Profile(Generic[TAlphabet, TState], ABC):
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._multiple_hits: bool = True
        self._special_transitions = SpecialTransitions()

    @property
    def alphabet(self):
        return self._alphabet

    @property
    def null_model(self) -> NullModel:
        raise NotImplementedError()

    @property
    def alt_model(self) -> AltModel:
        raise NotImplementedError()

    @property
    def multiple_hits(self) -> bool:
        return self._multiple_hits

    @multiple_hits.setter
    def multiple_hits(self, multiple_hits: bool):
        self._multiple_hits = multiple_hits

    @abstractmethod
    def search(
        self, sequence: Sequence, window_length: int = 0
    ) -> SearchResults[TAlphabet, TState]:
        del sequence
        del window_length
        raise NotImplementedError()

    def _set_special_transitions(self, length: int):
        L = length
        if L == 0:
            return

        if self._multiple_hits:
            l1q = lq = -log(2)
        else:
            lq = lprob_zero()
            l1q = log(1.0)

        q = exp(lq)
        lp = log(L) - log(L + 2 + q / (1 - q))
        l1p = log(2 + q / (1 - q)) - log(L + 2 + q / (1 - q))
        lr = log(L) - log(L + 1)

        t = self._special_transitions

        t.NN = t.CC = t.JJ = lp
        t.NB = t.CT = t.JB = l1p
        t.RR = lr
        t.EJ = lq
        t.EC = l1q
