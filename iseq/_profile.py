from abc import ABC, abstractmethod
from math import log
from typing import Generic, Tuple, TypeVar, List

from nmm.alphabet import Alphabet
from nmm.prob import lprob_zero
from nmm.sequence import Sequence

from ._model import AltModel, NullModel
from ._result import SearchResult
from ._type import MutablePath, TState

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

__all__ = ["Profile"]


class Profile(Generic[TAlphabet, TState], ABC):
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._multiple_hits: bool = True

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
        self, sequence: Sequence, window: int = 0
    ) -> List[SearchResult[TAlphabet, TState]]:
        del sequence
        del window
        raise NotImplementedError()

    def _search(
        self, sequence: Sequence, window: int = 0
    ) -> List[Tuple[float, MutablePath[TState]]]:

        self._set_target_length(len(sequence))
        results = self.alt_model.viterbi(sequence, window)
        score_path_list: List[Tuple[float, MutablePath[TState]]] = []

        for result in results:
            score0 = self.null_model.likelihood(result.sequence)
            score1 = result.loglikelihood
            score_path_list.append((score1 - score0, result.path))

        return score_path_list

    def _set_fragment_length(self):
        if self.alt_model.length == 0:
            return

        B = self.alt_model.special_node.B
        E = self.alt_model.special_node.E

        # Uniform local alignment fragment length distribution
        t = self.alt_model.special_transitions
        t.BM = log(2) - log(self.alt_model.length) - log(self.alt_model.length + 1)
        t.ME = 0.0
        for node in self.alt_model.core_nodes():
            self.alt_model.set_transition(B, node.M, t.BM)
            self.alt_model.set_transition(node.M, E, t.ME)

        for node in self.alt_model.core_nodes()[1:]:
            self.alt_model.set_transition(node.D, E, 0.0)

    def _set_target_length(self, length: int):
        from math import exp

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

        t = self.alt_model.special_transitions

        t.NN = t.CC = t.JJ = lp
        t.NB = t.CT = t.JB = l1p
        t.RR = lr
        t.EJ = lq
        t.EC = l1q

        node = self.alt_model.special_node

        self.alt_model.set_transition(node.S, node.B, t.NB)
        self.alt_model.set_transition(node.S, node.N, t.NN)
        self.alt_model.set_transition(node.N, node.N, t.NN)
        self.alt_model.set_transition(node.N, node.B, t.NB)

        self.alt_model.set_transition(node.E, node.T, t.EC + t.CT)
        self.alt_model.set_transition(node.E, node.C, t.EC + t.CC)
        self.alt_model.set_transition(node.C, node.C, t.CC)
        self.alt_model.set_transition(node.C, node.T, t.CT)

        self.alt_model.set_transition(node.E, node.B, t.EJ + t.JB)
        self.alt_model.set_transition(node.E, node.J, t.EJ + t.JJ)
        self.alt_model.set_transition(node.J, node.J, t.JJ)
        self.alt_model.set_transition(node.J, node.B, t.JB)

        self.null_model.set_transition(t.RR)
