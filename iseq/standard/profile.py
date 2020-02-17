from typing import Any, Dict, List, Sequence, Tuple

from nmm import LPROB_ZERO, Alphabet, CSequence, MuteState, NormalState, lprob_normalize

from hmmer_reader import HMMERProfile

from ..profile import Profile
from .model import (
    StandardAltModel,
    StandardNode,
    StandardNullModel,
    StandardSpecialNode,
    Transitions,
)
from .result import StandardSearchResult


class StandardProfile(Profile):
    def __init__(
        self,
        alphabet: Alphabet,
        null_lprobs: Sequence[float],
        nodes_trans: Sequence[Tuple[StandardNode, Transitions]],
    ):
        super().__init__(alphabet)
        R = NormalState(b"R", alphabet, null_lprobs)
        self._null_model = StandardNullModel(R)

        special_node = StandardSpecialNode(
            S=MuteState(b"S", alphabet),
            N=NormalState(b"N", alphabet, null_lprobs),
            B=MuteState(b"B", alphabet),
            E=MuteState(b"E", alphabet),
            J=NormalState(b"J", alphabet, null_lprobs),
            C=NormalState(b"C", alphabet, null_lprobs),
            T=MuteState(b"T", alphabet),
        )

        self._alt_model = StandardAltModel(special_node, nodes_trans)
        self._set_fragment_length()

    @property
    def null_model(self) -> StandardNullModel:
        return self._null_model

    @property
    def alt_model(self) -> StandardAltModel:
        return self._alt_model

    def search(self, seq: CSequence) -> StandardSearchResult:
        self._set_target_length(seq.length)
        score0 = self.null_model.likelihood(seq)
        score1, path = self.alt_model.viterbi(seq)
        score = score1 - score0
        return StandardSearchResult(score, seq, path)


def create_profile(reader: HMMERProfile) -> StandardProfile:

    alphabet = Alphabet(reader.alphabet.encode(), b"X")

    prob_list = _create_probability_list(alphabet.symbols)
    null_lprobs = prob_list(reader.insert(0))

    nodes_trans: List[Tuple[StandardNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        M = NormalState(f"M{m}".encode(), alphabet, prob_list(reader.match(m)))
        I = NormalState(f"I{m}".encode(), alphabet, prob_list(reader.insert(m)))
        D = MuteState(f"D{m}".encode(), alphabet)

        node = StandardNode(M, I, D)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return StandardProfile(alphabet, null_lprobs, nodes_trans)


def _create_probability_list(symbols: bytes):
    def probability_list(lprob_table: Dict[str, Any]):
        probs = []
        for i in range(len(symbols)):
            key = symbols[i : i + 1].decode()
            probs.append(lprob_table.get(key, LPROB_ZERO))

        return lprob_normalize(probs)

    return probability_list
