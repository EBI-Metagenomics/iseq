from typing import List, Tuple

from hmmer_reader import HMMERProfile
from nmm.alphabet import Alphabet
from nmm.path import Path, Step
from nmm.prob import lprob_normalize
from nmm.sequence import SequenceABC
from nmm.state import MuteState, NormalState
from nmm import Interval

from ._fragment import Fragment
from ._model import AltModel, Node, NullModel, SpecialNode, Transitions
from ._profile import Profile
from ._result import SearchResults
from ._typing import TAlphabet, MutableStep

StandardFragment = Fragment[TAlphabet, NormalState]
StandardStep = Step[MutableStep[NormalState]]
StandardSearchResults = SearchResults[TAlphabet, NormalState]

StandardNode = Node[NormalState]
StandardSpecialNode = SpecialNode[NormalState]
StandardNullModel = NullModel[NormalState]
StandardAltModel = AltModel[NormalState]

__all__ = [
    "StandardAltModel",
    "StandardFragment",
    "StandardNullModel",
    "StandardProfile",
    "StandardStep",
    "create_profile",
]


class StandardProfile(Profile[TAlphabet, NormalState]):
    def __init__(
        self,
        alphabet: TAlphabet,
        null_lprobs: List[float],
        nodes_trans: List[Tuple[StandardNode, Transitions]],
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

    def search(
        self, sequence: SequenceABC[TAlphabet], window_length: int = 0
    ) -> StandardSearchResults:

        self._set_target_length(len(sequence))
        alt_results = self.alt_model.viterbi(sequence, window_length)

        def create_fragment(
            seq: SequenceABC[TAlphabet], path: Path[StandardStep], homologous: bool
        ):
            return StandardFragment(seq, path, homologous)

        search_results = StandardSearchResults(sequence, create_fragment)

        for alt_result in alt_results:
            subseq = alt_result.sequence
            score0 = self.null_model.likelihood(subseq)
            score1 = alt_result.loglikelihood
            score = score1 - score0
            window = Interval(subseq.start, subseq.start + len(subseq))
            search_results.append(score, window, alt_result.path)

        return search_results


def create_profile(reader: HMMERProfile) -> StandardProfile:

    alphabet = Alphabet.create(reader.alphabet.encode(), b"X")

    null_lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()

    nodes_trans: List[Tuple[StandardNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = NormalState(f"M{m}".encode(), alphabet, lprobs)

        lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        I = NormalState(f"I{m}".encode(), alphabet, lprobs)
        D = MuteState(f"D{m}".encode(), alphabet)

        node = StandardNode(M, I, D)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return StandardProfile(alphabet, null_lprobs, nodes_trans)
