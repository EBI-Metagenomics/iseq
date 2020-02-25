from typing import List, Tuple, TypeVar, Union

from hmmer_reader import HMMERProfile
from nmm.alphabet import Alphabet
from nmm.path import Path, Step
from nmm.prob import lprob_normalize
from nmm.sequence import Sequence, SequenceABC
from nmm.state import MuteState, NormalState

from ._fragment import Fragment
from ._model import AltModel, Node, NullModel, SpecialNode, Transitions
from ._profile import Profile
from ._result import SearchResult

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

StandardFragment = Fragment[TAlphabet, NormalState]
StandardStep = Step[Union[NormalState, MuteState]]
StandardPath = Path[StandardStep]
StandardSearchResult = SearchResult[TAlphabet, NormalState]

StandardNode = Node[NormalState]
StandardSpecialNode = SpecialNode[NormalState]
StandardNullModel = NullModel[NormalState]
StandardAltModel = AltModel[NormalState]

__all__ = [
    "create_profile",
    "StandardProfile",
    "StandardFragment",
    "StandardStep",
    "StandardPath",
    "StandardSearchResult",
]


class StandardProfile(Profile[TAlphabet, NormalState]):
    def __init__(
        self,
        alphabet: Alphabet,
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
        self, sequence: Sequence[TAlphabet], window: int = 0
    ) -> List[StandardSearchResult]:
        score_path_list = self._search(sequence, window)
        return [
            StandardSearchResult(score, sequence, path, _create_fragment)
            for score, path in score_path_list
        ]


def _create_fragment(
    sequence: SequenceABC[TAlphabet], path: Path[StandardStep], homologous: bool
):
    return StandardFragment(sequence, path, homologous)


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
