from typing import TypeVar, Union, Tuple, Sequence, List

from hmmer_reader import HMMERProfile
from nmm.prob import lprob_normalize
from nmm.alphabet import CAlphabet, Alphabet
from nmm.path import Path, Step
from nmm.state import MuteState, NormalState
from nmm.sequence import CSequence

from .fragment import Fragment
from .model import AltModel, Node, NullModel, SpecialNode, Transitions
from .result import SearchResult
from .profile import Profile

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)

StandardFragment = Fragment[TAlphabet, Union[NormalState, MuteState]]
StandardStep = Step[Union[NormalState, MuteState]]
StandardPath = Path[StandardStep]
StandardSearchResult = SearchResult[TAlphabet, Union[NormalState, MuteState]]

StandardNode = Node[NormalState]
StandardSpecialNode = SpecialNode[NormalState]
StandardNullModel = NullModel[NormalState]
StandardAltModel = AltModel[NormalState]


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
