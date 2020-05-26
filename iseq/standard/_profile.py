from typing import List, TypeVar

from hmmer_reader import HMMERModel
from imm import (
    Alphabet,
    Interval,
    MuteState,
    NormalState,
    Path,
    SequenceABC,
    lprob_normalize,
)

from .._hmmdata import HMMData
from .._model import EntryDistr, Node, Transitions
from .._profile import Profile
from ._typing import (
    StandardAltModel,
    StandardFragment,
    StandardNode,
    StandardNullModel,
    StandardSearchResults,
    StandardSpecialNode,
    StandardStep,
)

__all__ = [
    "StandardProfile",
    "create_standard_profile",
]


TAlphabet = TypeVar("TAlphabet", bound=Alphabet)


class StandardProfile(Profile[TAlphabet, NormalState]):
    def __init__(
        self,
        alphabet: TAlphabet,
        null_log_odds: List[float],
        core_nodes: List[Node],
        core_trans: List[Transitions],
        entry_distr: EntryDistr,
        hmmer3_compat=False,
    ):
        R = NormalState.create(b"R", alphabet, null_log_odds)
        null_model = StandardNullModel(R)

        special_node = StandardSpecialNode(
            S=MuteState.create(b"S", alphabet),
            N=NormalState.create(b"N", alphabet, null_log_odds),
            B=MuteState.create(b"B", alphabet),
            E=MuteState.create(b"E", alphabet),
            J=NormalState.create(b"J", alphabet, null_log_odds),
            C=NormalState.create(b"C", alphabet, null_log_odds),
            T=MuteState.create(b"T", alphabet),
        )

        alt_model = StandardAltModel(special_node, core_nodes, core_trans, entry_distr)
        super().__init__(alphabet, null_model, alt_model, hmmer3_compat)

    def search(
        self, sequence: SequenceABC[TAlphabet], window_length: int = 0
    ) -> StandardSearchResults:

        self._set_target_length_model(len(sequence))

        alt_results = self._alt_model.viterbi(sequence, window_length)

        def create_fragment(
            seq: SequenceABC[TAlphabet], path: Path[StandardStep], homologous: bool
        ):
            return StandardFragment(seq, path, homologous)

        search_results = StandardSearchResults(sequence, create_fragment)

        for alt_result in alt_results:
            subseq = alt_result.sequence
            score0 = self._null_model.likelihood(subseq)
            score1 = alt_result.loglikelihood
            score = score1 - score0
            window = Interval(subseq.start, subseq.start + len(subseq))
            if self._hmmer3_compat:
                score1 -= 3
            search_results.append(score, window, alt_result.path, score1)

        return search_results


def create_standard_profile(reader: HMMERModel) -> StandardProfile:

    alphabet = Alphabet.create(reader.alphabet.encode(), b"X")

    null_lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()

    nodes: List[StandardNode] = []
    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = NormalState.create(f"M{m}".encode(), alphabet, lprobs)

        lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        I = NormalState.create(f"I{m}".encode(), alphabet, lprobs)
        D = MuteState.create(f"D{m}".encode(), alphabet)

        nodes.append(StandardNode(M, I, D))

    trans: List[Transitions] = []
    for m in range(0, reader.M + 1):
        t = Transitions(**reader.trans(m))
        t.normalize()
        trans.append(t)

    return StandardProfile(alphabet, null_lprobs, nodes, trans, EntryDistr.UNIFORM)


def create_hmmer3_profile(hmm: HMMData, hmmer3_compat: bool = False) -> StandardProfile:
    null_lprobs = hmm.null_lprobs
    null_log_odds = [0.0] * len(null_lprobs)

    nodes: List[StandardNode] = []
    for m in range(1, hmm.model_length + 1):
        lodds = [v0 - v1 for v0, v1 in zip(hmm.match_lprobs(m), null_lprobs)]
        M = NormalState.create(f"M{m}".encode(), hmm.alphabet, lodds)
        I = NormalState.create(f"I{m}".encode(), hmm.alphabet, null_log_odds)
        D = MuteState.create(f"D{m}".encode(), hmm.alphabet)

        nodes.append(StandardNode(M, I, D))

    trans = hmm.transitions

    return StandardProfile(
        hmm.alphabet, null_log_odds, nodes, trans, EntryDistr.OCCUPANCY, hmmer3_compat
    )
