from typing import List, TypeVar

from imm import Alphabet, Interval, MuteState, NormalState, Path, Sequence, SequenceABC
from nmm import DNAAlphabet, NTTranslator, NullTranslator, RNAAlphabet

from iseq.hmmer_model import HMMERModel
from iseq.model import EntryDistr, Node, Transitions
from iseq.profile import Profile, ProfileID

from .typing import (
    HMMER3AltModel,
    HMMER3Fragment,
    HMMER3Node,
    HMMER3NullModel,
    HMMER3SearchResults,
    HMMER3SpecialNode,
    HMMER3Step,
)

__all__ = [
    "HMMER3Profile",
    "create_profile",
]


TAlphabet = TypeVar("TAlphabet", bound=Alphabet)


class HMMER3Profile(Profile[TAlphabet, NormalState]):
    def __init__(
        self,
        profid: ProfileID,
        alphabet: TAlphabet,
        null_log_odds: List[float],
        core_nodes: List[Node],
        core_trans: List[Transitions],
        entry_distr: EntryDistr,
        hmmer3_compat=False,
    ):
        R = NormalState.create(b"R", alphabet, null_log_odds)
        null_model = HMMER3NullModel.create(R)

        special_node = HMMER3SpecialNode(
            S=MuteState.create(b"S", alphabet),
            N=NormalState.create(b"N", alphabet, null_log_odds),
            B=MuteState.create(b"B", alphabet),
            E=MuteState.create(b"E", alphabet),
            J=NormalState.create(b"J", alphabet, null_log_odds),
            C=NormalState.create(b"C", alphabet, null_log_odds),
            T=MuteState.create(b"T", alphabet),
        )

        alt_model = HMMER3AltModel.create(
            special_node, core_nodes, core_trans, entry_distr
        )
        super().__init__(profid, alphabet, null_model, alt_model, hmmer3_compat)

        if isinstance(alphabet, (DNAAlphabet, RNAAlphabet)):
            self._translator = NTTranslator()
        else:
            self._translator = NullTranslator()

    @property
    def window_length(self) -> int:
        return super().window_length

    @window_length.setter
    def window_length(self, length: int) -> None:
        if length < -1:
            raise ValueError("Length must be greater than or equal to -1.")

        if length == -1:
            length = 2 * self._alt_model.core_length
        self._window_length = length

    def create_sequence(self, sequence: bytes) -> Sequence:
        seq = self._translator.translate(sequence, self.alphabet)
        return Sequence.create(seq, self.alphabet)

    def search(self, sequence: SequenceABC[TAlphabet]) -> HMMER3SearchResults:

        self._set_target_length_model(len(sequence))

        alt_results = self._alt_model.viterbi(sequence, self.window_length)

        def create_fragment(
            seq: SequenceABC[TAlphabet], path: Path[HMMER3Step], homologous: bool
        ):
            return HMMER3Fragment(seq, path, homologous)

        search_results = HMMER3SearchResults(sequence, create_fragment)

        for alt_result in alt_results:
            subseq = alt_result.sequence
            viterbi_score0 = self._null_model.likelihood(subseq)
            viterbi_score1 = alt_result.loglikelihood
            score = viterbi_score1 - viterbi_score0
            window = Interval(subseq.start, subseq.start + len(subseq))
            if self._hmmer3_compat:
                viterbi_score1 -= 3
            search_results.append(
                score, window, alt_result.path, viterbi_score1, viterbi_score0
            )

        return search_results


def create_profile(
    hmm: HMMERModel,
    hmmer3_compat: bool = False,
    entry_distr: EntryDistr = EntryDistr.OCCUPANCY,
    window_length: int = 0,
) -> HMMER3Profile:
    null_lprobs = hmm.null_lprobs
    null_log_odds = [0.0] * len(null_lprobs)

    nodes: List[HMMER3Node] = []
    for m in range(1, hmm.model_length + 1):
        lodds = [v0 - v1 for v0, v1 in zip(hmm.match_lprobs(m), null_lprobs)]
        M = NormalState.create(f"M{m}".encode(), hmm.alphabet, lodds)
        I = NormalState.create(f"I{m}".encode(), hmm.alphabet, null_log_odds)
        D = MuteState.create(f"D{m}".encode(), hmm.alphabet)

        nodes.append(HMMER3Node(M, I, D))

    trans = hmm.transitions

    profid = ProfileID(hmm.model_id.name, hmm.model_id.acc)
    prof = HMMER3Profile(
        profid,
        hmm.alphabet,
        null_log_odds,
        nodes,
        trans,
        entry_distr,
        hmmer3_compat,
    )
    prof.window_length = window_length
    return prof
