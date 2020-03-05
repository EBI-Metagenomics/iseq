from math import log
from typing import List

from hmmer_reader import HMMERModel
from nmm import Interval
from nmm.alphabet import Alphabet, AminoAlphabet
from nmm.path import Path
from nmm.prob import lprob_normalize, lprob_zero
from nmm.sequence import SequenceABC
from nmm.state import MuteState, NormalState

from .._model import EntryDistr, Node, Transitions
from .._profile import Profile
from .._typing import TAlphabet
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
        R = NormalState(b"R", alphabet, null_log_odds)
        null_model = StandardNullModel(R)

        special_node = StandardSpecialNode(
            S=MuteState(b"S", alphabet),
            N=NormalState(b"N", alphabet, null_log_odds),
            B=MuteState(b"B", alphabet),
            E=MuteState(b"E", alphabet),
            J=NormalState(b"J", alphabet, null_log_odds),
            C=NormalState(b"C", alphabet, null_log_odds),
            T=MuteState(b"T", alphabet),
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
            search_results.append(score, window, alt_result.path, score1)

        return search_results


def create_standard_profile(reader: HMMERModel) -> StandardProfile:

    alphabet = Alphabet.create(reader.alphabet.encode(), b"X")

    null_lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()

    nodes: List[StandardNode] = []
    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = NormalState(f"M{m}".encode(), alphabet, lprobs)

        lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        I = NormalState(f"I{m}".encode(), alphabet, lprobs)
        D = MuteState(f"D{m}".encode(), alphabet)

        nodes.append(StandardNode(M, I, D))

    trans: List[Transitions] = []
    for m in range(0, reader.M + 1):
        t = Transitions(**reader.trans(m))
        t.normalize()
        trans.append(t)

    return StandardProfile(alphabet, null_lprobs, nodes, trans, EntryDistr.UNIFORM)


def create_hmmer3_profile(
    reader: HMMERModel, hmmer3_compat: bool = False
) -> StandardProfile:
    alphabet = AminoAlphabet.create(reader.alphabet.encode(), b"X")
    null_lprobs = _hmmer3_null_amino_frequences(alphabet)
    null_log_odds = [0.0] * len(null_lprobs)

    nodes: List[StandardNode] = []
    for m in range(1, reader.M + 1):
        lodds = [v0 - v1 for v0, v1 in zip(reader.match(m).values(), null_lprobs)]
        M = NormalState(f"M{m}".encode(), alphabet, lodds)
        I = NormalState(f"I{m}".encode(), alphabet, null_log_odds)
        D = MuteState(f"D{m}".encode(), alphabet)

        nodes.append(StandardNode(M, I, D))

    trans: List[Transitions] = []
    for m in range(0, reader.M + 1):
        t = Transitions(**reader.trans(m))
        trans.append(t)

    return StandardProfile(
        alphabet, null_log_odds, nodes, trans, EntryDistr.OCCUPANCY, hmmer3_compat
    )


def _hmmer3_null_amino_frequences(alphabet: AminoAlphabet):
    """
    Copy/paste from HMMER3 amino acid frequences infered form Swiss-Prot 50.8,
    (Oct 2006), counting over 85956127 (86.0M) residues.
    """
    lprobs = {
        b"A"[0]: log(0.0787945),
        b"C"[0]: log(0.0151600),
        b"D"[0]: log(0.0535222),
        b"E"[0]: log(0.0668298),
        b"F"[0]: log(0.0397062),
        b"G"[0]: log(0.0695071),
        b"H"[0]: log(0.0229198),
        b"I"[0]: log(0.0590092),
        b"K"[0]: log(0.0594422),
        b"L"[0]: log(0.0963728),
        b"M"[0]: log(0.0237718),
        b"N"[0]: log(0.0414386),
        b"P"[0]: log(0.0482904),
        b"Q"[0]: log(0.0395639),
        b"R"[0]: log(0.0540978),
        b"S"[0]: log(0.0683364),
        b"T"[0]: log(0.0540687),
        b"V"[0]: log(0.0673417),
        b"W"[0]: log(0.0114135),
        b"Y"[0]: log(0.0304133),
    }
    abc = list(alphabet.symbols)
    return [lprobs.get(sym, lprob_zero()) for sym in abc]
