from typing import List, Tuple
from math import log

from hmmer_reader import HMMERProfile
from nmm.alphabet import Alphabet, AminoAlphabet
from nmm.path import Path, Step
from nmm.prob import lprob_normalize, lprob_zero
from nmm.sequence import SequenceABC
from nmm.state import MuteState, NormalState
from nmm import Interval

from ._fragment import Fragment
from ._model import AltModel, Node, NullModel, SpecialNode, Transitions, MSVModel
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
StandardMSVModel = MSVModel[NormalState]

__all__ = [
    "StandardMSVModel",
    "StandardAltModel",
    "StandardFragment",
    "StandardNullModel",
    "StandardProfile",
    "StandardStep",
    "create_standard_profile",
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
        self._null_model = StandardNullModel(R, self._special_transitions)

        self._special_node = StandardSpecialNode(
            S=MuteState(b"S", alphabet),
            N=NormalState(b"N", alphabet, null_lprobs),
            B=MuteState(b"B", alphabet),
            E=MuteState(b"E", alphabet),
            J=NormalState(b"J", alphabet, null_lprobs),
            C=NormalState(b"C", alphabet, null_lprobs),
            T=MuteState(b"T", alphabet),
        )

        self._alt_model = StandardAltModel(
            self._special_node, nodes_trans, self._special_transitions
        )
        self._alt_model.set_fragment_length()

        self._msv_model = StandardMSVModel(
            self._special_node, nodes_trans, self._special_transitions
        )
        self._msv_model.set_fragment_length()

    @property
    def null_model(self) -> StandardNullModel:
        return self._null_model

    @property
    def alt_model(self) -> StandardAltModel:
        return self._alt_model

    def search(
        self, sequence: SequenceABC[TAlphabet], window_length: int = 0
    ) -> StandardSearchResults:

        from time import time

        self._set_special_transitions(len(sequence))
        self._alt_model.update_special_transitions()
        self._msv_model.update_special_transitions()
        self._null_model.update_special_transitions()

        # self._alt_model._hmm.view()
        # self._null_model._hmm.view()

        start = time()
        self._msv_model.viterbi(sequence, window_length)
        print(f"MSV elapsed: {time()-start} seconds")

        start = time()
        alt_results = self.alt_model.viterbi(sequence, window_length)
        print(f"ALT elapsed: {time()-start} seconds")

        def create_fragment(
            seq: SequenceABC[TAlphabet], path: Path[StandardStep], homologous: bool
        ):
            return StandardFragment(seq, path, homologous)

        search_results = StandardSearchResults(sequence, create_fragment)

        for alt_result in alt_results:
            subseq = alt_result.sequence
            score0 = self.null_model.likelihood(subseq)
            score1 = alt_result.loglikelihood
            print(f"ALT loglikelihood: {score1:.12f}")
            score = score1 - score0
            window = Interval(subseq.start, subseq.start + len(subseq))
            search_results.append(score, window, alt_result.path)

        return search_results


def create_standard_profile(reader: HMMERProfile) -> StandardProfile:

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


def create_hmmer3_profile(reader: HMMERProfile) -> StandardProfile:

    # alphabet = Alphabet.create(reader.alphabet.encode(), b"X")
    alphabet = AminoAlphabet.create(reader.alphabet.encode(), b"X")
    null_lprobs = lprob_normalize(_hmmer3_null_amino_frequences(alphabet))
    one_lprobs = [0.0] * len(null_lprobs)

    # if isinstance(alphabet, AminoAlphabet):
    #     _hmmer3_null_amino_frequences(alphabet)
    # else:
    #     raise NotImplementedError()

    # null_lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()

    nodes_trans: List[Tuple[StandardNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())) - null_lprobs
        # lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = NormalState(f"M{m}".encode(), alphabet, lprobs.tolist())

        # lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        # I = NormalState(f"I{m}".encode(), alphabet, lprobs)
        I = NormalState(f"I{m}".encode(), alphabet, one_lprobs)
        D = MuteState(f"D{m}".encode(), alphabet)

        node = StandardNode(M, I, D)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    prof = StandardProfile(alphabet, one_lprobs, nodes_trans)
    # prof.alt_model.calculate_occupancy()
    return prof


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
