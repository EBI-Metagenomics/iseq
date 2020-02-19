from typing import List, Sequence, Tuple

from nmm import GeneticCode
from nmm.alphabet import Alphabet, AminoAlphabet, BaseAlphabet
from nmm.prob import CodonTable, lprob_normalize, AminoTable
from nmm.sequence import CSequence
from nmm.state import FrameState, MuteState

from hmmer_reader import HMMERProfile

from ..profile import Profile
from .base_table import FrameBaseTable
from .codon_prob import FrameCodonProb
from .model import (
    FrameAltModel,
    FrameNode,
    FrameNullModel,
    FrameSpecialNode,
    Transitions,
)
from .result import FrameSearchResult


class FrameStateFactory:
    def __init__(
        self,
        base_abc: BaseAlphabet,
        amino_abc: Alphabet,
        gcode: GeneticCode,
        epsilon: float,
    ):
        self._base_abc = base_abc
        self._amino_abc = amino_abc
        self._gcode = gcode
        self._epsilon = epsilon

    def create(self, name: bytes, aminot: AminoTable) -> FrameState:
        codonp = FrameCodonProb(self._base_abc, aminot, self._gcode)
        baset = FrameBaseTable(codonp)
        codont = CodonTable(codonp)
        return FrameState(name, baset, codont, self._epsilon)

    @property
    def base_alphabet(self) -> BaseAlphabet:
        return self._base_abc

    @property
    def genetic_code(self) -> GeneticCode:
        return self._gcode

    @property
    def epsilon(self) -> float:
        return self._epsilon


class FrameProfile(Profile):
    def __init__(
        self,
        fstate_factory: FrameStateFactory,
        base_abc: BaseAlphabet,
        null_aminot: AminoTable,
        nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
    ):
        super().__init__(base_abc)

        R = fstate_factory.create(b"R", null_aminot)
        self._null_model = FrameNullModel(R)

        special_node = FrameSpecialNode(
            S=MuteState(b"S", fstate_factory.base_alphabet),
            N=fstate_factory.create(b"N", null_aminot),
            B=MuteState(b"B", fstate_factory.base_alphabet),
            E=MuteState(b"E", fstate_factory.base_alphabet),
            J=fstate_factory.create(b"J", null_aminot),
            C=fstate_factory.create(b"C", null_aminot),
            T=MuteState(b"T", fstate_factory.base_alphabet),
        )

        self._alt_model = FrameAltModel(special_node, nodes_trans)
        self._set_fragment_length()

    @property
    def null_model(self) -> FrameNullModel:
        return self._null_model

    @property
    def alt_model(self) -> FrameAltModel:
        return self._alt_model

    def search(self, seq: CSequence) -> FrameSearchResult:
        self._set_target_length(seq.length)
        score0 = self.null_model.likelihood(seq)
        score1, path = self.alt_model.viterbi(seq)
        score = score1 - score0
        return FrameSearchResult(self._alphabet, score, seq, path)


def create_profile(reader: HMMERProfile, epsilon: float = 0.1) -> FrameProfile:

    base_abc = BaseAlphabet(b"ACGU", b"X")
    amino_abc = AminoAlphabet(reader.alphabet.encode(), b"X")

    lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()
    null_aminot = AminoTable(amino_abc, lprobs)
    ffact = FrameStateFactory(
        base_abc, amino_abc, GeneticCode(base_abc, amino_abc), epsilon
    )

    nodes_trans: List[Tuple[FrameNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = ffact.create(f"M{m}".encode(), AminoTable(amino_abc, lprobs))

        lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        I = ffact.create(f"I{m}".encode(), AminoTable(amino_abc, lprobs))

        D = MuteState(f"D{m}".encode(), base_abc)

        node = FrameNode(M, I, D,)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return FrameProfile(ffact, base_abc, null_aminot, nodes_trans)
