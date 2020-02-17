from math import log
from typing import Any, Dict, List, Sequence, Tuple

from hmmer_reader import HMMERProfile

from nmm import (
    Alphabet,
    Base,
    BaseTable,
    CodonTable,
    GeneticCode,
    LPROB_ZERO,
    FrameState,
    MuteState,
    AlphabetTable,
    lprob_normalize,
)
from .result import FrameSearchResult
from .model import (
    FrameAltModel,
    FrameNode,
    FrameNullModel,
    FrameSpecialNode,
    Transitions,
)
from ..profile import Profile


class FrameStateFactory:
    def __init__(
        self, base: Base, prot_abc: Alphabet, gcode: GeneticCode, epsilon: float,
    ):
        self._base = base
        self._prot_abc = prot_abc
        self._gcode = gcode
        self._epsilon = epsilon

    def create(self, name: bytes, prot_abct: AlphabetTable) -> FrameState:
        codon_lprobs = _infer_codon_lprobs(prot_abct, self._gcode)
        base_lprobs = _infer_base_lprobs(codon_lprobs, self._alphabet)
        base_table = BaseTable.create(self._alphabet, base_lprobs)
        codon_table = CodonTable.create(self._alphabet, codon_lprobs)
        return FrameState(name, base_table, codon_table, self._epsilon)

    @property
    def bases(self) -> Alphabet:
        return self._alphabet

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
        aa_lprobs: Dict[bytes, float],
        nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
    ):
        super().__init__()

        R = fstate_factory.create(b"R", aa_lprobs)
        self._null_model = FrameNullModel(R)

        special_node = FrameSpecialNode(
            S=MuteState(b"S", fstate_factory.bases),
            N=fstate_factory.create(b"N", aa_lprobs),
            B=MuteState(b"B", fstate_factory.bases),
            E=MuteState(b"E", fstate_factory.bases),
            J=fstate_factory.create(b"J", aa_lprobs),
            C=fstate_factory.create(b"C", aa_lprobs),
            T=MuteState(b"T", fstate_factory.bases),
        )

        self._alt_model = FrameAltModel(special_node, nodes_trans)
        self._set_fragment_length()

    @property
    def null_model(self) -> FrameNullModel:
        return self._null_model

    @property
    def alt_model(self) -> FrameAltModel:
        return self._alt_model

    def search(self, seq: bytes) -> FrameSearchResult:
        self._set_target_length(len(seq))
        score0 = self.null_model.likelihood(seq)
        score1, path = self.alt_model.viterbi(seq)
        score = score1 - score0
        return FrameSearchResult(score, seq, path)


def create_profile(reader: HMMERProfile, epsilon: float = 0.1) -> FrameProfile:

    base = Base(Alphabet(b"ACGU", b"X"))
    breakpoint()
    prot = Alphabet(reader.alphabet.encode(), b"X")

    prob_list = _create_probability_list(prot.symbols)

    null_lprobs = prob_list(reader.insert(0))
    ffact = FrameStateFactory(base, prot, GeneticCode(base), epsilon)

    nodes_trans: List[Tuple[FrameNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        M = ffact.create(
            f"M{m}".encode(), AlphabetTable(prot, prob_list(reader.match(m)))
        )
        I = ffact.create(
            f"I{m}".encode(), AlphabetTable(prot, prob_list(reader.insert(m)))
        )
        D = MuteState(f"D{m}".encode(), base)

        node = FrameNode(M, I, D,)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return FrameProfile(ffact, null_lprobs, nodes_trans)


def _infer_codon_lprobs(aa_lprobs: Dict[bytes, float], gencode: GeneticCode):
    from numpy import logaddexp

    codon_lprobs = []
    lprob_norm = LPROB_ZERO
    for aa, lprob in aa_lprobs.items():

        codons = gencode.codons(aa)
        if len(codons) == 0:
            continue

        norm = log(len(codons))
        for codon in codons:
            codon_lprobs.append((codon, lprob - norm))
            lprob_norm = logaddexp(lprob_norm, codon_lprobs[-1][1])

    codon_lprobs = [(i[0], i[1] - lprob_norm) for i in codon_lprobs]
    return dict(codon_lprobs)


def _infer_base_lprobs(codon_lprobs, alphabet: Alphabet):
    from scipy.special import logsumexp

    lprobs: Dict[Base, list] = {Base(sym): [] for sym in alphabet.symbols}
    lprob_norm = log(3)
    for codon, lprob in codon_lprobs.items():
        lprobs[codon.base(0)] += [lprob - lprob_norm]
        lprobs[codon.base(1)] += [lprob - lprob_norm]
        lprobs[codon.base(2)] += [lprob - lprob_norm]

    return {b: logsumexp(lp) for b, lp in lprobs.items()}


def _create_probability_list(symbols: bytes):
    def probability_list(lprob_table: Dict[str, Any]):
        probs = []
        for i in range(len(symbols)):
            key = symbols[i : i + 1].decode()
            probs.append(lprob_table.get(key, LPROB_ZERO))

        return lprob_normalize(probs).tolist()

    return probability_list
