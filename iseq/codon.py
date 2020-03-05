from math import log
from typing import List, Tuple, Union

from hmmer_reader import HMMERProfile
from nmm import GeneticCode
from nmm.alphabet import AminoAlphabet, BaseAlphabet
from nmm.codon import Codon, codon_iter
from nmm.path import Path, Step
from nmm.prob import (
    AminoTable,
    BaseTable,
    CodonProb,
    CodonTable,
    lprob_normalize,
    lprob_zero,
)
from nmm.sequence import Sequence, SequenceABC
from nmm.state import CodonState, CodonState, MuteState
from nmm import Interval

from ._codon import CodonFragment, CodonPath, CodonStep
from ._fragment import Fragment
from ._model import AltModel, Node, NullModel, SpecialNode, Transitions
from ._profile import Profile
from ._result import SearchResults

CodonStep = Step[Union[CodonState, MuteState]]
CodonSearchResults = SearchResults[BaseAlphabet, CodonState]

CodonNode = Node[CodonState]
CodonSpecialNode = SpecialNode[CodonState]
CodonNullModel = NullModel[CodonState]
CodonAltModel = AltModel[CodonState]

__all__ = [
    "CodonAltModel",
    "CodonFragment",
    "CodonNullModel",
    "CodonProfile",
    "create_profile",
]


class CodonFragment(Fragment[BaseAlphabet, CodonState]):
    def __init__(
        self,
        sequence: SequenceABC[BaseAlphabet],
        path: Path[CodonStep],
        homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    def decode(self) -> CodonFragment:
        codons: List[Codon] = []
        steps: List[CodonStep] = []

        start = 0
        for step in self.path:
            if isinstance(step.state, MuteState):

                mstate = MuteState(step.state.name, step.state.alphabet)
                steps.append(CodonStep.create(mstate, 0))

            elif isinstance(step.state, CodonState):

                subseq = self.sequence[start : start + step.seq_len]
                codon = step.state.decode(subseq)[1]
                codons.append(codon)

                name = step.state.name
                cstate = CodonState(name, self.sequence.alphabet, {codon: log(1.0)})
                steps.append(CodonStep.create(cstate, 3))

            else:
                raise RuntimeError("Wrong state type.")

            start += step.seq_len

        CodonSequence = Sequence[BaseAlphabet]
        abc = self.sequence.alphabet
        sequence = CodonSequence.create(b"".join(c.symbols for c in codons), abc)
        return CodonFragment(sequence, CodonPath.create(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


def _create_base_table(codonp: CodonProb):
    from numpy import logaddexp

    base_abc = codonp.alphabet
    base_lprob = {base: lprob_zero() for base in base_abc.symbols}
    norm = log(3)
    for codon in codon_iter(base_abc):
        lprob = codonp.get_lprob(codon)
        triplet = codon.symbols

        base_lprob[triplet[0]] = logaddexp(base_lprob[triplet[0]], lprob - norm)
        base_lprob[triplet[1]] = logaddexp(base_lprob[triplet[1]], lprob - norm)
        base_lprob[triplet[2]] = logaddexp(base_lprob[triplet[2]], lprob - norm)

    return BaseTable.create(base_abc, [base_lprob[base] for base in base_abc.symbols])


def _create_codon_prob(aminot: AminoTable, gencode: GeneticCode) -> CodonProb:
    from numpy import logaddexp

    codonp = CodonProb.create(gencode.base_alphabet)

    codon_lprobs = []
    lprob_norm = lprob_zero()
    for i in range(len(aminot.alphabet.symbols)):
        aa = aminot.alphabet.symbols[i : i + 1]
        lprob = aminot.lprob(aa)

        codons = gencode.codons(aa)
        if len(codons) == 0:
            continue

        norm = log(len(codons))
        for codon in codons:
            codon_lprobs.append((codon, lprob - norm))
            lprob_norm = logaddexp(lprob_norm, codon_lprobs[-1][1])

    for codon, lprob in codon_lprobs:
        codonp.set_lprob(codon, lprob - lprob_norm)

    return codonp


class CodonStateFactory:
    def __init__(
        self, gcode: GeneticCode,
    ):
        self._gcode = gcode

    def create(self, name: bytes, aminot: AminoTable) -> CodonState:
        codonp = _create_codon_prob(aminot, self._gcode)
        return CodonState(name, codonp)

    @property
    def genetic_code(self) -> GeneticCode:
        return self._gcode


class CodonProfile(Profile[BaseAlphabet, CodonState]):
    def __init__(
        self,
        factory: CodonStateFactory,
        null_aminot: AminoTable,
        nodes_trans: List[Tuple[CodonNode, Transitions]],
    ):
        base_alphabet = factory.genetic_code.base_alphabet
        super().__init__(base_alphabet)

        R = factory.create(b"R", null_aminot)
        self._null_model = CodonNullModel(R, self._special_transitions)

        self._special_node = CodonSpecialNode(
            S=MuteState(b"S", base_alphabet),
            N=factory.create(b"N", null_aminot),
            B=MuteState(b"B", base_alphabet),
            E=MuteState(b"E", base_alphabet),
            J=factory.create(b"J", null_aminot),
            C=factory.create(b"C", null_aminot),
            T=MuteState(b"T", base_alphabet),
        )

        self._alt_model = CodonAltModel(
            self._special_node, nodes_trans, self._special_transitions
        )
        self._alt_model.set_fragment_length()

    @property
    def null_model(self) -> CodonNullModel:
        return self._null_model

    @property
    def alt_model(self) -> CodonAltModel:
        return self._alt_model

    def search(
        self, sequence: SequenceABC[BaseAlphabet], window_length: int = 0
    ) -> CodonSearchResults:

        self._special_transitions(len(sequence))
        self._alt_model.set_special_transitions()
        self._null_model.set_special_transitions()

        alt_results = self.alt_model.viterbi(sequence, window_length)

        def create_fragment(
            seq: SequenceABC[BaseAlphabet], path: Path[CodonStep], homologous: bool
        ):
            return CodonFragment(seq, path, homologous)

        search_results = CodonSearchResults(sequence, create_fragment)

        for alt_result in alt_results:
            subseq = alt_result.sequence
            score0 = self.null_model.likelihood(subseq)
            score1 = alt_result.loglikelihood
            score = score1 - score0
            window = Interval(subseq.start, subseq.start + len(subseq))
            search_results.append(score, window, alt_result.path)

        return search_results


def create_profile(
    reader: HMMERProfile, base_abc: BaseAlphabet, epsilon: float = 0.1
) -> CodonProfile:

    amino_abc = AminoAlphabet.create(reader.alphabet.encode(), b"X")

    lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()
    null_aminot = AminoTable.create(amino_abc, lprobs)
    factory = CodonStateFactory(GeneticCode(base_abc, amino_abc))

    nodes_trans: List[Tuple[CodonNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = factory.create(f"M{m}".encode(), AminoTable.create(amino_abc, lprobs))

        lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        I = factory.create(f"I{m}".encode(), AminoTable.create(amino_abc, lprobs))

        D = MuteState(f"D{m}".encode(), base_abc)

        node = CodonNode(M, I, D,)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return CodonProfile(factory, null_aminot, nodes_trans)
