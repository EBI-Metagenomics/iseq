from math import log
from typing import TypeVar, Union, Tuple, List

from hmmer_reader import HMMERProfile
from nmm import GeneticCode
from nmm.prob import (
    lprob_normalize,
    AminoTable,
    CodonProb,
    lprob_zero,
    BaseTable,
    CodonTable,
)
from nmm.alphabet import Alphabet, BaseAlphabet, AminoAlphabet
from nmm.path import Path, Step
from nmm.state import MuteState, FrameState, CodonState
from nmm.sequence import SequenceABC, Sequence
from nmm.codon import Codon, codon_iter

from .fragment import Fragment
from .model import AltModel, Node, NullModel, SpecialNode, Transitions
from .result import SearchResult
from .profile import Profile
from .codon import CodonFragment, CodonPath, CodonStep

FrameStep = Step[Union[FrameState, MuteState]]
FramePath = Path[FrameStep]
FrameSearchResult = SearchResult[BaseAlphabet, FrameState]

FrameNode = Node[FrameState]
FrameSpecialNode = SpecialNode[FrameState]
FrameNullModel = NullModel[FrameState]
FrameAltModel = AltModel[FrameState]


class FrameFragment(Fragment[BaseAlphabet, FrameState]):
    def __init__(
        self,
        sequence: SequenceABC[BaseAlphabet],
        path: Path[Step[FrameState]],
        homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    def decode(self) -> CodonFragment:
        codons: List[Codon] = []
        steps: List[CodonStep] = []

        start = 0
        seq = self.sequence
        abc = seq.alphabet
        for step in self._path:
            if isinstance(step.state, MuteState):

                mstate = MuteState(step.state.name, step.state.alphabet)
                steps.append(CodonStep.create(mstate, 0))

            elif isinstance(step.state, FrameState):

                subseq = seq[start : start + step.seq_len]
                codon = Codon.create(b"XXX", abc)
                step.state.decode(subseq, codon)
                codons.append(codon)

                name = step.state.name
                cstate = CodonState(name, abc, {codon: log(1.0)})
                steps.append(CodonStep.create(cstate, 3))

            else:
                raise RuntimeError("Wrong state type.")

            start += step.seq_len

        sequence = Sequence.create(b"".join(c.symbols for c in codons), abc)
        return CodonFragment(sequence, CodonPath.create(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


# #     def __init__(
# #         self,
# #         alphabet: CBaseAlphabet,
# #         loglik: float,
# #         sequence: SequenceABC,
# #         path: FramePath,
# #     ):
# #         self._loglik = loglik
# #         self._fragments: List[FrameFragment] = []
# #         self._intervals: List[Interval] = []

# #         steps = list(path)
# #         for fragi, stepi, homologous in self._create_fragments(path):
# #             substeps = steps[stepi.start : stepi.stop]
# #             fragment_path = FramePath([(s.state, s.seq_len) for s in substeps])
# #             seq = sequence[fragi]
# #             frag = FrameFragment(alphabet, seq, fragment_path, homologous)
# #             self._fragments.append(frag)
# #             self._intervals.append(fragi)

# #     @property
# #     def fragments(self) -> Sequence[FrameFragment]:
# #         return self._fragments

# #     @property
# #     def intervals(self) -> Sequence[Interval]:
# #         return self._intervals

# #     @property
# #     def loglikelihood(self) -> float:
# #         return self._loglik

# #     # def decode(self) -> CodonSearchResult:
# #     #     fragments: List[CodonFragment] = []
# #     #     intervals: List[Interval] = []

# #     #     start = end = 0
# #     #     for i, frag in enumerate(self._fragments):

# #     #         codon_frag = frag.decode()
# #     #         end += len(codon_frag.sequence)

# #     #         fragments.append(codon_frag)
# #     #         intervals.append(Interval(start, end))

# #     #         start = end

# #     #     return CodonSearchResult(self.score, fragments, intervals)


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


class FrameStateFactory:
    def __init__(
        self, gcode: GeneticCode, epsilon: float,
    ):
        self._gcode = gcode
        self._epsilon = epsilon

    def create(self, name: bytes, aminot: AminoTable) -> FrameState:
        codonp = _create_codon_prob(aminot, self._gcode)
        baset = _create_base_table(codonp)
        codont = CodonTable.create(codonp)
        return FrameState(name, baset, codont, self._epsilon)

    @property
    def genetic_code(self) -> GeneticCode:
        return self._gcode

    @property
    def epsilon(self) -> float:
        return self._epsilon


class FrameProfile(Profile[BaseAlphabet, FrameState]):
    def __init__(
        self,
        factory: FrameStateFactory,
        null_aminot: AminoTable,
        nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
    ):
        base_alphabet = factory.genetic_code.base_alphabet
        super().__init__(base_alphabet)

        R = factory.create(b"R", null_aminot)
        self._null_model = FrameNullModel(R)

        special_node = FrameSpecialNode(
            S=MuteState(b"S", base_alphabet),
            N=factory.create(b"N", null_aminot),
            B=MuteState(b"B", base_alphabet),
            E=MuteState(b"E", base_alphabet),
            J=factory.create(b"J", null_aminot),
            C=factory.create(b"C", null_aminot),
            T=MuteState(b"T", base_alphabet),
        )

        self._alt_model = FrameAltModel(special_node, nodes_trans)
        self._set_fragment_length()

    @property
    def null_model(self) -> FrameNullModel:
        return self._null_model

    @property
    def alt_model(self) -> FrameAltModel:
        return self._alt_model

    def search(self, sequence: Sequence[BaseAlphabet]) -> FrameSearchResult:
        score, path = self._search(sequence)
        return FrameSearchResult(score, sequence, path, _create_fragment)


def _create_fragment(
    sequence: SequenceABC[BaseAlphabet], path: Path[Step[FrameState]], homologous: bool
):
    return FrameFragment(sequence, path, homologous)


def create_profile(reader: HMMERProfile, epsilon: float = 0.1) -> FrameProfile:

    base_abc = BaseAlphabet.create(b"ACGU", b"X")
    amino_abc = AminoAlphabet.create(reader.alphabet.encode(), b"X")

    lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()
    null_aminot = AminoTable.create(amino_abc, lprobs)
    factory = FrameStateFactory(GeneticCode(base_abc, amino_abc), epsilon)

    nodes_trans: List[Tuple[FrameNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
        M = factory.create(f"M{m}".encode(), AminoTable.create(amino_abc, lprobs))

        lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
        I = factory.create(f"I{m}".encode(), AminoTable.create(amino_abc, lprobs))

        D = MuteState(f"D{m}".encode(), base_abc)

        node = FrameNode(M, I, D,)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return FrameProfile(factory, null_aminot, nodes_trans)
