from math import log
from typing import TypeVar, Union, Tuple, Sequence, List

from hmmer_reader import HMMERProfile
from nmm import GeneticCode
from nmm.prob import lprob_normalize, AminoTable
from nmm.alphabet import Alphabet, BaseAlphabet, AminoAlphabet
from nmm.path import Path, Step
from nmm.state import MuteState, FrameState, CodonState
from nmm.sequence import SequenceABC
from nmm.codon import Codon

from .fragment import Fragment
from .model import AltModel, Node, NullModel, SpecialNode, Transitions
from .result import SearchResult
from .profile import Profile
from .codon import CodonFragment, CodonPath

FrameStep = Step[Union[FrameState, MuteState]]
FramePath = Path[FrameStep]

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
        steps: List[CodonState] = []

        start: int = 0
        seq = self.sequence
        for step in self._path:
            if isinstance(step.state, MuteState):
                mstate = MuteState(step.state.name, step.state.alphabet)
                steps.append((mstate, 0))
            elif isinstance(step.state, FrameState):

                subseq = seq[start : start + step.seq_len]
                codon = Codon(b"XXX", self._alphabet)
                step.state.decode(subseq, codon)
                codons.append(codon)

                name = step.state.name
                abc = step.state.alphabet
                cstate = CodonState(name, abc, {codon: log(1.0)})
                steps.append((cstate, 3))

            else:
                raise RuntimeError("Wrong state type.")

            start += step.seq_len

        return CodonFragment(codons, CodonPath(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


# class FrameSearchResult(SearchResult):
#     pass


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


# class FrameBaseTable(BaseTable):
#     def __init__(self, codon_prob: CodonProb):
#         from numpy import logaddexp

#         base_abc = codon_prob.base
#         base_lprob = {base: LPROB_ZERO for base in base_abc.symbols}
#         norm = log(3)
#         for codon in codon_iter(base_abc):
#             lprob = codon_prob.get_lprob(codon)
#             triplet = codon.symbols

#             base_lprob[triplet[0]] = logaddexp(base_lprob[triplet[0]], lprob - norm)
#             base_lprob[triplet[1]] = logaddexp(base_lprob[triplet[1]], lprob - norm)
#             base_lprob[triplet[2]] = logaddexp(base_lprob[triplet[2]], lprob - norm)

#         super().__init__(base_abc, [base_lprob[base] for base in base_abc.symbols])


# class FrameCodonProb(CodonProb):
#     def __init__(
#         self, base_abc: CBaseAlphabet, aminot: AminoTable, gencode: GeneticCode
#     ):
#         from numpy import logaddexp

#         super().__init__(base_abc)
#         self._lprob: Dict[Codon, float] = {}

#         codon_lprobs = []
#         lprob_norm = LPROB_ZERO
#         for i in range(len(aminot.alphabet.symbols)):
#             aa = aminot.alphabet.symbols[i : i + 1]
#             lprob = aminot.lprob(aa)

#             codons = gencode.codons(aa)
#             if len(codons) == 0:
#                 continue

#             norm = log(len(codons))
#             for codon in codons:
#                 codon_lprobs.append((codon, lprob - norm))
#                 lprob_norm = logaddexp(lprob_norm, codon_lprobs[-1][1])

#         for codon, lprob in codon_lprobs:
#             self.set_lprob(codon, lprob - lprob_norm)


# class FrameStateFactory:
#     def __init__(
#         self,
#         base_abc: BaseAlphabet,
#         amino_abc: Alphabet,
#         gcode: GeneticCode,
#         epsilon: float,
#     ):
#         self._base_abc = base_abc
#         self._amino_abc = amino_abc
#         self._gcode = gcode
#         self._epsilon = epsilon

#     def create(self, name: bytes, aminot: AminoTable) -> FrameState:
#         codonp = FrameCodonProb(self._base_abc, aminot, self._gcode)
#         baset = FrameBaseTable(codonp)
#         codont = CodonTable(codonp)
#         return FrameState(name, baset, codont, self._epsilon)

#     @property
#     def base_alphabet(self) -> BaseAlphabet:
#         return self._base_abc

#     @property
#     def genetic_code(self) -> GeneticCode:
#         return self._gcode

#     @property
#     def epsilon(self) -> float:
#         return self._epsilon


# class FrameProfile(Profile[FrameSearchResult, Step[Union[FrameState, MuteState]]]):
#     def __init__(
#         self,
#         factory: FrameStateFactory,
#         base_abc: BaseAlphabet,
#         null_aminot: AminoTable,
#         nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
#     ):
#         super().__init__(base_abc)

#         R = factory.create(b"R", null_aminot)
#         self._null_model = FrameNullModel(R)

#         special_node = FrameSpecialNode(
#             S=MuteState(b"S", factory.base_alphabet),
#             N=factory.create(b"N", null_aminot),
#             B=MuteState(b"B", factory.base_alphabet),
#             E=MuteState(b"E", factory.base_alphabet),
#             J=factory.create(b"J", null_aminot),
#             C=factory.create(b"C", null_aminot),
#             T=MuteState(b"T", factory.base_alphabet),
#         )

#         self._alt_model = FrameAltModel(special_node, nodes_trans)
#         self._set_fragment_length()

#     @property
#     def null_model(self) -> FrameNullModel:
#         return self._null_model

#     @property
#     def alt_model(self) -> FrameAltModel:
#         return self._alt_model

#     def search(self, sequence: CSequence) -> FrameSearchResult:
#         score, path = self._search(sequence)
#         return FrameSearchResult(self._alphabet, score, sequence, path)


# def create_profile(reader: HMMERProfile, epsilon: float = 0.1) -> FrameProfile:

#     base_abc = BaseAlphabet(b"ACGU", b"X")
#     amino_abc = AminoAlphabet(reader.alphabet.encode(), b"X")

#     lprobs = lprob_normalize(list(reader.insert(0).values())).tolist()
#     null_aminot = AminoTable(amino_abc, lprobs)
#     ffact = FrameStateFactory(
#         base_abc, amino_abc, GeneticCode(base_abc, amino_abc), epsilon
#     )

#     nodes_trans: List[Tuple[FrameNode, Transitions]] = []

#     for m in range(1, reader.M + 1):
#         lprobs = lprob_normalize(list(reader.match(m).values())).tolist()
#         M = ffact.create(f"M{m}".encode(), AminoTable(amino_abc, lprobs))

#         lprobs = lprob_normalize(list(reader.insert(m).values())).tolist()
#         I = ffact.create(f"I{m}".encode(), AminoTable(amino_abc, lprobs))

#         D = MuteState(f"D{m}".encode(), base_abc)

#         node = FrameNode(M, I, D,)

#         trans = Transitions(**reader.trans(m - 1))
#         trans.normalize()

#         nodes_trans.append((node, trans))

#     return FrameProfile(ffact, base_abc, null_aminot, nodes_trans)
