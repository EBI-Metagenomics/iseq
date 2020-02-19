from typing import TypeVar, Union, Tuple, Sequence, List

from hmmer_reader import HMMERProfile
from nmm.prob import lprob_normalize
from nmm.alphabet import CAlphabet, Alphabet
from nmm.path import Path, Step
from nmm.state import MuteState, FrameState
from nmm.sequence import CSequence

from .fragment import Fragment
from .model import AltModel, Node, NullModel, SpecialNode, Transitions
from .result import SearchResult
from .profile import Profile

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)

FrameStep = Step[Union[FrameState, MuteState]]
FramePath = Path[FrameStep]

FrameNode = Node[FrameState]
FrameSpecialNode = SpecialNode[FrameState]
FrameNullModel = NullModel[FrameState]
FrameAltModel = AltModel[FrameState]


class FrameFragment(Fragment[TAlphabet, Union[FrameState, MuteState]]):
    def __init__(
        self, sequence: Sequence[TAlphabet], path: Path[FrameStep], homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    # @property
    # def sequence(self) -> SequenceABC:
    #     return self._sequence

    # def items(self) -> Iterator[Tuple[bytes, FrameStep]]:
    #     start = end = 0
    #     for step in self._path:
    #         end += step.seq_len
    #         yield (self._sequence.symbols[start:end], step)
    #         start = end

    # def decode(self) -> CodonFragment:
    #     codons: List[Codon] = []
    #     steps = []
    #     # npath = CodonPath()

    #     start: int = 0
    #     seq = self.sequence
    #     for step in self._path:
    #         if isinstance(step.state, MuteState):
    #             mstate = MuteState(step.state.name, step.state.alphabet)
    #             # npath.append_codon_step(mstate, 0)
    #             steps.append((mstate, 0))
    #         else:
    #             assert isinstance(step.state, FrameState)

    #             subseq = seq[Interval(start, start + step.seq_len)]
    #             # codon = step.state.decode(seq[start : start + step.seq_len])[0]
    #             codon = Codon(b"XXX", self._alphabet)
    #             step.state.decode(subseq, codon)
    #             codons.append(codon)

    #             name = step.state.name
    #             abc = step.state.alphabet
    #             cstate = CodonState(name, abc, {codon: log(1.0)})
    #             # npath.append_codon_step(cstate, 3)
    #             steps.append((cstate, 3))

    #         start += step.seq_len

    #     return CodonFragment(codons, CodonPath(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


# class FrameSearchResult(SearchResult):
#     def __init__(
#         self,
#         alphabet: CBaseAlphabet,
#         loglik: float,
#         sequence: SequenceABC,
#         path: FramePath,
#     ):
#         self._loglik = loglik
#         self._fragments: List[FrameFragment] = []
#         self._intervals: List[Interval] = []

#         steps = list(path)
#         for fragi, stepi, homologous in self._create_fragments(path):
#             substeps = steps[stepi.start : stepi.stop]
#             fragment_path = FramePath([(s.state, s.seq_len) for s in substeps])
#             seq = sequence[fragi]
#             frag = FrameFragment(alphabet, seq, fragment_path, homologous)
#             self._fragments.append(frag)
#             self._intervals.append(fragi)

#     @property
#     def fragments(self) -> Sequence[FrameFragment]:
#         return self._fragments

#     @property
#     def intervals(self) -> Sequence[Interval]:
#         return self._intervals

#     @property
#     def loglikelihood(self) -> float:
#         return self._loglik

#     # def decode(self) -> CodonSearchResult:
#     #     fragments: List[CodonFragment] = []
#     #     intervals: List[Interval] = []

#     #     start = end = 0
#     #     for i, frag in enumerate(self._fragments):

#     #         codon_frag = frag.decode()
#     #         end += len(codon_frag.sequence)

#     #         fragments.append(codon_frag)
#     #         intervals.append(Interval(start, end))

#     #         start = end

#     #     return CodonSearchResult(self.score, fragments, intervals)
