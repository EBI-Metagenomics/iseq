from typing import List, Sequence

from nmm import Interval, SequenceABC

from ..result import SearchResult
from .fragment import FrameFragment
from .path import FramePath


class FrameSearchResult(SearchResult):
    def __init__(self, loglik: float, sequence: SequenceABC, path: FramePath):
        self._loglik = loglik
        self._fragments: List[FrameFragment] = []
        self._intervals: List[Interval] = []

        steps = list(path)
        for fragi, stepi, homologous in self._create_fragments(path):
            substeps = steps[stepi.start : stepi.stop]
            fragment_path = FramePath([(s.state, s.seq_len) for s in substeps])
            seq = sequence.slice(fragi)
            frag = FrameFragment(seq, fragment_path, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def fragments(self) -> Sequence[FrameFragment]:
        return self._fragments

    @property
    def intervals(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def loglikelihood(self) -> float:
        return self._loglik

    # def decode(self) -> CodonSearchResult:
    #     fragments: List[CodonFragment] = []
    #     intervals: List[Interval] = []

    #     start = end = 0
    #     for i, frag in enumerate(self._fragments):

    #         codon_frag = frag.decode()
    #         end += len(codon_frag.sequence)

    #         fragments.append(codon_frag)
    #         intervals.append(Interval(start, end))

    #         start = end

    #     return CodonSearchResult(self.score, fragments, intervals)
