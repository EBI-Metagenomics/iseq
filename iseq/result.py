from abc import ABC, abstractmethod
from typing import Sequence

from nmm import Interval
from nmm.path import CPath

from .fragment import Fragment


class SearchResult(ABC):
    def _create_fragments(self, path: CPath):

        frag_start = frag_stop = 0
        step_start = step_stop = 0
        homologous = False

        for step_stop, step in enumerate(path):

            change = not homologous and step.state.name.startswith(b"M")
            change = change or homologous and step.state.name.startswith(b"E")
            change = change or not homologous and step.state.name.startswith(b"T")

            if change:
                if frag_start < frag_stop:
                    fragi = Interval(frag_start, frag_stop)
                    stepi = Interval(step_start, step_stop)
                    yield (fragi, stepi, homologous)

                frag_start = frag_stop
                step_start = step_stop
                homologous = not homologous

            frag_stop += step.seq_len

    # @property
    # def symbols(self) -> bytes:
    #     return b"".join(frag.subsequence.symbols for frag in self.fragments)

    @property
    @abstractmethod
    def fragments(self) -> Sequence[Fragment]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def intervals(self) -> Sequence[Interval]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def loglikelihood(self) -> float:
        raise NotImplementedError()

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
