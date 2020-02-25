from __future__ import annotations
from typing import Callable, Generic, Iterable, List, Tuple

from nmm import Interval
from nmm.path import Path, Step
from nmm.sequence import SequenceABC

from ._fragment import Fragment
from ._typing import MutablePath, TAlphabet, TState

create_fragment_type = Callable[
    [SequenceABC[TAlphabet], MutablePath[TState], bool], Fragment[TAlphabet, TState]
]

__all__ = ["SearchResults", "SearchResult"]


class SearchResults(Generic[TAlphabet, TState]):
    def __init__(
        self, sequence: SequenceABC[TAlphabet], create_fragment: create_fragment_type,
    ):
        self._sequence = sequence
        self._create_fragment = create_fragment
        self._results: List[SearchResult[TAlphabet, TState]] = []
        self._windows: List[Interval] = []

    def append(
        self, loglik: float, window: Interval, path: MutablePath[TState],
    ):
        subseq = self._sequence[window.start : window.stop]
        r = SearchResult[TAlphabet, TState](loglik, subseq, path, self._create_fragment)
        self._results.append(r)
        self._windows.append(window)

    @property
    def results(self):
        return self._results

    @property
    def windows(self):
        return self._windows

    @property
    def length(self) -> int:
        return len(self._results)

    def __str__(self) -> str:
        return f"{str(self._results)}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class SearchResult(Generic[TAlphabet, TState]):
    def __init__(
        self,
        loglik: float,
        sequence: SequenceABC[TAlphabet],
        path: MutablePath[TState],
        create_fragment: create_fragment_type,
    ):
        self._loglik = loglik
        self._fragments: List[Fragment[TAlphabet, TState]] = []
        self._intervals: List[Interval] = []

        steps = list(path)
        for fragi, stepi, homologous in _create_fragments(path):
            substeps = steps[stepi.start : stepi.stop]
            new_steps = [Step.create(s.state, s.seq_len) for s in substeps]
            new_path = Path.create(new_steps)
            seq = sequence[fragi]
            frag = create_fragment(seq, new_path, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def fragments(self) -> List[Fragment[TAlphabet, TState]]:
        return self._fragments

    @property
    def intervals(self) -> List[Interval]:
        return self._intervals

    @property
    def loglikelihood(self) -> float:
        return self._loglik

    def __str__(self) -> str:
        return f"{str(self.loglikelihood)},{str(self._fragments)}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


def _create_fragments(path: Path) -> Iterable[Tuple[Interval, Interval, bool]]:

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
