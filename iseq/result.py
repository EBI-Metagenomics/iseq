from typing import Generic, Iterable, List, Sequence, Tuple, TypeVar

from nmm import Interval
from nmm.alphabet import CAlphabet
from nmm.path import CPath, Path, Step
from nmm.sequence import SequenceABC
from nmm.state import CState

from .fragment import Fragment

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)
TState = TypeVar("TState", bound=CState)


class SearchResult(Generic[TAlphabet, TState]):
    def __init__(
        self, loglik: float, sequence: SequenceABC[TAlphabet], path: Path[Step[TState]]
    ):
        self._loglik = loglik
        self._fragments: List[Fragment[TAlphabet, TState]] = []
        self._intervals: List[Interval] = []

        steps = list(path)
        for fragi, stepi, homologous in _create_fragments(path):
            substeps = steps[stepi.start : stepi.stop]
            fragment_path = Path([Step(s.state, s.seq_len) for s in substeps])
            seq = sequence[fragi]
            frag: Fragment[TAlphabet, TState] = Fragment(seq, fragment_path, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def fragments(self) -> Sequence[Fragment[TAlphabet, TState]]:
        return self._fragments

    @property
    def intervals(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def loglikelihood(self) -> float:
        return self._loglik

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


def _create_fragments(path: CPath) -> Iterable[Tuple[Interval, Interval, bool]]:

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
