from __future__ import annotations

from dataclasses import astuple, dataclass
from typing import Callable, Generic, Iterable, List, NamedTuple, Tuple, TypeVar

from imm import Alphabet, Interval, Path, SequenceABC, State, Step

from .fragment import Fragment

__all__ = ["SearchResults", "SearchResult", "create_fragment_type"]

A = TypeVar("A", bound=Alphabet)
S = TypeVar("S", bound=State)

create_fragment_type = Callable[[SequenceABC[A], S, bool], Fragment[A, S]]


@dataclass
class IFragment(Generic[A, S]):
    interval: Interval
    fragment: Fragment[A, S]

    def __iter__(self):
        yield from astuple(self)


DebugRow = NamedTuple(
    "DebugRow",
    [
        ("window_number", int),
        ("window", Interval),
        ("alt_viterbi_score", float),
        ("null_viterbi_score", float),
    ],
)


class SearchResults(Generic[A, S]):
    def __init__(
        self,
        sequence: SequenceABC[A],
        create_fragment: create_fragment_type,
    ):
        self._sequence = sequence
        self._create_fragment = create_fragment
        self._results: List[SearchResult[A, S]] = []
        self._windows: List[Interval] = []

    def append(
        self,
        loglik: float,
        window: Interval,
        path: S,
        alt_viterbi_score: float,
        null_viterbi_score: float,
    ):
        subseq = self._sequence[window.start : window.stop]
        r = SearchResult[A, S](
            loglik,
            subseq,
            path,
            self._create_fragment,
            alt_viterbi_score,
            null_viterbi_score,
        )
        self._results.append(r)
        self._windows.append(window)

    @property
    def results(self) -> List[SearchResult[A, S]]:
        return self._results

    @property
    def windows(self) -> List[Interval]:
        return self._windows

    @property
    def length(self) -> int:
        return len(self._results)

    def debug_table(self) -> List[DebugRow]:
        windows = self.windows
        results = self.results
        rows: List[DebugRow] = []
        for win_num, (window, result) in enumerate(zip(windows, results)):

            rows.append(
                DebugRow(
                    win_num,
                    window,
                    result.alt_viterbi_score,
                    result.null_viterbi_score,
                )
            )

        return rows

    def ifragments(self) -> List[IFragment[A, S]]:
        waiting: List[IFragment[A, S]] = []

        windows = self.windows
        results = self.results
        ifragments: List[IFragment[A, S]] = []
        for window, result in zip(windows, results):

            candidates: List[IFragment[A, S]] = []

            for i, frag in zip(result.intervals, result.fragments):
                if not frag.homologous:
                    continue

                interval = Interval(window.start + i.start, window.start + i.stop)
                candidates.append(IFragment(interval, frag))

            ready, waiting = intersect_ifragments(waiting, candidates)
            ifragments.extend(ready)

        ifragments.extend(waiting)

        return ifragments

    def __str__(self) -> str:
        return f"{str(self._results)}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class SearchResult(Generic[A, S]):
    def __init__(
        self,
        loglik: float,
        sequence: SequenceABC[A],
        path: Path[S],
        create_fragment: create_fragment_type,
        alt_viterbi_score: float,
        null_viterbi_score: float,
    ):
        self._loglik = loglik
        self._fragments: List[Fragment[A, S]] = []
        self._intervals: List[Interval] = []
        self._alt_viterbi_score = alt_viterbi_score
        self._null_viterbi_score = null_viterbi_score

        steps = list(path)
        for fragi, stepi, homologous in create_fragments(path):
            substeps = steps[stepi.start : stepi.stop]
            new_steps = [Step.create(s.state, s.seq_len) for s in substeps]
            new_path = Path.create(new_steps)
            seq = sequence[fragi]
            frag = create_fragment(seq, new_path, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def alt_viterbi_score(self) -> float:
        return self._alt_viterbi_score

    @property
    def null_viterbi_score(self) -> float:
        return self._null_viterbi_score

    @property
    def fragments(self) -> List[Fragment[A, S]]:
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


def create_fragments(path: Path) -> Iterable[Tuple[Interval, Interval, bool]]:

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


def intersect_ifragments(
    waiting: List[IFragment[A, S]], candidates: List[IFragment[A, S]]
) -> Tuple[List[IFragment[A, S]], List[IFragment[A, S]]]:

    ready: List[IFragment[A, S]] = []
    new_waiting: List[IFragment[A, S]] = []

    i = 0
    j = 0

    curr_stop = 0
    while i < len(waiting) and j < len(candidates):

        if waiting[i].interval.start < candidates[j].interval.start:
            ready.append(waiting[i])
            curr_stop = waiting[i].interval.stop
            i += 1
        elif waiting[i].interval.start == candidates[j].interval.start:
            if waiting[i].interval.stop >= candidates[j].interval.stop:
                ready.append(waiting[i])
                curr_stop = waiting[i].interval.stop
            else:
                new_waiting.append(candidates[j])
                curr_stop = candidates[j].interval.stop
            i += 1
            j += 1
        else:
            new_waiting.append(candidates[j])
            curr_stop = candidates[j].interval.stop
            j += 1

        while i < len(waiting) and waiting[i].interval.stop <= curr_stop:
            i += 1

        while j < len(candidates) and candidates[j].interval.stop <= curr_stop:
            j += 1

    while i < len(waiting):
        ready.append(waiting[i])
        i += 1

    while j < len(candidates):
        new_waiting.append(candidates[j])
        j += 1

    return ready, new_waiting
