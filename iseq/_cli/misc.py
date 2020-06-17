from typing import List, NamedTuple, Tuple
from iseq.protein import ProteinFragment
from imm import Interval

IntFrag = NamedTuple("IntFrag", [("interval", Interval), ("fragment", ProteinFragment)])


def consolidate(search_results):
    waiting: List[IntFrag] = []

    windows = search_results.windows
    results = search_results.results
    intfrags = []
    debug_list = []
    for win_num, (window, result) in enumerate(zip(windows, results)):

        debug_list.append(
            (
                win_num,
                window.start + 1,
                window.stop,
                result.alt_viterbi_score,
                result.null_viterbi_score,
            )
        )

        # self._show_search_result(result, window)
        candidates: List[IntFrag] = []

        for i, frag in zip(result.intervals, result.fragments):
            if not frag.homologous:
                continue

            interval = Interval(window.start + i.start, window.start + i.stop)
            candidates.append(IntFrag(interval, frag))

        ready, waiting = intersect_fragments(waiting, candidates)
        intfrags += ready
        # intervals += write_fragments(ready)

    # intervals += write_fragments(waiting)
    intfrags += waiting
    return intfrags, debug_list


def write_fragments(ifragments: List[IntFrag]):
    # seqid = f"{target.defline.split()[0]}"
    intervals = []
    for ifrag in ifragments:
        start = ifrag.interval.start
        stop = ifrag.interval.stop
        intervals.append(Interval(start, stop))
        # self._output_writer.write_item(seqid, start, stop, profile.window_length)
    return intervals


def intersect_fragments(
    waiting: List[IntFrag], candidates: List[IntFrag]
) -> Tuple[List[IntFrag], List[IntFrag]]:

    ready: List[IntFrag] = []
    new_waiting: List[IntFrag] = []

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
