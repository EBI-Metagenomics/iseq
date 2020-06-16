from abc import ABC, abstractmethod
from typing import List, NamedTuple, Tuple

from click.utils import LazyFile
from fasta_reader import FASTAItem
from hmmer_reader import HMMERParser
from imm import Interval, Sequence
from nmm import DNAAlphabet, NTTranslator, NullTranslator, RNAAlphabet, Translator

from iseq.profile import Profile
from iseq.protein import ProteinFragment
from iseq.result import SearchResult

from .debug_writer import DebugWriter
from .output_writer import OutputWriter

IntFrag = NamedTuple("IntFrag", [("interval", Interval), ("fragment", ProteinFragment)])


class Scanner(ABC):
    def __init__(
        self, output_writer: OutputWriter, debug_writer: DebugWriter, stdout,
    ):
        self._output_writer = output_writer
        self._debug_writer = debug_writer
        self._stdout = stdout

    def finalize_stream(self, name: str, stream: LazyFile):
        if stream.name != "-":
            self._stdout.write(f"Writing {name} to <{stream.name}> file.\n")

        stream.close_intelligently()

    @abstractmethod
    def process_profile(
        self, profile_parser: HMMERParser, targets: List[FASTAItem], window_length: int
    ):
        del profile_parser
        del targets
        del window_length
        raise NotImplementedError()

    def show_profile_parser(self, profile_parser: HMMERParser):
        self._show_header("Profile")
        self._show_profile(profile_parser)

    def _show_header(self, title: str):
        self._stdout.write(title + "\n")
        self._stdout.write("=" * len(title) + "\n")
        self._stdout.write("\n")

    def _show_profile(self, hmmprof: HMMERParser):
        name = dict(hmmprof.metadata)["NAME"]
        acc = dict(hmmprof.metadata).get("ACC", "-")

        self._stdout.write(f"Header       {hmmprof.header}\n")
        self._stdout.write(f"Alphabet     {hmmprof.alphabet}\n")
        self._stdout.write(f"Model length {hmmprof.M}\n")
        self._stdout.write(f"Name         {name}\n")
        self._stdout.write(f"Accession    {acc}\n")
        self._stdout.write("\n")

    def _show_search_result(self, result: SearchResult, window: Interval):

        self._stdout.write("\n")

        start = window.start
        stop = window.stop
        n = sum(frag.homologous for frag in result.fragments)
        msg = f"Found {n} homologous fragment(s) within the range [{start+1}, {stop}]."
        self._stdout.write(msg + "\n")

        j = 0
        for interval, frag in zip(result.intervals, result.fragments):
            if not frag.homologous:
                continue

            start = window.start + interval.start
            stop = window.start + interval.stop
            msg = f"Fragment={j + 1}; Position=[{start + 1}, {stop}]\n"
            self._stdout.write(msg)
            states = []
            matches = []
            for frag_step in iter(frag):
                states.append(frag_step.step.state.name.decode())
                matches.append(str(frag_step.sequence))

            self._stdout.write("\t".join(states) + "\n")
            self._stdout.write("\t".join(matches) + "\n")
            j += 1

    def _scan_targets(self, profile: Profile, targets: List[FASTAItem]):
        self._show_header("Targets")
        if isinstance(profile.alphabet, (DNAAlphabet, RNAAlphabet)):
            translator = NTTranslator()
        else:
            translator = NullTranslator()

        for target in targets:
            self._scan_target(profile, target, translator)
            self._stdout.write("\n")

    def _scan_target(self, profile: Profile, target: FASTAItem, translator: Translator):

        self._stdout.write(">" + target.defline + "\n")
        self._stdout.write(sequence_summary(target.sequence) + "\n")

        tmp = translator.translate(target.sequence.encode(), profile.alphabet)
        seq = Sequence.create(tmp, profile.alphabet)
        search_results = profile.search(seq)
        seqid = f"{target.defline.split()[0]}"

        waiting: List[IntFrag] = []

        windows = search_results.windows
        results = search_results.results
        for win_num, (window, result) in enumerate(zip(windows, results)):

            self._debug_writer.write_row(
                seqid,
                win_num,
                window.start + 1,
                window.stop,
                result.alt_viterbi_score,
                result.null_viterbi_score,
            )

            self._show_search_result(result, window)
            candidates: List[IntFrag] = []

            for i, frag in zip(result.intervals, result.fragments):
                if not frag.homologous:
                    continue

                interval = Interval(window.start + i.start, window.start + i.stop)
                candidates.append(IntFrag(interval, frag))

            ready, waiting = intersect_fragments(waiting, candidates)

            self._write_fragments(profile, target, ready)

        self._write_fragments(profile, target, waiting)

    @abstractmethod
    def _write_fragments(
        self, profile: Profile, target: FASTAItem, ifragments: List[IntFrag]
    ):
        del profile
        del target
        del ifragments
        raise NotImplementedError()


def sequence_summary(sequence: str):
    max_nchars = 79
    if len(sequence) <= max_nchars:
        return sequence

    middle = " ... "

    begin_nchars = (max_nchars - len(middle)) // 2
    end_nchars = begin_nchars + (max_nchars - len(middle)) % 2

    return sequence[:begin_nchars] + middle + sequence[-end_nchars:]


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
