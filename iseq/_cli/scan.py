import os
import pathlib
from time import time
from typing import IO, List, NamedTuple, Union, Tuple, Iterator

import click
from click.utils import LazyFile

from fasta_reader import FASTAItem, FASTAWriter, open_fasta
from hmmer_reader import HMMERParser, open_hmmer
from nmm import GeneticCode, Interval
from nmm.alphabet import CanonicalAminoAlphabet
from nmm.sequence import Sequence

from .._fasta import infer_fasta_alphabet
from .._gff import GFFItem, GFFWriter
from .._result import SearchResult, SearchResults
from ..frame import FrameProfile, create_profile, FrameFragment


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option(
    "--epsilon", type=float, default=1e-2, help="Indel probability. Defaults to 1e-2."
)
@click.option(
    "--output",
    type=click.File("w"),
    help="Save results to OUTPUT (GFF format).",
    default=os.devnull,
)
@click.option(
    "--ocodon",
    type=click.File("w"),
    help="Save codon sequences to OCODON (FASTA format).",
    default=os.devnull,
)
@click.option(
    "--oamino",
    type=click.File("w"),
    help="Save amino acid sequences to OAMINO (FASTA format).",
    default=os.devnull,
)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
@click.option(
    "--window",
    type=int,
    help="Window length. Defaults to zero, which means no window.",
    default=0,
)
def scan(profile, target, epsilon: float, output, ocodon, oamino, quiet, window: int):
    """
    Search nucleotide sequence(s) against a protein profile database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """

    output_writer = OutputWriter(output, epsilon, window)
    codon_writer = FASTAWriter(ocodon)
    amino_writer = FASTAWriter(oamino)

    fasta = open_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    target.seek(0)

    gcode = GeneticCode(target_alphabet, CanonicalAminoAlphabet())
    with open_fasta(target) as fasta:
        targets = list(fasta)

    if quiet:
        stdout = click.open_file(os.devnull, "a")
    else:
        stdout = click.get_text_stream("stdout")

    scanner = Scanner(
        output_writer, codon_writer, amino_writer, gcode, epsilon, window, stdout
    )

    for i, hmmprof in enumerate(open_hmmer(profile)):
        # start = time()
        M = hmmprof.M
        scanner.process_profile(M, hmmprof, targets)
        assert len(targets) == 1
        # print("{},{},{:.12f}".format(M, len(targets[0].sequence), time() - start))
        if i == 176:
            break

    scanner.finalize_stream(output)
    scanner.finalize_stream(ocodon)
    scanner.finalize_stream(oamino)


class OutputWriter:
    def __init__(
        self, file: Union[str, pathlib.Path, IO[str]], epsilon: float, window: int
    ):
        self._gff = GFFWriter(file)
        self._profile = "NOTSET"
        self._epsilon = epsilon
        self._window = window
        self._item_idx = 1

    @property
    def profile(self) -> str:
        return self._profile

    @profile.setter
    def profile(self, profile: str):
        self._profile = profile

    def write_item(self, seqid: str, start: int, end: int):
        item_id = f"item{self._item_idx}"
        att = f"ID={item_id};Profile={self._profile};Epsilon={self._epsilon};Window={self._window}"
        item = GFFItem(seqid, "nmm", ".", start + 1, end, 0.0, "+", ".", att)
        self._gff.write_item(item)
        self._item_idx += 1
        return item_id

    def close(self):
        """
        Close the associated stream.
        """
        self._gff.close()


IntFrag = NamedTuple("IntFrag", [("interval", Interval), ("fragment", FrameFragment)])


class Scanner:
    def __init__(
        self,
        output_writer: OutputWriter,
        codon_writer: FASTAWriter,
        amino_writer: FASTAWriter,
        genetic_code: GeneticCode,
        epsilon: float,
        window_length: int,
        stdout,
    ):
        self._output_writer = output_writer
        self._codon_writer = codon_writer
        self._amino_writer = amino_writer
        self._genetic_code = genetic_code
        self._epsilon = epsilon
        self._window_length = window_length
        self._stdout = stdout

    def process_profile(self, M, profile_parser: HMMERParser, targets: List[FASTAItem]):

        self._output_writer.profile = dict(profile_parser.metadata)["ACC"]
        base_alphabet = self._genetic_code.base_alphabet
        start = time()
        prof = create_profile(profile_parser, base_alphabet, self._epsilon)
        elapsed_create = time() - start

        self._show_header("Profile")
        self._show_profile(profile_parser)

        self._show_header("Targets")
        start = time()
        for target in targets:
            self._scan_target(prof, target)
            self._stdout.write("\n")
        elapsed_scan = time() - start

        print(f"{M},{elapsed_create:.8f},{elapsed_scan:.8f}")

    def finalize_stream(self, stream: LazyFile):
        if stream.name != "-":
            self._stdout.write(f"Writing to <{stream.name}> file.\n")

        stream.close_intelligently()

    def _scan_target(self, profile: FrameProfile, target: FASTAItem):

        self._stdout.write(">" + target.defline + "\n")
        self._stdout.write(sequence_summary(target.sequence) + "\n")

        seq = Sequence.create(target.sequence.encode(), profile.alphabet)
        search_results = profile.search(seq, self._window_length)
        seqid = f"{target.defline.split()[0]}"

        waiting: List[IntFrag] = []

        for window, result in zip(search_results.windows, search_results.results):

            self._show_search_result(result, window)
            candidates: List[IntFrag] = []

            for i, frag in zip(result.intervals, result.fragments):
                if not frag.homologous:
                    continue

                interval = Interval(window.start + i.start, window.start + i.stop)
                candidates.append(IntFrag(interval, frag))

            ready, waiting = intersect_fragments(waiting, candidates)

            self._write_fragments(seqid, ready)

        self._write_fragments(seqid, waiting)

    def _write_fragments(self, seqid: str, ifragments: List[IntFrag]):
        for ifrag in ifragments:
            start = ifrag.interval.start
            stop = ifrag.interval.stop
            item_id = self._output_writer.write_item(seqid, start, stop)

            codon_result = ifrag.fragment.decode()
            self._codon_writer.write_item(item_id, str(codon_result.sequence))

            amino_result = codon_result.decode(self._genetic_code)
            self._amino_writer.write_item(item_id, str(amino_result.sequence))

    def _show_profile(self, hmmprof: HMMERParser):
        name = dict(hmmprof.metadata)["NAME"]
        acc = dict(hmmprof.metadata)["ACC"]

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

    def _show_header(self, title: str):
        self._stdout.write(title + "\n")
        self._stdout.write("=" * len(title) + "\n")
        self._stdout.write("\n")


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
