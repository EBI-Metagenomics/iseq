import os
from time import time

import click

from click.utils import LazyFile
from typing import List, NamedTuple, Tuple
from fasta_reader import FASTAItem, FASTAWriter, open_fasta
from imm import Sequence, Interval
from nmm import DNAAlphabet, Input, GeneticCode, CanonicalAminoAlphabet
from .._result import SearchResult
from .. import wrap
from ..frame import FrameFragment
from ..frame._typing import FrameAltModel, FrameNullModel
from ..frame._profile import FrameProfile
from .scan import _infer_target_alphabet, OutputWriter

IntFrag = NamedTuple("IntFrag", [("interval", Interval), ("fragment", FrameFragment)])


@click.command()
@click.argument("profile", type=str)
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
@click.option(
    "--hmmer3-compat/--no-hmmer3-compat",
    help="Enable full HMMER3 compatibility.",
    default=False,
)
def bscan(
    profile,
    target,
    epsilon: float,
    output,
    ocodon,
    oamino,
    quiet,
    window: int,
    hmmer3_compat: bool,
):
    """
    Press.

    TODO: fix this doc.
    Search nucleotide sequence(s) against a protein profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """
    owriter = OutputWriter(output, epsilon, window)
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)

    if quiet:
        stdout = click.open_file(os.devnull, "a")
    else:
        stdout = click.get_text_stream("stdout")

    alt_filepath = (profile + ".alt").encode()
    null_filepath = (profile + ".null").encode()
    meta_filepath = (profile + ".meta").encode()

    target_abc = _infer_target_alphabet(target)
    gcode = GeneticCode(target_abc, CanonicalAminoAlphabet())

    with open_fasta(target) as fasta:
        targets = list(fasta)

    ELAPSED = {"wrap": 0.0, "scan": 0.0}

    with Input.create(alt_filepath) as afile:
        with Input.create(null_filepath) as nfile:
            with open(meta_filepath, "r") as mfile:
                for alt, null, acc in zip(afile, nfile, mfile):
                    owriter.profile = acc.strip()
                    start = time()
                    special_node = wrap.special_node(alt.hmm)
                    core_nodes = wrap.core_nodes(alt.hmm)
                    alt_model = FrameAltModel.create2(
                        special_node, core_nodes, alt.hmm, alt.dp
                    )
                    # print(alt_model)

                    null_model = FrameNullModel.create2(null.hmm)
                    # print(null_model)

                    abc = alt_model.hmm.alphabet
                    prof = FrameProfile.create2(
                        abc, null_model, alt_model, hmmer3_compat
                    )
                    ELAPSED["wrap"] += time() - start

                    for tgt in targets:
                        stdout.write(">" + tgt.defline + "\n")
                        stdout.write(sequence_summary(tgt.sequence) + "\n")

                        seq = Sequence.create(tgt.sequence.encode(), prof.alphabet)
                        start = time()
                        search_results = prof.search(seq, window)
                        ELAPSED["scan"] += time() - start
                        seqid = f"{tgt.defline.split()[0]}"

                        waiting: List[IntFrag] = []

                        for win, result in zip(
                            search_results.windows, search_results.results
                        ):

                            _show_search_result(stdout, result, win)
                            candidates: List[IntFrag] = []

                            for i, frag in zip(result.intervals, result.fragments):
                                if not frag.homologous:
                                    continue

                                interval = Interval(
                                    win.start + i.start, win.start + i.stop
                                )
                                candidates.append(IntFrag(interval, frag))

                            ready, waiting = intersect_fragments(waiting, candidates)

                            _write_fragments(
                                gcode, owriter, cwriter, awriter, seqid, ready
                            )

                        _write_fragments(
                            gcode, owriter, cwriter, awriter, seqid, waiting
                        )

                        stdout.write("\n")

    finalize_stream(stdout, "output", output)
    finalize_stream(stdout, "ocodon", ocodon)
    finalize_stream(stdout, "oamino", oamino)
    print(ELAPSED)
    from ..frame._profile import _ELAPSED

    print(_ELAPSED)


def sequence_summary(sequence: str):
    max_nchars = 79
    if len(sequence) <= max_nchars:
        return sequence

    middle = " ... "

    begin_nchars = (max_nchars - len(middle)) // 2
    end_nchars = begin_nchars + (max_nchars - len(middle)) % 2

    return sequence[:begin_nchars] + middle + sequence[-end_nchars:]


def _show_search_result(stdout, result: SearchResult, window: Interval):

    stdout.write("\n")

    start = window.start
    stop = window.stop
    n = sum(frag.homologous for frag in result.fragments)
    msg = f"Found {n} homologous fragment(s) within the range [{start+1}, {stop}]."
    stdout.write(msg + "\n")

    j = 0
    for interval, frag in zip(result.intervals, result.fragments):
        if not frag.homologous:
            continue

        start = window.start + interval.start
        stop = window.start + interval.stop
        msg = f"Fragment={j + 1}; Position=[{start + 1}, {stop}]\n"
        stdout.write(msg)
        states = []
        matches = []
        for frag_step in iter(frag):
            states.append(frag_step.step.state.name.decode())
            matches.append(str(frag_step.sequence))

        stdout.write("\t".join(states) + "\n")
        stdout.write("\t".join(matches) + "\n")
        j += 1


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


def _write_fragments(
    genetic_code,
    output_writer,
    codon_writer,
    amino_writer,
    seqid: str,
    ifragments: List[IntFrag],
):
    for ifrag in ifragments:
        start = ifrag.interval.start
        stop = ifrag.interval.stop
        item_id = output_writer.write_item(seqid, start, stop)

        codon_result = ifrag.fragment.decode()
        codon_writer.write_item(item_id, str(codon_result.sequence))

        amino_result = codon_result.decode(genetic_code)
        amino_writer.write_item(item_id, str(amino_result.sequence))


def finalize_stream(stdout, name: str, stream: LazyFile):
    if stream.name != "-":
        stdout.write(f"Writing {name} to <{stream.name}> file.\n")

    stream.close_intelligently()
