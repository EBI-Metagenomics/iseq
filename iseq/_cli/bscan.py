import os
from typing import IO, List, NamedTuple, Tuple

import click
from click.utils import LazyFile
from fasta_reader import FASTAWriter, open_fasta
from imm import Interval, Sequence
from nmm import CanonicalAminoAlphabet, GeneticCode, Input

from iseq import wrap
from iseq.alphabet import infer_fasta_alphabet
from iseq.profile import ProfileID
from iseq.protein import ProteinFragment, ProteinProfile
from iseq.protein.typing import ProteinAltModel, ProteinNullModel
from iseq.result import SearchResult

from .output_writer import OutputWriter

IntFrag = NamedTuple("IntFrag", [("interval", Interval), ("fragment", ProteinFragment)])


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
    owriter = OutputWriter(output)
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

    with Input.create(alt_filepath) as afile:
        with Input.create(null_filepath) as nfile:
            with open(meta_filepath, "r") as mfile:
                for alt, null, line in zip(afile, nfile, mfile):
                    name, acc = line.strip().split("\t")
                    special_node = wrap.special_node(alt.hmm)
                    core_nodes = wrap.core_nodes(alt.hmm)
                    alt_model = ProteinAltModel.create2(
                        special_node, core_nodes, alt.hmm, alt.dp
                    )

                    null_model = ProteinNullModel.create2(null.hmm)

                    abc = alt_model.hmm.alphabet
                    profid = ProfileID(name, acc)
                    prof = ProteinProfile.create2(
                        profid, abc, null_model, alt_model, hmmer3_compat
                    )
                    prof.window_length = window

                    for tgt in targets:
                        # stdout.write(">" + tgt.defline + "\n")
                        # stdout.write(sequence_summary(tgt.sequence) + "\n")

                        seq = Sequence.create(tgt.sequence.encode(), prof.alphabet)
                        search_results = prof.search(seq)
                        seqid = f"{tgt.defline.split()[0]}"

                        waiting: List[IntFrag] = []

                        for win, result in zip(
                            search_results.windows, search_results.results
                        ):

                            # _show_search_result(stdout, result, win)
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
                                profid,
                                gcode,
                                owriter,
                                cwriter,
                                awriter,
                                seqid,
                                ready,
                                window,
                            )

                        _write_fragments(
                            profid,
                            gcode,
                            owriter,
                            cwriter,
                            awriter,
                            seqid,
                            waiting,
                            window,
                        )

                        # stdout.write("\n")

    finalize_stream(stdout, "output", output)
    finalize_stream(stdout, "ocodon", ocodon)
    finalize_stream(stdout, "oamino", oamino)


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
    profid,
    genetic_code,
    output_writer,
    codon_writer,
    amino_writer,
    seqid: str,
    ifragments: List[IntFrag],
    window_length: int,
):
    for ifrag in ifragments:
        start = ifrag.interval.start
        stop = ifrag.interval.stop
        item_id = output_writer.write_item(seqid, profid, start, stop, window_length)

        codon_result = ifrag.fragment.decode()
        codon_writer.write_item(item_id, str(codon_result.sequence))

        amino_result = codon_result.decode(genetic_code)
        amino_writer.write_item(item_id, str(amino_result.sequence))


def finalize_stream(stdout, name: str, stream: LazyFile):
    if stream.name != "-":
        stdout.write(f"Writing {name} to <{stream.name}> file.\n")

    stream.close_intelligently()


def _infer_target_alphabet(target: IO[str]):
    fasta = open_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    target.seek(0)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    return target_alphabet
