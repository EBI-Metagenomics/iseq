import os
from typing import IO, List, NamedTuple, Tuple

import click
from click.utils import LazyFile
from fasta_reader import FASTAWriter, open_fasta
from imm import Interval, Sequence
from nmm import CanonicalAminoAlphabet, GeneticCode, Input

from iseq import wrap
from iseq.alphabet import infer_fasta_alphabet
from iseq.hmmsearch import HMMSearch
from iseq.profile import ProfileID
from iseq.protein import ProteinFragment, ProteinProfile
from iseq.protein.typing import ProteinAltModel, ProteinNullModel
from iseq.result import SearchResult
from pathlib import Path

from .debug_writer import DebugWriter
from .output_writer import OutputWriter
from .pscan import update_gff_file

IntFrag = NamedTuple("IntFrag", [("interval", Interval), ("fragment", ProteinFragment)])


@click.command()
@click.argument(
    "profile",
    type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True),
)
@click.argument("target", type=click.File("r"))
@click.option(
    "--output",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save results to OUTPUT (GFF format).",
    default="output.gff",
)
@click.option(
    "--ocodon",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save codon sequences to OCODON (FASTA format).",
    default="ocodon.fasta",
)
@click.option(
    "--oamino",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save amino acid sequences to OAMINO (FASTA format).",
    default="oamino.fasta",
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
    "--odebug",
    type=click.File("w"),
    help="Save debug info into a tab-separated values file.",
    default=os.devnull,
)
@click.option(
    "--e-value/--no-e-value",
    help="Enable E-value computation. Defaults to True.",
    default=True,
)
def bscan(
    profile, target, output, ocodon, oamino, quiet, window: int, odebug, e_value: bool,
):
    """
    Binary scan.
    """

    owriter = OutputWriter(output)
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)
    dwriter = DebugWriter(odebug)

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
                    prof = ProteinProfile.create2(profid, abc, null_model, alt_model)
                    prof.window_length = window

                    for tgt in targets:
                        seq = prof.create_sequence(tgt.sequence.encode())
                        search_results = prof.search(seq)
                        ifragments = search_results.ifragments()
                        seqid = f"{tgt.defline.split()[0]}"

                        for ifrag in ifragments:
                            start = ifrag.interval.start
                            stop = ifrag.interval.stop
                            item_id = owriter.write_item(
                                seqid,
                                prof.profid,
                                start,
                                stop,
                                prof.window_length,
                                {"Epsilon": 0.01},
                            )
                            codon_frag = ifrag.fragment.decode()
                            cwriter.write_item(item_id, str(codon_frag.sequence))
                            amino_frag = codon_frag.decode(gcode)
                            awriter.write_item(item_id, str(amino_frag.sequence))

                        if odebug is not os.devnull:
                            for i in search_results.debug_table():
                                dwriter.write_row(seqid, i)

    owriter.close()
    cwriter.close()
    awriter.close()
    odebug.close_intelligently()

    if e_value:
        hmmsearch = HMMSearch()
        tbldata = hmmsearch.search(Path(profile), Path(oamino))
        update_gff_file(output, tbldata)


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
