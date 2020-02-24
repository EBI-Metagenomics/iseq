import os
import pathlib
from typing import Any, IO, Union

import click
from click.utils import LazyFile

from fasta_reader import FASTAParser, FASTAWriter, open_fasta
from hmmer_reader import HMMERParser, open_hmmer
from nmm import GeneticCode
from nmm.alphabet import CanonicalAminoAlphabet, DNAAlphabet
from nmm.sequence import Sequence

from .._gff import GFFItem, GFFWriter
from .._result import SearchResult
from ..frame import FrameProfile, create_profile
from .._fasta import infer_fasta_alphabet


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
def scan(profile, target, epsilon: float, output, ocodon, oamino):
    """
    Search nucleotide sequence(s) against a protein profile database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """

    output_writer = OutputWriter(output, epsilon)
    codon_writer = FASTAWriter(ocodon)
    amino_writer = FASTAWriter(oamino)

    fasta = open_fasta(target)
    alphabet = infer_fasta_alphabet(fasta)
    if alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    target.seek(0)

    gcode = GeneticCode(alphabet, CanonicalAminoAlphabet())
    with open_fasta(target) as fasta:
        targets = list(fasta)

    for hmmprof in open_hmmer(profile):
        output_writer.profile = dict(hmmprof.metadata)["ACC"]
        prof = create_profile(hmmprof, gcode.base_alphabet, epsilon=epsilon)

        show_header1("Profile")
        print()
        show_profile(hmmprof)
        print()

        show_header1("Targets")

        process_sequence(
            prof, targets, output_writer, codon_writer, amino_writer, gcode
        )

        print()

    if output is not None:
        finalize_stream(output)

    if codon_writer is not None:
        finalize_stream(ocodon)

    if amino_writer is not None:
        finalize_stream(oamino)


class OutputWriter:
    def __init__(self, file: Union[str, pathlib.Path, IO[str]], epsilon: float):
        self._gff = GFFWriter(file)
        self._profile = "NOTSET"
        self._epsilon = epsilon
        self._item_idx = 1

    @property
    def profile(self) -> str:
        return self._profile

    @profile.setter
    def profile(self, profile: str):
        self._profile = profile

    def write_item(self, seqid: str, start: int, end: int):
        item_id = f"item{self._item_idx}"
        att = f"ID={item_id};Profile={self._profile};Epsilon={self._epsilon}"
        item = GFFItem(seqid, "nmm", ".", start + 1, end, 0.0, "+", ".", att)
        self._gff.write_item(item)
        self._item_idx += 1
        return item_id

    def close(self):
        """
        Close the associated stream.
        """
        self._gff.close()


def process_sequence(
    prof: FrameProfile,
    fasta: FASTAParser,
    record: OutputWriter,
    codon_writer: Union[FASTAWriter, None],
    amino_writer: Union[FASTAWriter, None],
    gcode,
):

    for tgt in fasta:
        print()
        print(">" + tgt.defline)
        print(sequence_summary(tgt.sequence))

        # seq = tgt.sequence.encode().replace(b"T", b"U")
        seq = Sequence[DNAAlphabet].create(tgt.sequence.encode(), prof.alphabet)
        frame_result = prof.search(seq)
        # codon_result = frame_result.decode()
        seqid = f"{tgt.defline.split()[0]}"

        show_search_result(frame_result)

        for i, frag in enumerate(frame_result.fragments):
            if not frag.homologous:
                continue

            start = frame_result.intervals[i].start
            stop = frame_result.intervals[i].stop
            item_id = record.write_item(seqid, start, stop)
            codon_result = frag.decode()

            if codon_writer is not None:
                codon_writer.write_item(item_id, str(codon_result.sequence))

            if amino_writer is not None:
                amino_result = codon_result.decode(gcode)
                amino_writer.write_item(item_id, str(amino_result.sequence))


def finalize_stream(stream: LazyFile):
    if stream.name != "-":
        print(f"Writing to <{stream.name}> file.")

    stream.close_intelligently()


def show_profile(hmmprof: HMMERParser):
    name = dict(hmmprof.metadata)["NAME"]
    acc = dict(hmmprof.metadata)["ACC"]

    print(f"Header       {hmmprof.header}")
    print(f"Alphabet     {hmmprof.alphabet}")
    print(f"Model length {hmmprof.M}")
    print(f"Name         {name}")
    print(f"Accession    {acc}")


def show_search_result(result: SearchResult):
    frags = result.fragments
    nhomo = sum(frag.homologous for frag in result.fragments)

    print()
    print(f"Found {nhomo} homologous fragments ({len(frags)} in total).")

    for i, frag in enumerate(frags):
        if not frag.homologous:
            continue
        start = result.intervals[i].start
        stop = result.intervals[i].stop
        print(f"Homologous fragment={i}; Position=[{start + 1}, {stop}]")
        states = []
        matches = []
        for frag_step in iter(frag):
            states.append(frag_step.step.state.name.decode())
            matches.append(str(frag_step.sequence))

        print("\t".join(states))
        print("\t".join(matches))


def show_header1(title: str):
    print(title)
    print("=" * len(title))


def show_header2(title: str):
    print(title)
    print("-" * len(title))


def sequence_summary(sequence: str):
    max_nchars = 79
    if len(sequence) <= max_nchars:
        return sequence

    middle = " ... "

    begin_nchars = (max_nchars - len(middle)) // 2
    end_nchars = begin_nchars + (max_nchars - len(middle)) % 2

    return sequence[:begin_nchars] + middle + sequence[-end_nchars:]
