import os
import pathlib
from typing import IO, List, Union

import click
from click.utils import LazyFile

from fasta_reader import FASTAItem, FASTAWriter, open_fasta
from hmmer_reader import HMMERParser, open_hmmer
from nmm import GeneticCode
from nmm.alphabet import CanonicalAminoAlphabet
from nmm.sequence import Sequence

from .._fasta import infer_fasta_alphabet
from .._gff import GFFItem, GFFWriter
from .._result import SearchResult
from ..frame import FrameProfile, create_profile


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
    target_alphabet = infer_fasta_alphabet(fasta)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    target.seek(0)

    gcode = GeneticCode(target_alphabet, CanonicalAminoAlphabet())
    with open_fasta(target) as fasta:
        targets = list(fasta)

    s = Scan(output_writer, codon_writer, amino_writer, gcode, epsilon)

    for hmmprof in open_hmmer(profile):
        s.scan_profile(hmmprof, targets)

    finalize_stream(output)
    finalize_stream(ocodon)
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


class Scan:
    def __init__(
        self,
        output_writer: OutputWriter,
        codon_writer: FASTAWriter,
        amino_writer: FASTAWriter,
        genetic_code: GeneticCode,
        epsilon: float,
    ):
        self._output_writer = output_writer
        self._codon_writer = codon_writer
        self._amino_writer = amino_writer
        self._genetic_code = genetic_code
        self._epsilon = epsilon

    def scan_profile(self, profile: HMMERParser, targets: List[FASTAItem]):

        self._output_writer.profile = dict(profile.metadata)["ACC"]
        prof = create_profile(
            profile, self._genetic_code.base_alphabet, epsilon=self._epsilon
        )

        show_header("Profile")
        show_profile(profile)

        show_header("Targets")
        for target in targets:
            self._scan_target(prof, target)
            print()

    def _scan_target(self, profile: FrameProfile, target: FASTAItem):

        print(">" + target.defline)
        print(sequence_summary(target.sequence))

        seq = Sequence.create(target.sequence.encode(), profile.alphabet)
        result = profile.search(seq)
        seqid = f"{target.defline.split()[0]}"

        show_search_result(result)

        for i, frag in zip(result.intervals, result.fragments):
            if not frag.homologous:
                continue

            start = i.start
            stop = i.stop
            item_id = self._output_writer.write_item(seqid, start, stop)

            codon_result = frag.decode()
            self._codon_writer.write_item(item_id, str(codon_result.sequence))

            amino_result = codon_result.decode(self._genetic_code)
            self._amino_writer.write_item(item_id, str(amino_result.sequence))


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
    print()


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


def show_header(title: str):
    print(title)
    print("=" * len(title))
    print()


def sequence_summary(sequence: str):
    max_nchars = 79
    if len(sequence) <= max_nchars:
        return sequence

    middle = " ... "

    begin_nchars = (max_nchars - len(middle)) // 2
    end_nchars = begin_nchars + (max_nchars - len(middle)) % 2

    return sequence[:begin_nchars] + middle + sequence[-end_nchars:]
