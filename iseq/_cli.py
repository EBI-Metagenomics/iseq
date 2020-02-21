from typing import Any, Union

import click
from click.utils import LazyFile

from fasta_reader import FASTAWriter, FASTAParser, open_fasta
from hmmer_reader import open_hmmer
from nmm import GeneticCode
from ._gff import GFFItem, GFFWriter
from ._result import SearchResult
from .frame import FrameProfile


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option("--epsilon", type=float, default=1e-2)
@click.option("--output", type=click.File("w"))
@click.option("--ocodon", type=click.File("w"))
@click.option("--oamino", type=click.File("w"))
def search(profile, target, epsilon: float, output, ocodon, oamino):
    """
    Search nucleotide sequences against a HMMER3 Protein profile.
    """
    from .frame import create_profile
    from nmm.alphabet import DNAAlphabet, CanonicalAminoAlphabet

    record_writer = RecordWriter(epsilon)

    codon_writer: Union[FASTAWriter, None] = None
    if ocodon is not None:
        codon_writer = FASTAWriter(ocodon)

    amino_writer: Union[FASTAWriter, None] = None
    if oamino is not None:
        amino_writer = FASTAWriter(oamino)

    gcode = GeneticCode(DNAAlphabet(), CanonicalAminoAlphabet())
    with open_fasta(target) as fasta:
        targets = list(fasta)

    for hmmprof in open_hmmer(profile):
        record_writer.profile = dict(hmmprof.metadata)["ACC"]
        prof = create_profile(hmmprof, gcode.base_alphabet, epsilon=epsilon)

        show_header1("Profile")
        print()
        print(hmmprof)
        print()

        show_header1("Targets")

        process_sequence(
            prof, targets, record_writer, codon_writer, amino_writer, gcode
        )

        print()

    if output is not None:
        record_writer.dump(output)
        finalize_stream(output)

    if codon_writer is not None:
        finalize_stream(ocodon)

    if amino_writer is not None:
        finalize_stream(oamino)


cli.add_command(search)


class RecordWriter:
    def __init__(self, epsilon: float):
        self._gff = GFFWriter()
        self._profile = "NOTSET"
        self._epsilon = epsilon
        self._item_idx = 1

    @property
    def profile(self):
        return self._profile

    @profile.setter
    def profile(self, profile: str):
        self._profile = profile

    # def add_items(self, result: SearchResult, seqid: str):
    #     for i, frag in enumerate(result.fragments):
    #         if not frag.homologous:
    #             continue

    #         start = result.intervals[i].start
    #         end = result.intervals[i].end

    #         ID = f"item{self._item_idx}"
    #         att = f"ID={ID};Profile={self._profile};Epsilon={self._epsilon}"
    #         item = GFFItem(seqid, "nmm", ".", start + 1, end, 0.0, "+", ".", att)
    #         self._gff.append(item)
    #         self._item_idx += 1

    def add_item(self, seqid: str, start: int, end: int):
        item_id = f"item{self._item_idx}"
        att = f"ID={item_id};Profile={self._profile};Epsilon={self._epsilon}"
        item = GFFItem(seqid, "nmm", ".", start + 1, end, 0.0, "+", ".", att)
        self._gff.append(item)
        self._item_idx += 1
        return item_id

    def dump(self, fp):
        self._gff.dump(fp)


class TargetWriter:
    def __init__(self, ocodon=Union[LazyFile, Any], oamino=Union[LazyFile, Any]):
        self._ocodon = ocodon
        self._oamino = oamino

    @property
    def ocodon(self):
        return self._ocodon

    @property
    def oamino(self):
        return self._oamino

    @property
    def has_ocodon(self):
        return self._ocodon is not None

    @property
    def has_oamino(self):
        return self._oamino is not None

    def add_ocodon_target(self, seqid: str, sequence: str):
        if self._ocodon is None:
            return

        self._ocodon.write(">" + seqid + "_codon\n")
        self._ocodon.write(sequence + "\n")

    def add_oamino_target(self, seqid: str, sequence: str):
        if self._oamino is None:
            return

        self._oamino.write(">" + seqid + "_amino\n")
        self._oamino.write(sequence + "\n")


def process_sequence(
    prof: FrameProfile,
    fasta: FASTAParser,
    record: RecordWriter,
    codon_writer: Union[FASTAWriter, None],
    amino_writer: Union[FASTAWriter, None],
    gcode,
):
    from nmm.sequence import Sequence
    from nmm.alphabet import DNAAlphabet

    for ti, tgt in enumerate(fasta):
        print()
        show_header2(f"Target {ti}")
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
            item_id = record.add_item(seqid, start, stop)
            codon_result = frag.decode()

            if codon_writer is not None:
                seq = codon_result.sequence.decode()
                codon_writer.write_item(item_id, seq)

            if amino_writer is not None:
                amino_result = codon_result.decode(gcode)
                seq = amino_result.sequence.decode()
                amino_writer.write_item(item_id, seq)


def finalize_stream(stream: LazyFile):
    if stream.name != "-":
        print(f"Writing to <{stream.name}> file.")

    stream.close_intelligently()


def write_target(file, defline: str, sequence: str):
    file.write(">" + defline + "\n")
    file.write(sequence + "\n")


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
