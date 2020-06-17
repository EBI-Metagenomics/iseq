import os
import re
from collections import OrderedDict
from pathlib import Path
from typing import IO

import click
from fasta_reader import FASTAWriter, open_fasta
from hmmer_reader import open_hmmer
from nmm import AminoAlphabet, BaseAlphabet, CanonicalAminoAlphabet, GeneticCode

from iseq.alphabet import infer_fasta_alphabet, infer_hmmer_alphabet
from iseq.hmmer_model import HMMERModel
from iseq.hmmsearch import HMMSearch
from iseq.protein import create_profile
from iseq.tblout import TBLData

from .debug_writer import DebugWriter
from .misc import consolidate
from .output_writer import OutputWriter


@click.command()
@click.argument(
    "profile",
    type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True),
)
@click.argument("target", type=click.File("r"))
@click.option(
    "--epsilon", type=float, default=1e-2, help="Indel probability. Defaults to 1e-2."
)
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
def pscan(
    profile,
    target,
    epsilon: float,
    output,
    ocodon,
    oamino,
    quiet,
    window: int,
    odebug,
    e_value: bool,
):
    """
    Search nucleotide sequence(s) against a protein profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """

    owriter = OutputWriter(output)
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)
    dwriter = DebugWriter(odebug)

    if quiet:
        click.open_file(os.devnull, "a")
    else:
        click.get_text_stream("stdout")

    with open(profile, "r") as file:
        profile_abc = _infer_profile_alphabet(file)
    target_abc = _infer_target_alphabet(target)

    assert isinstance(target_abc, BaseAlphabet) and isinstance(
        profile_abc, AminoAlphabet
    )

    gcode = GeneticCode(target_abc, CanonicalAminoAlphabet())

    with open_fasta(target) as fasta:
        targets = list(fasta)

    for plain_model in open_hmmer(profile):
        model = HMMERModel(plain_model)
        prof = create_profile(model, gcode.base_alphabet, window, epsilon)
        for tgt in targets:
            seq = prof.create_sequence(tgt.sequence.encode())
            search_results = prof.search(seq)
            intfrags, debug_list = consolidate(search_results)
            seqid = f"{tgt.defline.split()[0]}"
            for intfrag in intfrags:
                start = intfrag.interval.start
                stop = intfrag.interval.stop
                item_id = owriter.write_item(
                    seqid,
                    model.model_id,
                    start,
                    stop,
                    prof.window_length,
                    {"Epsilon": epsilon},
                )
                codon_result = intfrag.fragment.decode()
                cwriter.write_item(item_id, str(codon_result.sequence))
                amino_result = codon_result.decode(gcode)
                awriter.write_item(item_id, str(amino_result.sequence))
            for debug_item in debug_list:
                dwriter.write_row(seqid, *debug_item)

    owriter.close()
    cwriter.close()
    awriter.close()
    odebug.close_intelligently()

    if e_value:
        hmmsearch = HMMSearch()
        tbldata = hmmsearch.search(Path(profile), Path(oamino))
        update_gff_file(output, tbldata)


def _infer_profile_alphabet(profile: IO[str]):
    hmmer = open_hmmer(profile)
    hmmer_alphabet = infer_hmmer_alphabet(hmmer)
    profile.seek(0)
    if hmmer_alphabet is None:
        raise click.UsageError("Could not infer alphabet from PROFILE.")
    return hmmer_alphabet


def _infer_target_alphabet(target: IO[str]):
    fasta = open_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    target.seek(0)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    return target_alphabet


def update_gff_file(filepath, tbldata: TBLData):
    import in_place

    with in_place.InPlace(filepath) as file:
        for row in file:
            row = row.rstrip()
            if row.startswith("#"):
                file.write(row)
                file.write("\n")
                continue

            match = re.match(r"^(.+\t)([^\t]+)$", row)
            if match is None:
                file.write(row)
                file.write("\n")
                continue

            left = match.group(1)
            right = match.group(2)

            if right == ".":
                file.write(row)
                file.write("\n")
                continue

            fields_list = []
            for v in right.split(";"):
                name, value = v.split("=", 1)
                fields_list.append((name, value))

            attr = OrderedDict(fields_list)

            data = tbldata.find(
                target_name=attr["ID"],
                query_name=attr["Profile_name"],
                query_acc=attr["Profile_acc"],
            )
            if data is None:
                file.write(row)
                file.write("\n")
                continue

            attr["E-value"] = data.full_sequence.e_value
            file.write(left + ";".join(k + "=" + v for k, v in attr.items()))
            file.write("\n")
