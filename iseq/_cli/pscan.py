import os
import re
import sys
from collections import OrderedDict
from typing import IO, List

import click
from fasta_reader import FASTAWriter, read_fasta
from hmmer import HMMER
from hmmer.typing import TBLRow
from hmmer_reader import num_models, open_hmmer
from nmm import AminoAlphabet, BaseAlphabet, IUPACAminoAlphabet
from tqdm import tqdm

from iseq.alphabet import alphabet_name, infer_fasta_alphabet, infer_hmmer_alphabet
from iseq.codon_table import CodonTable
from iseq.hmmer_model import HMMERModel
from iseq.protein import create_profile, create_profile2

from .debug_writer import DebugWriter
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
    "--quiet/--no-quiet",
    "-q/-nq",
    help="Disable standard output.",
    default=False,
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
@click.option(
    "--model",
    type=click.Choice(["1", "2"]),
    help="Model 1 or 2. Defaults 1 for now",
    default="1",
)
@click.option(
    "--hit-prefix",
    help="Hit prefix. Defaults to `item`.",
    default="item",
    type=str,
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
    model: str,
    hit_prefix: str,
):
    """
    Search nucleotide sequence(s) against a protein profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """

    owriter = OutputWriter(output, item_prefix=hit_prefix)
    cwriter = FASTAWriter(ocodon, sys.maxsize)
    awriter = FASTAWriter(oamino, sys.maxsize)
    dwriter = DebugWriter(odebug)

    with open(profile, "r") as file:
        profile_abc = _infer_profile_alphabet(file)
    target_abc = infer_target_alphabet(target)

    assert isinstance(target_abc, BaseAlphabet) and isinstance(
        profile_abc, AminoAlphabet
    )

    gcode = CodonTable(target_abc, IUPACAminoAlphabet())

    with read_fasta(target) as fasta:
        targets = list(fasta)

    total = num_models(profile)
    for plain_model in tqdm(
        open_hmmer(profile), desc="Models", total=total, disable=quiet
    ):
        hmodel = HMMERModel(plain_model)
        if model == "1":
            prof = create_profile(hmodel, gcode.base_alphabet, window, epsilon)
        else:
            prof = create_profile2(hmodel, gcode.base_alphabet, window, epsilon)
            assert model == "2"

        for tgt in tqdm(targets, desc="Targets", leave=False, disable=quiet):
            seq = prof.create_sequence(tgt.sequence.encode())
            search_results = prof.search(seq)
            ifragments = search_results.ifragments()
            seqid = f"{tgt.id}"

            for ifrag in ifragments:
                start = ifrag.interval.start
                stop = ifrag.interval.stop
                item_id = owriter.write_item(
                    seqid,
                    alphabet_name(seq.alphabet),
                    prof.profid,
                    alphabet_name(prof.alphabet),
                    start,
                    stop,
                    prof.window_length,
                    {"Epsilon": epsilon},
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
        hmmer = HMMER(profile)
        result = hmmer.search(oamino, "/dev/null", tblout=True)
        update_gff_file(output, result.tbl)


def _infer_profile_alphabet(profile: IO[str]):
    hmmer = open_hmmer(profile)
    hmmer_alphabet = infer_hmmer_alphabet(hmmer)
    profile.seek(0)
    if hmmer_alphabet is None:
        raise click.UsageError("Could not infer alphabet from PROFILE.")
    return hmmer_alphabet


def infer_target_alphabet(target: IO[str]):
    fasta = read_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    target.seek(0)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    return target_alphabet


def update_gff_file(filepath, tbldata: List[TBLRow]):
    import in_place

    tbl = {}
    for row in iter(tbldata):
        key = (row.target.name, row.query.name, row.query.accession)
        tbl[key] = row

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

            key = (attr["ID"], attr["Profile_name"], attr["Profile_acc"])
            if key not in tbl:
                file.write(row)
                file.write("\n")
                continue

            attr["E-value"] = tbl[key].full_sequence.e_value
            file.write(left + ";".join(k + "=" + v for k, v in attr.items()))
            file.write("\n")
