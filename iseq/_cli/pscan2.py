import os
import re
from collections import OrderedDict
from typing import IO, Dict, List, Set, Tuple

import click
from fasta_reader import FASTAWriter, read_fasta
from hmmer import HMMER
from hmmer.typing import DomTBLRow
from hmmer_reader import num_models, open_hmmer
from nmm import AminoAlphabet, BaseAlphabet, IUPACAminoAlphabet
from tqdm import tqdm

from iseq.alphabet import alphabet_name, infer_fasta_alphabet, infer_hmmer_alphabet
from iseq.codon_table import CodonTable
from iseq.gff import read as read_gff
from iseq.hmmer_model import HMMERModel
from iseq.protein import create_profile2

from .debug_writer import DebugWriter
from .output_writer import OutputWriter


@click.command()
@click.argument(
    "profile",
    type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=False),
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
    "--max-e-value",
    type=float,
    default=10.0,
    help="Filter out items for which E-value > --max-e-evalue.",
)
@click.option(
    "--hit-prefix",
    help="Hit prefix. Defaults to `item`.",
    default="item",
    type=str,
)
@click.option(
    "--cut-ga/--no-cut-ga",
    help="Enable use of profile's GA gathering cutoffs to set all thresholding. Defaults to True.",
    default=True,
)
def pscan2(
    profile,
    target,
    epsilon: float,
    output,
    ocodon,
    oamino,
    quiet,
    window: int,
    odebug,
    max_e_value: float,
    hit_prefix: str,
    cut_ga: bool,
):
    """
    Search nucleotide sequence(s) against a protein profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """

    owriter = OutputWriter(output, item_prefix=hit_prefix)
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)
    dwriter = DebugWriter(odebug)

    with open(profile, "r") as file:
        profile_abc = infer_profile_alphabet(file)
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
        prof = create_profile2(hmodel, gcode.base_alphabet, window, epsilon)

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

    if not quiet:
        click.echo("Computing e-values... ", nl=False)
    hmmer = HMMER(profile)
    if not hmmer.is_pressed:
        hmmer.press()
    result = hmmer.scan(oamino, "/dev/null", domtblout=True, cut_ga=cut_ga)
    score_table = ScoreTable(result.domtbl)

    update_gff_file(output, score_table, max_e_value)
    target_set = target_set_from_gff_file(output)
    update_fasta_file(ocodon, target_set)
    update_fasta_file(oamino, target_set)
    if not quiet:
        click.echo("done.")

    def rep(k: str):
        return int(k.replace(hit_prefix, ""))

    if not quiet:
        click.echo("Translating files... ", nl=False)
    items = list(sorted(list(target_set), key=rep))
    target_map = {item: f"{hit_prefix}{i+1}" for i, item in enumerate(items)}
    translate_gff_file(output, target_map)
    translate_fasta_file(ocodon, target_map)
    translate_fasta_file(oamino, target_map)
    if not quiet:
        click.echo("done.")


def infer_profile_alphabet(profile: IO[str]):
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


class ScoreTable:
    def __init__(self, domtbldata: List[DomTBLRow]):
        self._tbldata: Dict[Tuple[str, str, str], DomTBLRow] = {}
        for line in domtbldata:
            key = (line.query.name, line.target.name, line.target.accession)
            if key in self._tbldata:
                prev_score = float(self._tbldata[key].full_sequence.score)
                if float(line.full_sequence.score) > prev_score:
                    self._tbldata[key] = line
            else:
                self._tbldata[key] = line

    def e_value(self, target_id: str, profile_name: str, profile_acc: str) -> str:
        key = (target_id, profile_name, profile_acc)
        return self._tbldata[key].full_sequence.e_value

    def has(self, target_id: str, profile_name: str, profile_acc: str) -> bool:
        key = (target_id, profile_name, profile_acc)
        return key in self._tbldata


def update_gff_file(filepath, score_table: ScoreTable, max_e_value: float):
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

            key = (attr["ID"], attr["Profile_name"], attr["Profile_acc"])
            if not score_table.has(*key):
                continue

            attr["E-value"] = score_table.e_value(*key)
            if float(attr["E-value"]) > max_e_value:
                continue

            file.write(left + ";".join(k + "=" + v for k, v in attr.items()))
            file.write("\n")


def translate_gff_file(filepath, target_map: Dict[str, str]):
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

            attr["ID"] = target_map[attr["ID"]]
            file.write(left + ";".join(k + "=" + v for k, v in attr.items()))
            file.write("\n")


def target_set_from_gff_file(filepath) -> Set[str]:
    gff = read_gff(filepath)
    gff.ravel()
    df = gff.dataframe
    if len(df) == 0:
        return set()
    return set(df["att_ID"].tolist())


def update_fasta_file(filepath, target_set: Set[str]):
    with read_fasta(filepath) as fasta:
        targets = list(fasta)

    with FASTAWriter(filepath) as writer:
        for target in targets:
            tgt_id = target.id
            if tgt_id in target_set:
                writer.write_item(target.defline, target.sequence)


def translate_fasta_file(filepath, target_map: Dict[str, str]):
    with read_fasta(filepath) as fasta:
        targets = list(fasta)

    with FASTAWriter(filepath) as writer:
        for target in targets:
            old_tgt_id = target.id
            tgt_id = target_map[old_tgt_id]
            if target.has_desc:
                defline = tgt_id + " " + target.desc
            else:
                defline = tgt_id
            writer.write_item(defline, target.sequence)
