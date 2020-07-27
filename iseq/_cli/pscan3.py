import re
from collections import OrderedDict
from io import StringIO
from pathlib import Path
from typing import Dict, List, Set, TextIO, Tuple

import click
from fasta_reader import FASTAItem, FASTAWriter, open_fasta
from hmmer import HMMER
from hmmer_reader import num_models, open_hmmer
from imm import Alphabet
from nmm import IUPACAminoAlphabet
from tqdm import tqdm

from iseq.alphabet import alphabet_name, infer_fasta_alphabet, infer_hmmer_alphabet
from iseq.codon_table import CodonTable
from iseq.domtblout import DomTBLData, DomTBLRow
from iseq.gff import read as read_gff
from iseq.hmmer_model import HMMERModel
from iseq.profile import ProfileID
from iseq.protein import ProteinProfile, create_profile2

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
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
@click.option(
    "--window",
    type=int,
    help="Window length. Defaults to zero, which means no window.",
    default=0,
)
@click.option(
    "--hit-prefix", help="Hit prefix. Defaults to `item`.", default="item", type=str,
)
def pscan3(
    profile: str,
    target: TextIO,
    epsilon: float,
    output: str,
    ocodon: str,
    oamino: str,
    quiet: bool,
    window: int,
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
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)

    target_abc = infer_target_alphabet(target)

    gcode = CodonTable(target_abc, IUPACAminoAlphabet())
    scan = PScan3(Path(profile), owriter, cwriter, awriter, gcode)
    scan.scan(target, window, epsilon, quiet)
    scan.close()


class PScan3:
    def __init__(
        self,
        profile: Path,
        output: OutputWriter,
        ocodon: FASTAWriter,
        oamino: FASTAWriter,
        codon_table: CodonTable,
    ):
        self._profile = profile
        self._output = output
        self._ocodon = ocodon
        self._oamino = oamino
        self._codon_table = codon_table

        hmmer = HMMER(profile)
        if not hmmer.is_indexed:
            hmmer.index()

        self._hmmer = hmmer

    @property
    def num_models(self) -> int:
        return num_models(self._profile)

    def scan(self, target: TextIO, window: int, epsilon: float, quiet: bool):
        with open_fasta(target) as fasta:
            targets = list(fasta)

        base_alphabet = self._codon_table.base_alphabet
        with open_hmmer(self._profile) as pfile:
            total = self.num_models
            profiles = tqdm(iter(pfile), desc="Models", total=total, disable=quiet)
            for plain_model in profiles:
                hmodel = HMMERModel(plain_model)
                prof = create_profile2(hmodel, base_alphabet, window, epsilon)
                self._scan_targets(prof, targets, epsilon, quiet)

    def _scan_targets(
        self,
        prof: ProteinProfile,
        targets: List[FASTAItem],
        epsilon: float,
        quiet: bool,
    ):
        wlen = prof.window_length
        for tgt in tqdm(targets, desc="Targets", leave=False, disable=quiet):
            seq = prof.create_sequence(tgt.sequence.encode())
            search_results = prof.search(seq)
            ifragments = search_results.ifragments()

            tgt_abc = seq.alphabet
            prof_abc = prof.alphabet
            for ifrag in ifragments:
                self._process_fragment(
                    ifrag, f"{tgt.id}", tgt_abc, prof.profid, prof_abc, wlen, epsilon
                )

    def _process_fragment(
        self,
        frag,
        target_id: str,
        target_abc: Alphabet,
        prof_id: ProfileID,
        prof_abc: Alphabet,
        window_length: int,
        epsilon: float,
    ):
        start = frag.interval.start
        stop = frag.interval.stop
        codon_frag = frag.fragment.decode()
        amino_frag = codon_frag.decode(self._codon_table)
        e_value, score, bias = target_score(
            self._hmmer, prof_id.acc, str(amino_frag.sequence)
        )
        if score == "NAN":
            return
        item_id = self._output.write_item(
            target_id,
            alphabet_name(target_abc),
            prof_id,
            alphabet_name(prof_abc),
            start,
            stop,
            window_length,
            {"Epsilon": epsilon, "E-value": e_value, "Score": score, "Bias": bias},
        )
        self._ocodon.write_item(item_id, str(codon_frag.sequence))
        self._oamino.write_item(item_id, str(amino_frag.sequence))

    def close(self):
        self._output.close()
        self._ocodon.close()
        self._oamino.close()


def infer_profile_alphabet(profile: TextIO):
    hmmer = open_hmmer(profile)
    hmmer_alphabet = infer_hmmer_alphabet(hmmer)
    profile.seek(0)
    if hmmer_alphabet is None:
        raise click.UsageError("Could not infer alphabet from PROFILE.")
    return hmmer_alphabet


def infer_target_alphabet(target: TextIO):
    fasta = open_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    target.seek(0)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    return target_alphabet


class ScoreTable:
    def __init__(self, domtbldata: DomTBLData):
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
    df = gff.to_dataframe()
    if len(df) == 0:
        return set()
    return set(df["att_ID"].tolist())


def update_fasta_file(filepath, target_set: Set[str]):
    with open_fasta(filepath) as fasta:
        targets = list(fasta)

    with FASTAWriter(filepath) as writer:
        for target in targets:
            tgt_id = target.id
            if tgt_id in target_set:
                writer.write_item(target.defline, target.sequence)


def translate_fasta_file(filepath, target_map: Dict[str, str]):
    with open_fasta(filepath) as fasta:
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


def target_score(hmmer: HMMER, acc: str, seq: str, heuristic=True, cut_ga=False):
    result = hmmer.search(
        StringIO(">Unknown\n" + seq),
        "/dev/null",
        tblout=True,
        heuristic=heuristic,
        cut_ga=cut_ga,
        hmmkey=acc,
    )

    rows = result.tbl
    if len(rows) == 0:
        e_value = "NAN"
        score = "NAN"
        bias = "NAN"
    else:
        assert len(rows) == 1
        e_value = rows[0].full_sequence.e_value
        score = rows[0].full_sequence.score
        bias = rows[0].full_sequence.bias

    return e_value, score, bias
