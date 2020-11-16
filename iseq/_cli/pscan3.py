from io import StringIO
from pathlib import Path
from typing import Iterable, NamedTuple, TextIO

import click
from fasta_reader import FASTAItem, FASTAWriter, read_fasta
from hmmer import HMMER
from hmmer_reader import num_models, open_hmmer
from imm import Alphabet
from nmm import IUPACAminoAlphabet
from tqdm import tqdm

from iseq.alphabet import alphabet_name, infer_fasta_alphabet
from iseq.codon_table import CodonTable
from iseq.hmmer_model import HMMERModel
from iseq.profile import ProfileID
from iseq.protein import ProteinProfile, create_profile2

from .output_writer import OutputWriter

HMMEROptions = NamedTuple("HMMEROptions", [("heuristic", bool), ("cut_ga", bool)])


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
    "--hit-prefix",
    help="Hit prefix. Defaults to `item`.",
    default="item",
    type=str,
)
@click.option(
    "--heuristic/--no-heuristic",
    help="Enable HMMER heuristics. Defaults to True.",
    default=True,
)
@click.option(
    "--cut-ga/--no-cut-ga",
    help="Enable use of profile's GA gathering cutoffs to set all thresholding. Defaults to False.",
    default=False,
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
    heuristic: bool,
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

    target_abc = infer_target_alphabet(target)

    gcode = CodonTable(target_abc, IUPACAminoAlphabet())
    opts = HMMEROptions(heuristic, cut_ga)
    scan = PScan3(Path(profile), owriter, cwriter, awriter, gcode, opts)
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
        hmmer_options: HMMEROptions,
    ):
        self._profile = profile
        self._output = output
        self._ocodon = ocodon
        self._oamino = oamino
        self._codon_table = codon_table

        hmmer = HMMER(profile)
        hmmer.timeout = 60
        if not hmmer.is_indexed:
            hmmer.index()

        self._hmmer = hmmer
        self._hmmer_options = hmmer_options

    @property
    def num_models(self) -> int:
        return num_models(self._profile)

    def scan(self, target: TextIO, window: int, epsilon: float, quiet: bool):
        with read_fasta(target) as fasta:
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
        targets: Iterable[FASTAItem],
        epsilon: float,
        quiet: bool,
    ):
        wlen = prof.window_length
        prof_abc = prof.alphabet
        frags = {}
        targets_abc = {}
        for tgt in tqdm(targets, desc="Targets", leave=False, disable=quiet):
            seq = prof.create_sequence(tgt.sequence.encode())
            search_results = prof.search(seq)
            ifragments = search_results.ifragments()

            targets_abc[tgt.id] = seq.alphabet
            for i, ifrag in enumerate(ifragments):
                frags[(i, f"{tgt.id}")] = ifrag

        scores = self._score_fragments(prof, frags)
        for (i, tgt_id), frag in frags.items():
            score = scores.get((i, tgt_id), None)
            if score is None:
                continue
            self._process_fragment(
                frag,
                tgt_id,
                targets_abc[tgt_id],
                prof.profid,
                prof_abc,
                wlen,
                epsilon,
                score[0],
                score[1],
                score[2],
            )

    def _score_fragments(self, prof: ProteinProfile, frags):
        targets = []
        for (i, tgt_id), frag in frags.items():
            codon_frag = frag.fragment.decode()
            amino_frag = codon_frag.decode(self._codon_table)
            seqid = "_".join([str(i), tgt_id])
            targets.append(f">{seqid}\n" + str(amino_frag.sequence))

        result = self._hmmer.search(
            StringIO("\n".join(targets)),
            "/dev/null",
            tblout=True,
            heuristic=self._hmmer_options.heuristic,
            cut_ga=self._hmmer_options.cut_ga,
            hmmkey=prof.profid.acc,
            Z=1,
        )

        scores = {}
        for row in result.tbl:
            e_value = row.full_sequence.e_value
            score = row.full_sequence.score
            if score.lower() == "nan":
                continue
            bias = row.full_sequence.bias
            istr, tgt_id = row.target.name.partition("_")[::2]
            scores[(int(istr), tgt_id)] = (e_value, score, bias)

        return scores

    def _process_fragment(
        self,
        frag,
        target_id: str,
        target_abc: Alphabet,
        prof_id: ProfileID,
        prof_abc: Alphabet,
        window_length: int,
        epsilon: float,
        e_value: str,
        score: str,
        bias: str,
    ):
        start = frag.interval.start
        stop = frag.interval.stop
        codon_frag = frag.fragment.decode()
        amino_frag = codon_frag.decode(self._codon_table)
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


def infer_target_alphabet(target: TextIO):
    fasta = read_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    target.seek(0)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    return target_alphabet
