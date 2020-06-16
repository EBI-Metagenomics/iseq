from typing import List

from fasta_reader import FASTAItem, FASTAWriter
from hmmer_reader import HMMERParser
from nmm import GeneticCode

from iseq.protein import create_profile, ProteinProfile

from .debug_writer import DebugWriter
from .scanner import IntFrag, OutputWriter, Scanner


class ProteinScanner(Scanner):
    def __init__(
        self,
        output_writer: OutputWriter,
        debug_writer: DebugWriter,
        codon_writer: FASTAWriter,
        amino_writer: FASTAWriter,
        genetic_code: GeneticCode,
        epsilon: float,
        stdout,
    ):
        super().__init__(output_writer, debug_writer, stdout)
        self._codon_writer = codon_writer
        self._amino_writer = amino_writer
        self._genetic_code = genetic_code
        self._epsilon = epsilon

    def process_profile(
        self, profile_parser: HMMERParser, targets: List[FASTAItem], window_length: int
    ):

        mt = dict(profile_parser.metadata)
        self._output_writer.set_profile(
            name=mt.get("NAME", "-"), acc=mt.get("ACC", "-")
        )
        base_alphabet = self._genetic_code.base_alphabet
        # breakpoint()
        prof = create_profile(
            profile_parser, base_alphabet, window_length, self._epsilon
        )
        # print(prof.alt_model)
        # TODO: remove it
        # import sys

        # sys.exit(1)
        # breakpoint()
        self._scan_targets(prof, targets)

    def _write_fragments(
        self, profile: ProteinProfile, target: FASTAItem, ifragments: List[IntFrag]
    ):
        seqid = f"{target.defline.split()[0]}"
        for ifrag in ifragments:
            start = ifrag.interval.start
            stop = ifrag.interval.stop
            item_id = self._output_writer.write_item(
                seqid, start, stop, profile.window_length, {"Epsilon": self._epsilon}
            )

            codon_result = ifrag.fragment.decode()
            self._codon_writer.write_item(item_id, str(codon_result.sequence))

            amino_result = codon_result.decode(self._genetic_code)
            self._amino_writer.write_item(item_id, str(amino_result.sequence))
