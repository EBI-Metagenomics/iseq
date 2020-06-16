from typing import List

from fasta_reader import FASTAItem
from hmmer_reader import HMMERParser

from iseq.hmmdata import HMMData
from iseq.hmmer3 import create_profile, HMMER3Profile
from iseq.model import EntryDistr

from .debug_writer import DebugWriter
from .scanner import IntFrag, OutputWriter, Scanner


class HMMER3Scanner(Scanner):
    def __init__(
        self,
        output_writer: OutputWriter,
        debug_writer: DebugWriter,
        stdout,
        hmmer3_compat: bool,
        entry_distr: EntryDistr,
    ):
        self._hmmer3_compat = hmmer3_compat
        self._entry_distr = entry_distr
        super().__init__(output_writer, debug_writer, stdout)

    def process_profile(
        self, profile_parser: HMMERParser, targets: List[FASTAItem], window_length: int
    ):

        mt = dict(profile_parser.metadata)
        self._output_writer.set_profile(
            name=mt.get("NAME", "-"), acc=mt.get("ACC", "-")
        )
        hmmdata = HMMData(profile_parser)
        prof = create_profile(
            hmmdata, self._hmmer3_compat, self._entry_distr, window_length
        )
        self._scan_targets(prof, targets)

    def _write_fragments(
        self, profile: HMMER3Profile, target: FASTAItem, ifragments: List[IntFrag]
    ):
        seqid = f"{target.defline.split()[0]}"
        for ifrag in ifragments:
            start = ifrag.interval.start
            stop = ifrag.interval.stop
            self._output_writer.write_item(seqid, start, stop, profile.window_length)
