from typing import List

from fasta_reader import FASTAItem
from hmmer_reader import HMMERParser

from iseq.hmmdata import HMMData
from iseq.hmmer3 import create_profile
from iseq.model import EntryDistr

from .scanner import IntFrag, OutputWriter, Scanner


class HMMER3Scanner(Scanner):
    def __init__(
        self,
        output_writer: OutputWriter,
        window_length: int,
        stdout,
        hmmer3_compat: bool,
        entry_distr: EntryDistr,
    ):
        self._hmmer3_compat = hmmer3_compat
        self._entry_distr = entry_distr
        super().__init__(output_writer, window_length, stdout)

    def process_profile(self, profile_parser: HMMERParser, targets: List[FASTAItem]):

        mt = dict(profile_parser.metadata)
        self._output_writer.profile = mt.get("ACC", mt.get("NAME", "UNKNOWN"))
        # breakpoint()
        hmmdata = HMMData(profile_parser)
        prof = create_profile(hmmdata, self._hmmer3_compat, self._entry_distr)
        self._scan_targets(prof, targets)

    def _write_fragments(self, seqid: str, ifragments: List[IntFrag]):
        for ifrag in ifragments:
            start = ifrag.interval.start
            stop = ifrag.interval.stop
            self._output_writer.write_item(seqid, start, stop)
