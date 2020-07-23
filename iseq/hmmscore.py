import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path

from iseq.file import fetch_file, make_executable
from iseq.tblout import TBLData

__all__ = ["HMMScore"]


_filemap = {
    "hmmsearch_macosx_10_9_x86_64": "6ac1afefae642933b19f4efd87bf869b58c5429e48645ebf0743ef2f03d2fd93",
    "hmmsearch_manylinux2010_x86_64": "74f19eadbb3b43747b569c5e708ed32b1191261c4e332c89c4fd9100fdfaf520",
    "hmmfetch_macosx_10_9_x86_64": "b3a5b96273114b39848778f89f629a5b27c58a10121eecd25afa8ec3792b5425",
    "hmmfetch_manylinux2010_x86_64": "05ef1eb01c01a0b61d8ba924ca46bb43051dca942d60aa5caab5f880a25728d3",
}


class HMMScore:
    def __init__(self, profile: Path):

        if sys.platform not in ["linux", "darwin"]:
            raise RuntimeError(f"Unsupported platform: {sys.platform}.")

        hmmsearch = "hmmsearch_manylinux2010_x86_64"
        if sys.platform == "darwin":
            hmmsearch = "hmmsearch_macosx_10_9_x86_64"

        hmmfetch = "hmmfetch_manylinux2010_x86_64"
        if sys.platform == "darwin":
            hmmfetch = "hmmfetch_macosx_10_9_x86_64"

        url_base = "https://iseq-py.s3.amazonaws.com/binhouse"
        hmmsearch_path = fetch_file(hmmsearch, "bin", url_base, _filemap)
        hmmfetch_path = fetch_file(hmmfetch, "bin", url_base, _filemap)

        make_executable(hmmsearch_path)
        make_executable(hmmfetch_path)
        self._hmmsearch_path = hmmsearch_path
        self._hmmfetch_path = hmmfetch_path

        self._profile = profile.absolute()
        if not self._is_indexed():
            self._index_it()

    def _index_it(self):
        cmd = [self._hmmfetch_path, "--index", self._profile]
        subprocess.check_call(cmd)

    def _is_indexed(self) -> bool:
        prof = self._profile
        if prof.with_suffix(prof.suffix + ".ssi").exists():
            return True
        return prof.with_suffix(prof.suffix + ".h3m.ssi").exists()

    def score(self, accession: str, seq: str, heuristics=True, cut_ga=False):

        with tempfile.TemporaryDirectory() as tmpdir:
            target = Path(tmpdir) / "target.fasta"
            with open(target, "w") as file:
                file.write(">Unknown\n")
                file.write(seq)

            cmd = f"{self._hmmfetch_path} {self._profile} {accession} | "
            cmd += f"{self._hmmsearch_path} -o /dev/null "
            max_flag = ""
            if not heuristics:
                max_flag = "--max"
            cut_ga_flag = ""
            if cut_ga:
                cut_ga_flag = "--cut_ga"
            cmd += f"--tblout /dev/stdout {max_flag} {cut_ga_flag} - {target}"

            output = subprocess.check_output(cmd, shell=True)
            tbldata = TBLData(StringIO(output.decode()))
            rows = list(tbldata.rows.values())
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
