import sys
from io import StringIO
from pathlib import Path

from iseq.domtblout import DomTBLData
from iseq.file import fetch_file, make_executable, tmp_cwd

__all__ = ["HMMScan"]


_filemap = {
    "hmmscan_macosx_10_9_x86_64": "b83e438e307f9b0fbac0426dc0f7d6ca74d93969362764f90de1f6c9044025b3",
    "hmmscan_manylinux2010_x86_64": "c78a830d5019b27c178725900e11a79acd951ab088428ad102ebbf091628149a",
}


class HMMScan:
    def __init__(self):

        filename = "hmmscan"

        if sys.platform == "linux":
            filename = "hmmscan_manylinux2010_x86_64"
        elif sys.platform == "darwin":
            filename = "hmmscan_macosx_10_9_x86_64"
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}.")

        url_base = "https://iseq-py.s3.amazonaws.com/binhouse"
        prog_path = fetch_file(filename, "bin", url_base, _filemap)

        make_executable(prog_path)
        self._prog_path = prog_path

    def scan(self, profile: Path, target: Path) -> DomTBLData:
        import subprocess

        profile = profile.absolute()
        target = target.absolute()

        with tmp_cwd():
            cmd = [
                self._prog_path,
                "--cut_ga",
                "--domtblout",
                "domtblout",
                str(profile),
                str(target),
            ]
            subprocess.check_output(cmd)
            with open("domtblout", "r") as file:
                return DomTBLData(StringIO(file.read()))
