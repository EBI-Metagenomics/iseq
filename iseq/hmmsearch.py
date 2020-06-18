import sys
from pathlib import Path

from iseq.file import fetch_file, tmp_cwd, make_executable
from iseq.tblout import TBLData

__all__ = ["HMMSearch"]


_filemap = {
    "hmmsearch_macosx_10_9_x86_64": "6ac1afefae642933b19f4efd87bf869b58c5429e48645ebf0743ef2f03d2fd93",
    "hmmsearch_manylinux2010_x86_64": "74f19eadbb3b43747b569c5e708ed32b1191261c4e332c89c4fd9100fdfaf520",
}


class HMMSearch:
    def __init__(self):

        filename = "hmmsearch"

        if sys.platform == "linux":
            filename = "hmmsearch_manylinux2010_x86_64"
        elif sys.platform == "darwin":
            filename = "hmmsearch_macosx_10_9_x86_64"
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}.")

        url_base = "https://iseq-py.s3.amazonaws.com/binhouse"
        prog_path = fetch_file(filename, "bin", url_base, _filemap)

        make_executable(prog_path)
        self._prog_path = prog_path

    def search(self, profile: Path, target: Path) -> TBLData:
        import subprocess

        profile = profile.absolute()
        target = target.absolute()

        with tmp_cwd():
            cmd = [self._prog_path, "--tblout", "tblout", str(profile), str(target)]
            subprocess.check_output(cmd)
            with open("tblout", "r") as file:
                return TBLData(file)
