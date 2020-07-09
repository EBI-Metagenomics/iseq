import sys
from io import StringIO
from pathlib import Path

from iseq.file import fetch_file, make_executable, tmp_cwd

__all__ = ["HMMPress"]

_filemap = {
    "hmmpress_macosx_10_9_x86_64": "0e1683cd4556382a02b0c49b3d60126d2c1da13fae858e60d83faac4906637ec",
    "hmmpress_manylinux2010_x86_64": "e23676c6425b0368d1ffd8de6f7897903ac32c3ad61be40ef60fa55712b3b222",
}


class HMMPress:
    def __init__(self):

        filename = "hmmpress"

        if sys.platform == "linux":
            filename = "hmmpress_manylinux2010_x86_64"
        elif sys.platform == "darwin":
            filename = "hmmpress_macosx_10_9_x86_64"
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}.")

        url_base = "https://iseq-py.s3.amazonaws.com/binhouse"
        prog_path = fetch_file(filename, "bin", url_base, _filemap)

        make_executable(prog_path)
        self._prog_path = prog_path

    def press(self, profile: Path):
        import subprocess

        profile = profile.absolute()

        with tmp_cwd():
            cmd = [
                self._prog_path,
                str(profile),
            ]
            subprocess.check_output(cmd)

    @staticmethod
    def pressed(profile: Path) -> bool:
        profile = profile.absolute()
        exts = [".h3f", ".h3i", ".h3m", ".h3p"]
        return all([Path(str(profile) + ext).exists() for ext in exts])
