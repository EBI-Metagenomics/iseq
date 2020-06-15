from pathlib import Path

from iseq.file import tmp_cwd
from iseq.tblout import TBLData

__all__ = ["HMMSearch"]


class HMMSearch:
    def __init__(self):
        import distutils.spawn

        program = "hmmsearch"
        prog_path = distutils.spawn.find_executable(program)
        if prog_path is None:
            raise RuntimeError(f"Could not find the `{program}` program.")

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
