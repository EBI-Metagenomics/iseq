from pathlib import Path
from typing import Dict, NamedTuple

from iseq.file import tmp_cwd
from iseq.tblout import tblout_reader, TBLData

__all__ = ["HMMSearch"]


# NamedTuple("Ro",
# target name        accession  query name           accession

# target_name=row[0],
# target_accession=row[1],
# query_name=row[2],
# query_accession=row[3],


class HMMSearch:
    def __init__(self):
        import distutils.spawn

        program = "hmmsearch"
        prog_path = distutils.spawn.find_executable(program)
        if prog_path is None:
            raise RuntimeError(f"Could not find the `{program}` program.")

        self._prog_path = prog_path

    def compute_scores(self, profile: Path, target: Path):
        import subprocess

        profile = profile.absolute()
        target = target.absolute()

        with tmp_cwd():
            cmd = [self._prog_path, "--tblout", "tblout", str(profile), str(target)]
            subprocess.check_output(cmd)
            scores: Dict[str, str] = {}
            with open("tblout", "r") as file:
                for row in tblout_reader(file):
                    breakpoint()
                    scores[row.target_name] = row.full_sequence.e_value

        return scores

    def search(self, profile: Path, target: Path) -> TBLData:
        import subprocess

        profile = profile.absolute()
        target = target.absolute()

        with tmp_cwd():
            cmd = [self._prog_path, "--tblout", "tblout", str(profile), str(target)]
            subprocess.check_output(cmd)
            with open("tblout", "r") as file:
                return TBLData(file)
