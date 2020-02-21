from contextlib import contextmanager

import click


@click.command()
@click.argument("profile", type=click.Path(exists=True))
@click.argument("target", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
def score(profile, target, gff):
    """
    Score sequences against a database of HMMER3 profiles.
    """

    hmmsearch = HMMSearch()
    scores = hmmsearch.compute_scores(profile, target)
    hmmsearch.update_gff_file(gff, scores)


class HMMSearch:
    def __init__(self):
        import distutils.spawn

        program = "hmmsearch"
        prog_path = distutils.spawn.find_executable(program)
        if prog_path is None:
            raise click.ClickException(f"Could not find the `{program}` program.")

        self._prog_path = prog_path

    def compute_scores(self, profile, target):
        import os
        import subprocess
        from ..io import tblout_reader

        profile = os.path.abspath(profile)
        target = os.path.abspath(target)

        with tmp_cwd():
            cmd = [self._prog_path, "--tblout", "tblout", profile, target]
            subprocess.check_output(cmd)
            scores = {}
            with open("tblout", "r") as file:
                for row in tblout_reader(file):
                    scores[row.target_name] = row.full_sequence.e_value

        return scores

    def update_gff_file(self, filepath, scores):
        import fileinput
        import re
        from collections import OrderedDict

        for row in fileinput.input(filepath, inplace=True, backup=".bak"):
            row = row.rstrip()
            if row.startswith("#"):
                print(row)
                continue

            match = re.match(r"^(.+\t)([^\t]+)$", row)
            if match is None:
                print(row)
                continue

            left = match.group(1)
            right = match.group(2)

            if right == ".":
                print(row)
                continue

            attr = OrderedDict(v.split("=") for v in right.split(";"))
            if "ID" not in attr:
                print(row)
                continue

            if attr["ID"] not in scores:
                if "E-value" in attr:
                    del attr["E-value"]
            else:
                attr["E-value"] = scores[attr["ID"]]

            print(left + ";".join(k + "=" + v for k, v in attr.items()))


@contextmanager
def tmp_cwd():
    import os
    from tempfile import TemporaryDirectory

    oldpwd = os.getcwd()
    with TemporaryDirectory() as tmpdir:

        os.chdir(tmpdir)
        try:
            yield
        finally:
            os.chdir(oldpwd)
