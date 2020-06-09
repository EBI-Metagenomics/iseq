import fileinput
import os
import re
from collections import OrderedDict
from typing import Dict, Mapping

import click

from iseq.file import tmp_cwd
from iseq.tblout import tblout_reader


@click.command()
@click.argument("profile", type=click.Path(exists=True))
@click.argument("target", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
def score(profile, target, gff):
    """
    Score sequence(s) against a protein profiles database.

    GFF is a file defining potential homologous sequences found in TARGET
    while matching it against PROFILE. The GFF file will be in-place edited
    to contain the additional field `E-value` estimated by HMMER3.
    """

    hmmsearch = HMMSearch()
    scores = hmmsearch.compute_scores(profile, target)
    update_gff_file(gff, scores)


class HMMSearch:
    def __init__(self):
        import distutils.spawn

        program = "hmmsearch"
        prog_path = distutils.spawn.find_executable(program)
        if prog_path is None:
            raise click.ClickException(f"Could not find the `{program}` program.")

        self._prog_path = prog_path

    def compute_scores(self, profile, target):
        import subprocess

        profile = os.path.abspath(profile)
        target = os.path.abspath(target)

        with tmp_cwd():
            cmd = [self._prog_path, "--tblout", "tblout", profile, target]
            subprocess.check_output(cmd)
            scores: Dict[str, float] = {}
            with open("tblout", "r") as file:
                for row in tblout_reader(file):
                    scores[row.target_name] = row.full_sequence.e_value

        return scores


def update_gff_file(filepath, scores: Mapping[str, float]):

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

        fields_list = []
        for v in right.split(";"):
            name, value = v.split("=", 1)
            fields_list.append((name, value))

        attr = OrderedDict(fields_list)
        if "ID" not in attr:
            print(row)
            continue

        if attr["ID"] not in scores:
            if "E-value" in attr:
                del attr["E-value"]
        else:
            attr["E-value"] = str(scores[attr["ID"]])

        print(left + ";".join(k + "=" + v for k, v in attr.items()))
