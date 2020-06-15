import os
import re
from collections import OrderedDict
from pathlib import Path
from typing import Mapping, Optional

import click
from fasta_reader import open_fasta
from hmmer_reader import open_hmmer

from iseq.hmmsearch import HMMSearch
from iseq.model import EntryDistr
from iseq.tblout import TBLData

from .debug_writer import DebugWriter


@click.command()
@click.argument("profile", type=click.Path(exists=True, dir_okay=False))
@click.argument("target", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--output",
    type=click.Path(exists=False, writable=True, dir_okay=False),
    help="Save results to OUTPUT (GFF format).",
    default=os.devnull,
)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
@click.option(
    "--window",
    type=int,
    help="Window length. Defaults to zero, which means no window.",
    default=0,
)
@click.option(
    "--hmmer3-compat/--no-hmmer3-compat",
    help="Enable full HMMER3 compatibility. Defaults to False.",
    default=False,
)
@click.option(
    "--entry-distr",
    type=click.Choice(["uniform", "occupancy"], case_sensitive=False),
    help="Set the entry distribution. Defaults to occupancy.",
    default="occupancy",
)
@click.option(
    "--odebug",
    type=click.File("w"),
    help="Save debug info into a tab-separated values file.",
    default=os.devnull,
)
@click.option(
    "--e-value/--no-e-value",
    help="Enable E-value computation. Defaults to False.",
    default=False,
)
def hscan(
    profile: str,
    target: str,
    output,
    quiet,
    window: int,
    hmmer3_compat: bool,
    entry_distr: str,
    odebug,
    e_value: bool,
):
    """
    Search nucleotide sequence(s) against a profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """
    from .scanner import OutputWriter
    from .hmmer3_scanner import HMMER3Scanner

    owriter = OutputWriter(output, window)
    dwriter = DebugWriter(odebug)

    if entry_distr == "occupancy":
        edistr = EntryDistr.OCCUPANCY
    else:
        edistr = EntryDistr.UNIFORM

    if quiet:
        stdout = click.open_file(os.devnull, "a")
    else:
        stdout = click.get_text_stream("stdout")

    scanner: Optional[HMMER3Scanner] = None

    # if profile_abc.symbols != target_abc.symbols:
    #     raise click.UsageError("Alphabets mismatch.")

    scanner = HMMER3Scanner(owriter, dwriter, window, stdout, hmmer3_compat, edistr)

    with open_fasta(target) as fasta:
        targets = list(fasta)

    for hmmprof in open_hmmer(profile):
        scanner.show_profile_parser(hmmprof)
        scanner.process_profile(hmmprof, targets)

    # breakpoint()

    # scanner.finalize_stream("output", output)
    owriter.close()

    # breakpoint()
    if e_value:
        hmmsearch = HMMSearch()
        tbldata = hmmsearch.search(Path(profile), Path(target))
        update_gff_file(output, tbldata)

    scanner.finalize_stream("odebug", odebug)


def update_gff_file(filepath, tbldata: TBLData):
    import in_place

    with in_place.InPlace(filepath, backup_ext=".bak") as file:
        # for row in fileinput.input(filepath, inplace=True, backup=".bak"):
        for row in file:
            row = row.rstrip()
            if row.startswith("#"):
                file.write(row)
                file.write("\n")
                continue

            match = re.match(r"^(.+\t)([^\t]+)$", row)
            if match is None:
                file.write(row)
                file.write("\n")
                continue

            left = match.group(1)
            right = match.group(2)

            if right == ".":
                file.write(row)
                file.write("\n")
                continue

            fields_list = []
            for v in right.split(";"):
                name, value = v.split("=", 1)
                fields_list.append((name, value))

            attr = OrderedDict(fields_list)
            if "ID" not in attr:
                file.write(row)
                file.write("\n")
                continue

            if attr["ID"] not in scores:
                if "E-value" in attr:
                    del attr["E-value"]
            else:
                attr["E-value"] = str(scores[attr["ID"]])

            file.write(left + ";".join(k + "=" + v for k, v in attr.items()))
            file.write("\n")
