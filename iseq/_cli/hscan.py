import os

import click
from fasta_reader import read_fasta
from hmmer_reader import open_hmmer

from iseq.alphabet import alphabet_name
from iseq.hmmer3 import create_profile
from iseq.hmmer_model import HMMERModel
from iseq.model import EntryDistr

from .debug_writer import DebugWriter
from .output_writer import OutputWriter


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option(
    "--output",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save results to OUTPUT (GFF format).",
    default="output.gff",
)
@click.option(
    "--quiet/--no-quiet",
    "-q/-nq",
    help="Disable standard output.",
    default=False,
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
def hscan(
    profile,
    target,
    output,
    quiet,
    window: int,
    hmmer3_compat: bool,
    entry_distr: str,
    odebug,
):
    """
    Search nucleotide sequence(s) against a profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """
    owriter = OutputWriter(output)
    dwriter = DebugWriter(odebug)

    if entry_distr == "occupancy":
        edistr = EntryDistr.OCCUPANCY
    else:
        edistr = EntryDistr.UNIFORM

    if quiet:
        click.open_file(os.devnull, "a")
    else:
        click.get_text_stream("stdout")

    with read_fasta(target) as fasta:
        targets = list(fasta)

    for plain_model in open_hmmer(profile):
        model = HMMERModel(plain_model)
        prof = create_profile(model, hmmer3_compat, edistr, window)
        for tgt in targets:
            seq = prof.create_sequence(tgt.sequence.encode())
            search_results = prof.search(seq)
            ifragments = search_results.ifragments()
            seqid = f"{tgt.id}"
            for interval in [i.interval for i in ifragments]:
                start = interval.start
                stop = interval.stop
                owriter.write_item(
                    seqid,
                    alphabet_name(seq.alphabet),
                    prof.profid,
                    alphabet_name(prof.alphabet),
                    start,
                    stop,
                    prof.window_length,
                )

            if odebug is not os.devnull:
                for i in search_results.debug_table():
                    dwriter.write_row(seqid, i)

    owriter.close()
    odebug.close_intelligently()
