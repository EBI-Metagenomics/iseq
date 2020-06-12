import os
from typing import IO, Optional

import click
from fasta_reader import FASTAWriter, open_fasta
from hmmer_reader import open_hmmer
from nmm import AminoAlphabet, BaseAlphabet, CanonicalAminoAlphabet, GeneticCode

from iseq.alphabet import infer_fasta_alphabet, infer_hmmer_alphabet


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option(
    "--epsilon", type=float, default=1e-2, help="Indel probability. Defaults to 1e-2."
)
@click.option(
    "--output",
    type=click.File("w"),
    help="Save results to OUTPUT (GFF format).",
    default=os.devnull,
)
@click.option(
    "--ocodon",
    type=click.File("w"),
    help="Save codon sequences to OCODON (FASTA format).",
    default=os.devnull,
)
@click.option(
    "--oamino",
    type=click.File("w"),
    help="Save amino acid sequences to OAMINO (FASTA format).",
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
def pscan(
    profile, target, epsilon: float, output, ocodon, oamino, quiet, window: int,
):
    """
    Search nucleotide sequence(s) against a protein profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """
    from .scanner import OutputWriter
    from .protein_scanner import ProteinScanner

    owriter = OutputWriter(output, window)
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)

    if quiet:
        stdout = click.open_file(os.devnull, "a")
    else:
        stdout = click.get_text_stream("stdout")

    profile_abc = _infer_profile_alphabet(profile)
    target_abc = _infer_target_alphabet(target)

    scanner: Optional[ProteinScanner] = None

    assert isinstance(target_abc, BaseAlphabet) and isinstance(
        profile_abc, AminoAlphabet
    )

    gcode = GeneticCode(target_abc, CanonicalAminoAlphabet())

    scanner = ProteinScanner(owriter, cwriter, awriter, gcode, epsilon, window, stdout)

    with open_fasta(target) as fasta:
        targets = list(fasta)

    for hmmprof in open_hmmer(profile):
        scanner.show_profile_parser(hmmprof)
        scanner.process_profile(hmmprof, targets)

    scanner.finalize_stream("output", output)
    scanner.finalize_stream("ocodon", ocodon)
    scanner.finalize_stream("oamino", oamino)


def _infer_profile_alphabet(profile: IO[str]):
    hmmer = open_hmmer(profile)
    hmmer_alphabet = infer_hmmer_alphabet(hmmer)
    profile.seek(0)
    if hmmer_alphabet is None:
        raise click.UsageError("Could not infer alphabet from PROFILE.")
    return hmmer_alphabet


def _infer_target_alphabet(target: IO[str]):
    fasta = open_fasta(target)
    target_alphabet = infer_fasta_alphabet(fasta)
    target.seek(0)
    if target_alphabet is None:
        raise click.UsageError("Could not infer alphabet from TARGET.")
    return target_alphabet
