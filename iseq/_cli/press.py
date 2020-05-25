import os
from typing import IO

import click
from tqdm import tqdm

from hmmer_reader import open_hmmer
from nmm import CanonicalAminoAlphabet, DNAAlphabet, GeneticCode, Model, Output

from .._alphabet import infer_hmmer_alphabet
from ..frame import create_frame_profile
from ..standard import create_standard_profile


@click.command()
@click.argument("profile", type=click.File("r"))
@click.option(
    "--epsilon", type=float, default=1e-2, help="Indel probability. Defaults to 1e-2."
)
@click.option(
    "--output",
    type=click.File("w"),
    help="Save results to OUTPUT (NMM format).",
    default=None,
)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
def press(
    profile, epsilon: float, output, quiet,
):
    """
    Press.

    TODO: fix this doc.
    Search nucleotide sequence(s) against a protein profiles database.

    An OUTPUT line determines an association between a TARGET subsequence and
    a PROFILE protein profile. An association maps a target subsequence to a
    profile and represents a potential homology. Expect many false positive
    associations as we are not filtering out by statistical significance.
    """
    if quiet:
        click.open_file(os.devnull, "a")
    else:
        click.get_text_stream("stdout")

    if output is None:
        output = (profile.name + ".nmm").encode()

    base_abc = DNAAlphabet()
    # gcode = GeneticCode(base_abc, CanonicalAminoAlphabet())

    with Output.create(output) as file:
        for hmmprof in tqdm(open_hmmer(profile)):
            # scanner.show_profile_parser(hmmprof)
            # scanner.process_profile(hmmprof, targets)
            prof = create_frame_profile(hmmprof, base_abc, epsilon)
            hmm = prof.alt_model._hmm
            dp = hmm.create_dp(prof.alt_model.special_node.T)
            model = Model.create(hmm, dp)
            file.write(model)


# def _infer_profile_alphabet(profile: IO[str]):
#     hmmer = open_hmmer(profile)
#     hmmer_alphabet = infer_hmmer_alphabet(hmmer)
#     profile.seek(0)
#     if hmmer_alphabet is None:
#         raise click.UsageError("Could not infer alphabet from PROFILE.")
#     return hmmer_alphabet
