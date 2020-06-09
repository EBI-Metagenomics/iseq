import os

import click
from hmmer_reader import open_hmmer
from nmm import DNAAlphabet, Model, Output

from iseq.protein import create_profile


@click.command()
@click.argument("profile", type=click.File("r"))
@click.option(
    "--epsilon", type=float, default=1e-2, help="Indel probability. Defaults to 1e-2."
)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
def press(
    profile, epsilon: float, quiet,
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
    from tqdm import tqdm

    if quiet:
        click.open_file(os.devnull, "a")
    else:
        click.get_text_stream("stdout")

    alt_filepath = (profile.name + ".alt").encode()
    null_filepath = (profile.name + ".null").encode()
    meta_filepath = (profile.name + ".meta").encode()

    base_abc = DNAAlphabet()

    with Output.create(alt_filepath) as afile:
        with Output.create(null_filepath) as nfile:
            with open(meta_filepath, "w") as mfile:
                for hmmprof in tqdm(open_hmmer(profile)):
                    prof = create_profile(hmmprof, base_abc, epsilon)

                    hmm = prof.alt_model.hmm
                    dp = hmm.create_dp(prof.alt_model.special_node.T)
                    model = Model.create(hmm, dp)
                    afile.write(model)

                    hmm = prof.null_model.hmm
                    dp = hmm.create_dp(prof.null_model.state)
                    model = Model.create(hmm, dp)
                    nfile.write(model)

                    mfile.write(dict(hmmprof.metadata)["ACC"] + "\n")
