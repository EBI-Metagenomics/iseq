import os

import click

from fasta_reader import open_fasta
from imm import Sequence
from nmm import DNAAlphabet, Input, GeneticCode, CanonicalAminoAlphabet
from .. import wrap
from ..frame._typing import FrameAltModel, FrameNullModel
from ..frame._profile import FrameProfile
from .scan import _infer_target_alphabet


@click.command()
@click.argument("alt_filepath", type=str)
@click.argument("null_filepath", type=str)
@click.argument("target", type=click.File("r"))
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
@click.option(
    "--hmmer3-compat/--no-hmmer3-compat",
    help="Enable full HMMER3 compatibility.",
    default=False,
)
def unpress(alt_filepath, null_filepath, target, quiet, hmmer3_compat: bool):
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

    # DNAAlphabet()
    # gcode = GeneticCode(base_abc, CanonicalAminoAlphabet())

    # Unpress
    # [-18.5960979922376]
    # [-19.328720503452924]
    # [-19.23670894701474]

    # Scan
    # [-18.227717053999672]

    # Unpress
    # [-18.55713731433567]
    # <Path:<S,0>,<B,0>,<M17,3>,<M18,2>,<E,0>,<T,0>>
    # self._hmmer3_compat = False
    # SpecialTransitions(NN=-0.47000362924573547, NB=-0.980829253011726,
    #         EC=-0.6931471805599453, CC=-0.47000362924573547,
    #         CT=-0.980829253011726, EJ=-0.6931471805599453,
    #         JJ=-0.47000362924573547, JB=-0.980829253011726,
    #         RR=-0.18232155679395468, BM=0.0, ME=0.0)

    # <Sequence:CACGT>
    # [-18.227717053999672]
    # <Path:<S,0>,<B,0>,<M8,3>,<M9,2>,<E,0>,<T,0>>
    # self._hmmer3_compat = False
    # SpecialTransitions(NN=-0.47000362924573547, NB=-0.980829253011726,
    #     EC=-0.6931471805599453, CC=-0.47000362924573547,
    #     CT=-0.980829253011726, EJ=-0.6931471805599453,
    #     JJ=-0.47000362924573547, JB=-0.980829253011726,
    #     RR=-0.18232155679395468, BM=0.0, ME=0.0)
    target_abc = _infer_target_alphabet(target)

    # gcode = GeneticCode(target_abc, CanonicalAminoAlphabet())
    with open_fasta(target) as fasta:
        targets = list(fasta)

    with Input.create(alt_filepath.encode()) as alt_file:
        with Input.create(null_filepath.encode()) as null_file:
            for alt, null in tqdm(zip(alt_file, null_file)):
                special_node = wrap.special_node(alt.hmm)
                core_nodes = wrap.core_nodes(alt.hmm)
                alt_model = FrameAltModel.create2(
                    special_node, core_nodes, alt.hmm, alt.dp
                )
                print(alt_model)

                null_model = FrameNullModel.create2(null.hmm)
                print(null_model)

                abc = alt_model.hmm.alphabet
                prof = FrameProfile.create2(abc, null_model, alt_model, hmmer3_compat)
                seq = Sequence.create(targets[0].sequence.encode(), prof.alphabet)
                breakpoint()
                search_results = prof.search(seq, 0)

                return
