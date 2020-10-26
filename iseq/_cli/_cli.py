import click

from ._plot import plot
from .amino_decode import amino_decode
from .bscan import bscan
from .gff_filter import gff_filter
from .hscan import hscan
from .press import press
from .pscan import pscan
from .pscan2 import pscan2
from .pscan3 import pscan3


@click.group(name="iseq", context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option()
def cli():
    """
    Find nucleotide sequences against protein profiles.
    """


cli.add_command(amino_decode)
cli.add_command(bscan)
cli.add_command(gff_filter)
cli.add_command(hscan)
cli.add_command(plot)
cli.add_command(press)
cli.add_command(pscan)
cli.add_command(pscan2)
cli.add_command(pscan3)
