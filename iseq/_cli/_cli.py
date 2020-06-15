import click

from .._version import __version__
from .bscan import bscan
from .hscan import hscan
from .press import press
from .pscan import pscan


@click.group(name="iseq", context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
def cli():
    """
    Find nucleotide sequences against protein profiles.
    """


cli.add_command(bscan)
cli.add_command(hscan)
cli.add_command(press)
cli.add_command(pscan)
