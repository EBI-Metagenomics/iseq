import click

from .gff import gff_dedup
from .press import press
from .scan import scan
from .score import score
from .unpress import unpress
from .._version import __version__


@click.group(name="iseq", context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
def cli():
    """
    Find nucleotide sequences against protein profiles.
    """


cli.add_command(scan)
cli.add_command(score)
cli.add_command(press)
cli.add_command(unpress)
cli.add_command(gff_dedup)
