import click

from ._misc import get_version
from .scan import scan
from .score import score
from .gff import gff_dedup


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(get_version())
def cli():
    """
    Find nucleotide sequences against protein profiles.
    """


cli.add_command(scan)
cli.add_command(score)
cli.add_command(gff_dedup)
