import click

from ._misc import get_version
from .search import search
from .score import score


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(get_version())
def cli():
    """
    Find nucleotide sequences against protein profiles.
    """


cli.add_command(search)
cli.add_command(score)
