import click

from .gff import gff_dedup
from .press import press
from .scan import scan
from .score import score
from .unpress import unpress


def get_version():
    import re
    import iseq
    import importlib_resources as pkg_resources

    content = pkg_resources.read_text(iseq, "__init__.py")
    c = re.compile(r"__version__ *= *('[^']+'|\"[^\"]+\")")
    m = c.search(content)
    return m.groups()[0][1:-1]


@click.group(name="iseq", context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(get_version())
def cli():
    """
    Find nucleotide sequences against protein profiles.
    """


cli.add_command(scan)
cli.add_command(score)
cli.add_command(press)
cli.add_command(unpress)
cli.add_command(gff_dedup)
