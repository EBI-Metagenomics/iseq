import sys

import click

from iseq.gff import read as read_gff

# iseq plot e-value-hist --input tmp/output.gff
from .e_value_dist import e_value_dist


@click.group(name="plot", context_settings=dict(help_option_names=["-h", "--help"]))
def plot():
    """
    Plot.
    """


plot.add_command(e_value_dist)
