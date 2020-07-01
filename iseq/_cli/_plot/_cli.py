import click

from .e_value_dist import e_value_dist
from .profile_dist import profile_dist


@click.group(name="plot", context_settings=dict(help_option_names=["-h", "--help"]))
def plot():
    """
    Plot.
    """


plot.add_command(e_value_dist)
plot.add_command(profile_dist)
