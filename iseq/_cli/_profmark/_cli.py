import click

from .roc import plot_roc


@click.group(name="profmark", context_settings=dict(help_option_names=["-h", "--help"]))
def profmark():
    """
    Profmark.
    """


profmark.add_command(plot_roc)
