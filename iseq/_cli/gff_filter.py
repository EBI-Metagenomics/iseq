import sys

import click

from iseq.gff import read as read_gff


@click.command()
@click.argument("gff_file", type=click.File("r"))
@click.option(
    "--max-e-value",
    type=float,
    default=None,
    help="Filter out items for which E-value > --max-e-evalue.",
)
def gff_filter(gff_file, max_e_value):
    """
    Filter out items from a GFF_FILE file.

    The resulting file will be written to the standard output.
    """

    gff = read_gff(gff_file)

    kwargs = {}
    if max_e_value is not None:
        kwargs["max_e_value"] = float(max_e_value)

    gff = gff.filter(**kwargs)

    gff.write_file(sys.stdout)
