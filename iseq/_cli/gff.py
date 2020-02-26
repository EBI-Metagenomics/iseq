import click


@click.command()
@click.argument("gff", type=click.File("rw"))
def gff_dedup(gff):
    """
    Deduplicate a GFF file.
    """
    pass
