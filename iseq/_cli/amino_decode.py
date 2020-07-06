import click
from fasta_reader import FASTAWriter, open_fasta


@click.command()
@click.argument("fasta", type=click.File("r"))
@click.option(
    "--output",
    type=click.File("w", lazy=True),
    help="Save results to OUTPUT (FASTA format).",
    default="-",
)
def amino_decode(fasta, output):
    owriter = FASTAWriter(output, 60)

    for target in open_fasta(fasta):
        owriter.write_item(target.defline, target.sequence)

    output.close_intelligently()
