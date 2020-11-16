import click
import nmm
from fasta_reader import FASTAWriter, read_fasta

from iseq.codon_table import CodonTable
from iseq.gencode import GeneticCode


@click.command()
@click.argument("fasta", type=click.File("r"))
@click.option(
    "--output",
    type=click.File("w", lazy=True),
    help="Save results to OUTPUT (FASTA format).",
    default="-",
)
@click.option(
    "--seed",
    type=int,
    help="Random seed used to initialize the pseudo-random number generator.",
    default=0,
)
@click.option(
    "--transl-table",
    type=int,
    help="Translation table number. Defaults to 1 (The Standard Code).",
    default=1,
)
@click.option(
    "--alph",
    type=click.Choice(["dna", "rna"]),
    help="Nucleotide alphabet.",
    default="dna",
)
def amino_decode(fasta, output, seed: int, transl_table: int, alph: str):
    from numpy.random import RandomState

    owriter = FASTAWriter(output, 60)

    if alph == "dna":
        base_abc = nmm.DNAAlphabet()
    else:
        base_abc = nmm.RNAAlphabet()

    amino_abc = nmm.IUPACAminoAlphabet()
    codon_table = CodonTable(base_abc, amino_abc, GeneticCode(id=transl_table))

    random = RandomState(seed)
    for target in read_fasta(fasta):
        amino_seq = target.sequence
        nucls = []
        for amino in amino_seq:
            codons = codon_table.codons(amino.encode())
            codon = random.choice(codons)
            nucls.append(codon.symbols.decode())

        owriter.write_item(target.defline, "".join(nucls))

    output.close_intelligently()
