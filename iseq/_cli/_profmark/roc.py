from pathlib import Path

import click

from iseq.profmark import ProfMark


@click.command()
@click.argument(
    "dir",
    type=click.Path(exists=True, dir_okay=True, readable=True, resolve_path=True),
)
@click.argument(
    "accession", type=str,
)
@click.option(
    "--output",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save figure to OUTPUT.",
    default=None,
)
@click.option(
    "--solution-space",
    type=click.Choice(["prof-target", "prof", "target"]),
    help="Solution space.",
    default="prof-target",
)
@click.option(
    "--solution-space-idx/--no-solution-space-idx",
    help="Enable index for solution space.",
    default=True,
)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable verbosity.", default=False,
)
def plot_roc(
    dir: str,
    accession: str,
    output: str,
    solution_space: str,
    solution_space_idx: bool,
    quiet: bool,
):
    """
    Plot profiles distribution.
    """
    # import seaborn as sns
    from matplotlib import pyplot as plt

    folder = Path(dir)

    acc_dir = folder / accession
    if not acc_dir.is_dir():
        raise click.UsageError(f"{acc_dir} must be a directory.")

    hmm_file = acc_dir / "dbspace.hmm"

    target_file = acc_dir / "cds_amino.fasta"
    if not target_file.is_file():
        raise click.UsageError(f"{target_file} must exist and be a file.")

    domtbl_file = acc_dir / "domtblout.txt"
    if not domtbl_file.is_file():
        raise click.UsageError(f"{domtbl_file} must exist and be a file.")

    gff_file = acc_dir / "output.gff"
    if not gff_file.is_file():
        raise click.UsageError(f"{gff_file} must exist and be a file.")

    if not quiet:
        click.echo("Plotting... ", nl=False)

    pm = ProfMark(hmm_file, target_file, domtbl_file, gff_file)
    cm = pm.confusion_matrix(solution_space, solution_space_idx)
    P = cm.P
    N = cm.N
    title = f"{accession} (P/N: {P}/{N})"
    if output is None:
        ax = cm.roc.plot()
        ax.set_title(title)
        plt.show()
    else:
        f, ax = plt.subplots()
        cm.roc.plot(ax=ax)
        ax.set_title(title)
        f.savefig(output)

    if not quiet:
        click.echo("done.", nl=True)
        if output is not None:
            click.echo(f"Saved to {output}.")
