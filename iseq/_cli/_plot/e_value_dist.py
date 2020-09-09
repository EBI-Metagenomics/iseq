import click

from iseq.gff import read as read_gff


@click.command()
@click.argument("gff_file", type=click.File("r"))
@click.option(
    "--output",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save figure to OUTPUT.",
    default=None,
)
@click.option(
    "--quiet/--no-quiet",
    "-q/-nq",
    help="Disable verbosity.",
    default=False,
)
def e_value_dist(gff_file, output, quiet: bool):
    """
    Plot E-values distribution.
    """
    import seaborn as sns
    from matplotlib import pyplot as plt
    from numpy import finfo, log10

    gff = read_gff(gff_file, verbose=not quiet)
    df = gff.to_dataframe()

    if "att_E-value" not in df.columns:
        click.echo("No E-value found.")
        return

    df = df[~df["att_E-value"].isna()]

    if len(df) == 0:
        click.echo("Only NaN E-values.")
        return

    min_possible = 2 ** finfo(float).minexp
    df.loc[df["att_E-value"] == 0.0, "att_E-value"] = min_possible
    df["-log10(E-value)"] = -log10(df["att_E-value"])

    sns.set(color_codes=True)

    ax = df["-log10(E-value)"].plot(kind="hist", logx=False, logy=True)
    ax.set_title("Distribution of E-values")
    ax.set_xlabel("E-value (-log10(x))")

    if not quiet:
        click.echo("Plotting... ", nl=False)
    if output is None:
        plt.show()
    else:
        plt.savefig(output)

    if not quiet:
        click.echo("done.", nl=True)
        if output is not None:
            click.echo(f"Saved to {output}.")
