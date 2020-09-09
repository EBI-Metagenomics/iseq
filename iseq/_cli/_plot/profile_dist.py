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
@click.option("--width", help="Figure width.", default=None, type=float)
@click.option("--height", help="Figure height.", default=None, type=float)
@click.option(
    "--quiet/--no-quiet",
    "-q/-nq",
    help="Disable verbosity.",
    default=False,
)
def profile_dist(gff_file, output, width, height, quiet: bool):
    """
    Plot profiles distribution.
    """
    import seaborn as sns
    from matplotlib import pyplot as plt

    gff = read_gff(gff_file, verbose=not quiet)
    df = gff.to_dataframe()

    sns.set(color_codes=True)
    figsize = list(plt.rcParams.get("figure.figsize"))
    if width is not None:
        figsize[0] = float(width)
    if height is not None:
        figsize[1] = float(height)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title("Distribution of Profiles")
    df["Profile"] = df["att_Profile_name"]
    prof_size = df.groupby("Profile").size()
    prof_size.plot.bar(axes=ax)
    ax.set_ylabel("Number of subsequences")
    fig.tight_layout()

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
