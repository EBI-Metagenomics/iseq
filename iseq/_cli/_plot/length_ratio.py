import click

from iseq.gff import read as read_gff


@click.command()
@click.argument("gff_file", type=click.File("r"))
@click.argument(
    "hmmer_file", type=click.Path(exists=True, dir_okay=False, resolve_path=True)
)
@click.argument("prof_leng_mult", type=int)
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
def length_ratio(
    gff_file, hmmer_file, prof_leng_mult, output, width, height, quiet: bool
):
    """
    Plot length-ratio distribution.
    """
    import seaborn as sns
    from hmmer_reader import fetch_metadata
    from matplotlib import pyplot as plt

    gff = read_gff(gff_file, verbose=not quiet)
    df = gff.to_dataframe()

    meta = fetch_metadata(hmmer_file)
    meta = meta.set_index(["NAME", "ACC"])
    idx = list(zip(df["att_Profile_name"].tolist(), df["att_Profile_acc"].tolist()))
    df["Profile_length"] = meta.loc[idx, "LENG"].values * prof_leng_mult

    df["Length"] = df["end"] - df["start"] + 1
    df["Length_ratio"] = df["Length"] / df["Profile_length"]

    sns.set(color_codes=True)
    figsize = list(plt.rcParams.get("figure.figsize"))
    if width is not None:
        figsize[0] = float(width)
    if height is not None:
        figsize[1] = float(height)

    fig, ax = plt.subplots(figsize=figsize)
    sns.distplot(df["Length_ratio"], ax=ax)
    ax.set_title("Subsequence length / Profile length")
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
