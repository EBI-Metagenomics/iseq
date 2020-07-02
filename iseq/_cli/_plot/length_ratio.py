import click

from iseq.gff import read as read_gff


@click.command()
@click.argument("gff_file", type=click.File("r"))
@click.argument(
    "hmmer_file", type=click.Path(exists=True, dir_okay=False, resolve_path=True)
)
@click.option(
    "--output",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save figure to OUTPUT.",
    default=None,
)
@click.option("--width", help="Figure width.", default=None, type=float)
@click.option("--height", help="Figure height.", default=None, type=float)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable verbosity.", default=False,
)
def length_ratio(gff_file, hmmer_file, output, width, height, quiet: bool):
    """
    Plot length-ratio distribution.
    """
    import seaborn as sns
    from matplotlib import pyplot as plt
    from hmmer_reader import fetch_metadata

    gff = read_gff(gff_file, verbose=not quiet)
    df = gff.to_dataframe()

    meta = fetch_metadata(hmmer_file)
    meta = meta.set_index(["NAME", "ACC"])
    idx = list(zip(df["att_Profile_name"].tolist(), df["att_Profile_acc"].tolist()))
    df["Profile_length"] = meta.loc[idx, "LENG"].values
    breakpoint()

    # profile_length = pd.Series(profile_length["length"].values, index=profile_length["name"].values)
    # gff_hits["Length"] = gff_hits["end"] - gff_hits["start"] + 1
    # for idx in gff_hits.index:
    #     gff_hits.loc[idx, "Profile_length"] = profile_length[gff_hits.loc[idx, "Profile_name"]]
    # gff_hits["Length_ratio"] = gff_hits["Length"] / (3 * gff_hits["Profile_length"])

    sns.set(color_codes=True)
    figsize = list(plt.rcParams.get("figure.figsize"))
    if width is not None:
        figsize[0] = float(width)
    if height is not None:
        figsize[1] = float(height)

    pass

    # fig, ax = plt.subplots(figsize=figsize)
    # ax.set_title("Distribution of Profiles")
    # df["Profile"] = df["att_Profile_name"]
    # prof_size = df.groupby("Profile").size()
    # prof_size.plot.bar(axes=ax)
    # ax.set_ylabel("Number of subsequences")
    # fig.tight_layout()

    # if not quiet:
    #     click.echo("Plotting... ", nl=False)
    # if output is None:
    #     plt.show()
    # else:
    #     plt.savefig(output)

    # if not quiet:
    #     click.echo("done.", nl=True)
    #     if output is not None:
    #         click.echo(f"Saved to {output}.")
