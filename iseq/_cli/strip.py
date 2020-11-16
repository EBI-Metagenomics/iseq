# It might be useful in the future.
# import click
# from fasta_reader import FASTAWriter, read_fasta

# from iseq.gff import read as read_gff


# @click.command()
# @click.argument("gff", type=click.Path(exists=True, dir_okay=False))
# @click.argument("target", type=click.Path(exists=True, dir_okay=False))
# @click.option(
#     "--gff-output",
#     type=click.Path(exists=False, dir_okay=False),
#     help="Save the surviving GFF items to GFF_OUTPUT.",
# )
# @click.option(
#     "--target-output",
#     type=click.Path(exists=False, dir_okay=False),
#     help="Save the surviving sequences to TARGET_OUTPUT.",
# )
# @click.option(
#     "--e-value", type=float, default=1e-60, help="E-value threshold. Defaults to 1e-60."
# )
# def strip(gff, target, gff_output, target_output, e_value: float):
#     """
#     Filter-out sequences having E-value above a given threshold.
#     """

#     if gff_output is None:
#         gff_output = gff + ".strip"

#     if target_output is None:
#         target_output = target + ".strip"

#     ok_lines = set()
#     ok_items = set()
#     for i, gff_item in enumerate(read_gff(gff).items()):
#         fields = create_fields(gff_item.attributes)
#         if "E-value" in fields:
#             if float(fields["E-value"]) <= e_value:
#                 ok_lines.add(i)
#                 ok_items.add(fields["ID"])

#     with open(gff, "r") as r:
#         with open(gff_output, "w") as w:
#             w.write(r.readline())
#             for i, line in enumerate(r):
#                 if i in ok_lines:
#                     w.write(line)

#     with FASTAWriter(target_output) as writer:
#         with read_fasta(target) as fasta:
#             for tgt in fasta:
#                 if tgt.defline in ok_items:
#                     writer.write_item(tgt.defline, tgt.sequence)

#     print(f"Saved GFF to {gff_output}.")
#     print(f"Saved targets to {target_output}.")


# def create_fields(att: str):
#     fields = {}
#     for v in att.split(";"):
#         name, value = v.split("=", 1)
#         fields[name] = value
#     return fields
