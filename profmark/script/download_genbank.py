from pathlib import Path
from iseq.profmark import download_genbank

acc_filepath = Path(snakemake.input[0])
assert acc_filepath.exists()

accession = acc_filepath.name

assert len(snakemake.output) == 2

folder = Path(snakemake.output[0]).parent

download_genbank(folder, accession)
