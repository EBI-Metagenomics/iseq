from pathlib import Path

from Bio import Entrez
from iseq.file import cleanup_invalid_filepath, assert_file_hash


__all__ = ["download_genbank"]

_accessions = {
    "AE014075.1": {
        "genbank": "c26ac7e4208de948a0d3644423b48d9eaab3f4a3bf0e3ce1b64af6ded8e2d7af",
        "fasta": "d7e170b7758aa367b45f99f1c87e31977c1f1309b687a9e36d4402b3ccf65d82",
    }
}


def download_genbank(folder: Path, accession: str):

    filepath = folder / f"{accession}.gbk"
    cleanup_invalid_filepath(filepath, _accessions[accession]["genbank"])
    if not filepath.exists():
        download_efetch(filepath, accession, "gb")
    assert_file_hash(filepath, _accessions[accession]["genbank"])

    filepath = folder / f"{accession}.fasta"
    cleanup_invalid_filepath(filepath, _accessions[accession]["fasta"])
    if not filepath.exists():
        download_efetch(filepath, accession, "fasta")
    assert_file_hash(filepath, _accessions[accession]["fasta"])


def download_efetch(filepath: Path, accession: str, rettype: str):
    Entrez.email = "horta@ebi.ac.uk"
    efetch = Entrez.efetch

    with efetch(db="nuccore", id=accession, rettype=rettype, retmode="text") as handle:
        with open(filepath, "w") as file:
            file.write(handle.read())
