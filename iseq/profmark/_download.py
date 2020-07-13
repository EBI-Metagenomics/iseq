from pathlib import Path

from Bio import Entrez
from iseq.file import cleanup_invalid_filepath, assert_file_hash
from ._accession import accession_hash


__all__ = ["download_genbank"]


def download_genbank(folder: Path, accession: str):

    filepath = folder / f"{accession}.gbk"
    cleanup_invalid_filepath(filepath, accession_hash(accession)["genbank"])
    if not filepath.exists():
        download_efetch(filepath, accession, "gb")
    assert_file_hash(filepath, accession_hash(accession)["genbank"])

    filepath = folder / f"{accession}.fasta"
    cleanup_invalid_filepath(filepath, accession_hash(accession)["fasta"])
    if not filepath.exists():
        download_efetch(filepath, accession, "fasta")
    assert_file_hash(filepath, accession_hash(accession)["fasta"])


def download_efetch(filepath: Path, accession: str, rettype: str):
    Entrez.email = "horta@ebi.ac.uk"
    efetch = Entrez.efetch

    with efetch(db="nuccore", id=accession, rettype=rettype, retmode="text") as handle:
        with open(filepath, "w") as file:
            file.write(handle.read())
