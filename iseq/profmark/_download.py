from pathlib import Path

from Bio import Entrez

from iseq.file import assert_file_hash, cleanup_invalid_filepath

from ._accession import accession_hash

__all__ = ["download_genbank"]


def download_genbank(accession: str, rettype: str):
    """
    Parameters
    ----------
    accession
        Accession number.
    rettype
        Accepted values are ``"gb"`` and ``"fasta"``.
    """

    filepath = Path(f"{accession}.{rettype}")
    cleanup_invalid_filepath(filepath, accession_hash(accession)[rettype])
    if not filepath.exists():
        _download_efetch(filepath, accession, rettype)
    assert_file_hash(filepath, accession_hash(accession)[rettype])


def _download_efetch(filepath: Path, accession: str, rettype: str):
    Entrez.email = "horta@ebi.ac.uk"
    efetch = Entrez.efetch

    with efetch(db="nuccore", id=accession, rettype=rettype, retmode="text") as handle:
        with open(filepath, "w") as file:
            file.write(handle.read())
