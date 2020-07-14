__all__ = ["accession_hash"]


def accession_hash(accession: str):
    accessions = {
        "AE014075.1": {
            "amino": "d14d16c7455caeb1411101bf54967643adfe7619a12b8816e3c755613f337f38",
            "nucl": "96308d9d2da397993f596f28559935ac0c1f98ab817fc92627e2217826ba1cef",
            "genbank": "c26ac7e4208de948a0d3644423b48d9eaab3f4a3bf0e3ce1b64af6ded8e2d7af",
            "fasta": "d7e170b7758aa367b45f99f1c87e31977c1f1309b687a9e36d4402b3ccf65d82",
        }
    }

    return accessions[accession]
