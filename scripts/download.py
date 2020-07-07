from Bio import Entrez

Entrez.email = "horta@ebi.ac.uk"

accession = "AE014075.1"
# accession = "AE005673.1"
with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
    with open(f"{accession}.gb", "w") as file:
        file.write(handle.read())
