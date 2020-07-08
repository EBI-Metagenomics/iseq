from Bio import Entrez

Entrez.email = "horta@ebi.ac.uk"
efetch = Entrez.efetch

accession = "AE014075.1"

# with efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
#     filename = f"{accession}.gbk"
#     with open(filename, "w") as file:
#         file.write(handle.read())

with efetch(db="nuccore", id=accession, rettype="fasta", retmode="text") as handle:
    filename = f"{accession}.fasta"
    with open(filename, "w") as file:
        file.write(handle.read())
