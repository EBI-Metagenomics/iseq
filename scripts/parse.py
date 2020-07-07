from Bio import SeqIO

# from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import DNAAlphabet, RNAAlphabet

filename = "AE014075.1.gb"
for seq_record in SeqIO.parse(filename, "genbank"):
    for feature in seq_record.features:
        if feature.type != "CDS":
            continue

        print(feature)
        breakpoint()
        continue

        seq = feature.extract(seq_record).seq
        if isinstance(seq.alphabet, DNAAlphabet):
            remains = len(set(str(seq)) - set(IUPAC.unambiguous_dna.letters))
            if remains > 0:
                continue

        elif isinstance(seq.alphabet, RNAAlphabet):
            remains = len(set(str(seq)) - set(IUPAC.unambiguous_rna.letters))
            if remains > 0:
                continue

        else:
            raise ValueError("Unkown alphabet.")

        msg = ""
        if "function" in feature.qualifiers:
            msg += str(feature.qualifiers["function"])
        if "protein_id" in feature.qualifiers:
            msg += str(feature.qualifiers["protein_id"])
        print(msg)
    # print(seq_record.id)
    # print(repr(seq_record.seq))
    # print(len(seq_record))

