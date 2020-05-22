from typing import Optional, Union

from nmm import CanonicalAminoAlphabet, DNAAlphabet, RNAAlphabet

Alphabets = Union[DNAAlphabet, RNAAlphabet, CanonicalAminoAlphabet]


def infer_alphabet(sequence: bytes) -> Optional[Alphabets]:
    dna = DNAAlphabet()
    rna = RNAAlphabet()
    amino = CanonicalAminoAlphabet()

    abc = set(sequence)

    if len(abc - set(dna.symbols)) == 0:
        return dna

    if len(abc - set(rna.symbols)) == 0:
        return rna

    if len(abc - set(amino.symbols)) == 0:
        return amino

    return None
