from typing import Union

from nmm.alphabet import DNAAlphabet, RNAAlphabet


def infer_alphabet(sequence: bytes) -> Union[DNAAlphabet, RNAAlphabet, None]:
    dna = DNAAlphabet()
    rna = RNAAlphabet()

    abc = sorted(list(set(sequence)))

    if sorted(list(dna.symbols)) == abc:
        return dna

    if sorted(list(rna.symbols)) == abc:
        return rna

    return None
