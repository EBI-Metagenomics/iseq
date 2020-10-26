import warnings
from typing import Dict, List, Optional, Union

from nmm import AminoAlphabet, Codon, DNAAlphabet, RNAAlphabet

from .gencode import GeneticCode

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=PendingDeprecationWarning)
    from Bio.Data.CodonTable import (
        NCBICodonTable,
        unambiguous_dna_by_id,
        unambiguous_dna_by_name,
    )

__all__ = ["CodonTable"]


class CodonTable:
    """
    Codon table.

    Parameters
    ----------
    base_abc
        Base alphabet.
    amino_abc
        Amino acid alphabet.
    gencode
        NCBI `genetic code`_. Defaults to `GeneticCode("Standard")`.

    .. _genetic code: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """

    def __init__(
        self,
        base_abc: Union[DNAAlphabet, RNAAlphabet],
        amino_abc: AminoAlphabet,
        gencode: Optional[GeneticCode] = None,
    ):

        self._base_alphabet = base_abc
        self._amino_alphabet = amino_abc

        if gencode is None:
            gencode = GeneticCode("Standard")
        table = translation_table(gencode)

        def replace(seq: bytes):
            if isinstance(base_abc, RNAAlphabet):
                seq = seq.replace(b"T", b"U")
            return seq

        self._codons: Dict[bytes, List[Codon]] = {}

        self._start_codons = []
        for t, a in table.forward_table.items():
            triplet = replace(t.encode())
            aa = a.encode()
            if aa not in self._codons:
                self._codons[aa] = []

            codon = Codon.create(triplet, base_abc)
            if t in table.start_codons:
                self._start_codons.append(codon)
            self._codons[aa].append(codon)

        self._amino_acid: Dict[Codon, bytes] = {}
        for aa, codons in self._codons.items():
            for codon in codons:
                self._amino_acid[codon] = aa

        self._stop_codons = []
        for t in table.stop_codons:
            triplet = replace(t.encode())

            self._stop_codons.append(Codon.create(triplet, base_abc))

        assert len(self._amino_acid) <= 64
        assert len(self._amino_acid) + len(self._stop_codons) == 64

    @property
    def start_codons(self) -> List[Codon]:
        return self._start_codons

    @property
    def stop_codons(self) -> List[Codon]:
        return self._stop_codons

    def codons(self, amino_acid: Optional[bytes] = None) -> List[Codon]:

        if amino_acid is None:
            all_codons: List[Codon] = []
            for aa in self.amino_acids:
                all_codons += self.codons(aa)
            all_codons += self.stop_codons
            return list(sorted(all_codons))

        amino_acid = amino_acid.upper()
        return self._codons.get(amino_acid, [])

    def amino_acid(self, codon: Codon) -> bytes:
        return self._amino_acid[codon]

    def decode(self, codon: Codon) -> bytes:
        if codon in self._amino_acid:
            return self._amino_acid[codon]
        assert codon in self._stop_codons
        return self._amino_alphabet.stop_symbol

    @property
    def amino_acids(self) -> List[bytes]:
        return list(set(self._amino_acid.values()))

    @property
    def base_alphabet(self) -> Union[DNAAlphabet, RNAAlphabet]:
        return self._base_alphabet

    @property
    def amino_alphabet(self) -> AminoAlphabet:
        return self._amino_alphabet


def translation_table(gencode: GeneticCode) -> NCBICodonTable:

    if gencode.id not in unambiguous_dna_by_id:
        names = str(list(unambiguous_dna_by_name.keys()))
        msg = f"Unknown translation table {gencode.name}. Possible names are: {names}."
        raise ValueError(msg)

    return unambiguous_dna_by_id[gencode.id]
