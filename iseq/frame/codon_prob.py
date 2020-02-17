from typing import Dict
from math import log


from nmm import GeneticCode, AminoTable, LPROB_ZERO, Codon, CodonProb, CBaseAlphabet


class FrameCodonProb(CodonProb):
    def __init__(
        self, base_abc: CBaseAlphabet, aminot: AminoTable, gencode: GeneticCode
    ):
        from numpy import logaddexp

        super().__init__(base_abc)
        self._lprob: Dict[Codon, float] = {}

        codon_lprobs = []
        lprob_norm = LPROB_ZERO
        for i in range(len(aminot.alphabet.symbols)):
            aa = aminot.alphabet.symbols[i : i + 1]
            lprob = aminot.lprob(aa)

            codons = gencode.codons(aa)
            if len(codons) == 0:
                continue

            norm = log(len(codons))
            for codon in codons:
                codon_lprobs.append((codon, lprob - norm))
                lprob_norm = logaddexp(lprob_norm, codon_lprobs[-1][1])

        for codon, lprob in codon_lprobs:
            self.set_lprob(codon, lprob - lprob_norm)
