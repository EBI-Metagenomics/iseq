from typing import Dict
from math import log


from nmm import (
    GeneticCode,
    AminoTable,
    LPROB_ZERO,
    Codon,
    BaseTable,
    CBaseAlphabet,
    codon_iter,
    LPROB_ZERO,
)
from .codon_prob import CodonProb


class FrameBaseTable(BaseTable):
    def __init__(self, codon_prob: CodonProb):
        from numpy import logaddexp

        base_abc = codon_prob.base
        base_lprob = {base: LPROB_ZERO for base in base_abc.symbols}
        norm = log(3)
        for codon in codon_iter(base_abc):
            lprob = codon_prob.get_lprob(codon)
            triplet = codon.symbols

            base_lprob[triplet[0]] = logaddexp(base_lprob[triplet[0]], lprob - norm)
            base_lprob[triplet[1]] = logaddexp(base_lprob[triplet[1]], lprob - norm)
            base_lprob[triplet[2]] = logaddexp(base_lprob[triplet[2]], lprob - norm)

        super().__init__(base_abc, [base_lprob[base] for base in base_abc.symbols])

        # lprobs: Dict[BaseAlphabet, list] = {
        #     BaseAlphabet(sym): [] for sym in alphabet.symbols
        # }
        # for codon, lprob in codon_lprobs.items():
        #     lprobs[codon.base(0)] += [lprob - lprob_norm]
        #     lprobs[codon.base(1)] += [lprob - lprob_norm]
        #     lprobs[codon.base(2)] += [lprob - lprob_norm]

        # return {b: logsumexp(lp) for b, lp in lprobs.items()}

        # from numpy import logaddexp

        # super().__init__(base_abc)
        # self._lprob: Dict[Codon, float] = {}

        # codon_lprobs = []
        # lprob_norm = LPROB_ZERO
        # for i in range(len(aminot.alphabet.symbols)):
        #     aa = aminot.alphabet.symbols[i : i + 1]
        #     lprob = aminot.lprob(aa)

        #     codons = gencode.codons(aa)
        #     if len(codons) == 0:
        #         continue

        #     norm = log(len(codons))
        #     for codon in codons:
        #         codon_lprobs.append((codon, lprob - norm))
        #         lprob_norm = logaddexp(lprob_norm, codon_lprobs[-1][1])

        # for codon, lprob in codon_lprobs:
        #     self.set_lprob(codon, lprob - lprob_norm)
