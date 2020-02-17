from math import log


from nmm import (
    BaseTable,
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
