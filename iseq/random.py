from typing import List, Optional

from nmm import Codon

from iseq.codon_usage import CodonUsage

__all__ = ["RandomState"]


class RandomState:
    def __init__(self, seed: Optional[int] = None):
        from numpy.random import RandomState as RndState

        self._random = RndState(seed)

    def codon(
        self, codon_usage: CodonUsage, n: int = 1, include_stop=True
    ) -> List[Codon]:
        codons = codon_usage.codon_table.codons()
        if include_stop:
            p = [codon_usage.codon_prob(c) for c in codons]
        else:
            codons = list(set(codons) - set(codon_usage.codon_table.stop_codons))
            norm = 1 - codon_usage.stop_prob
            p = [codon_usage.codon_prob(c) / norm for c in codons]
        return self._random.choice(codons, n, p=p).tolist()

    def amino_acid(
        self, codon_usage: CodonUsage, n: int = 1, include_stop=True
    ) -> List[bytes]:
        codons = self.codon(codon_usage, n, include_stop)
        return [codon_usage.codon_table.decode(c) for c in codons]
