from typing import Dict

import importlib_resources as pkg_resources
from nmm import Codon, RNAAlphabet

from .codon_table import CodonTable

__all__ = ["CodonUsage"]


class CodonUsage:
    def __init__(self, codon_table: CodonTable, taxonomy_id: str):
        from pandas import read_pickle

        import iseq

        buffer = pkg_resources.open_binary(iseq, "codon_usage_20200706.pkl.xz")

        df = read_pickle(buffer, compression="xz")
        df = df.loc[taxonomy_id]

        if len(df) == 0:
            raise ValueError(f"Unkown taxonomy ID {taxonomy_id}.")

        if isinstance(codon_table.base_alphabet, RNAAlphabet):
            df["Codon"] = df["Codon"].str.replace("T", "U")

        df["Prob"] = df["Number"] / df["Number"].sum()

        self._codon_table = codon_table
        self._map: Dict[Codon, float] = {}

        for _, row in df.iterrows():
            codon = Codon.create(row["Codon"].encode(), codon_table.base_alphabet)
            self._map[codon] = row["Prob"]

    def codon_prob(self, codon: Codon) -> float:
        return self._map[codon]

    def amino_acid_prob(self, amino_acid: bytes) -> float:
        return sum(self.codon_prob(c) for c in self._codon_table.codons(amino_acid))

    @property
    def start_prob(self) -> float:
        return sum(self.codon_prob(c) for c in self._codon_table.start_codons)

    @property
    def stop_prob(self) -> float:
        return sum(self.codon_prob(c) for c in self._codon_table.stop_codons)

    @property
    def codon_table(self) -> CodonTable:
        return self._codon_table
