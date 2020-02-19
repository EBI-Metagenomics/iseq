from math import log
from typing import Iterator, List, Sequence, Tuple

from nmm import GeneticCode
from nmm.codon import Codon
from nmm.sequence import SequenceABC
from nmm.state import CodonState, MuteState, NormalState

from ..amino import AminoFragment, AminoPath
from ..fragment import Fragment
from .path import CodonPath
from .step import CodonStep


class CodonFragment(Fragment):
    def __init__(
        self, codons: Sequence[Codon], path: CodonPath, homologous: bool,
    ):
        super().__init__(homologous)
        self._codons = codons
        self._path = path

    @property
    def sequence(self) -> SequenceABC:
        from nmm import Sequence as Seq

        codon_seq = b"".join(c.symbols for c in self._codons)
        seq = Seq(codon_seq, self._path[0].state.alphabet)
        return seq

    def items(self) -> Iterator[Tuple[bytes, CodonStep]]:
        # start = end = 0
        idx = 0
        for step in self._path:
            # end += step.seq_len
            if isinstance(step.state, MuteState):
                yield (b"", step)
            else:
                assert isinstance(step.state, CodonState)
                yield (self._codons[idx].symbols, step)
                idx += 1
            # start = end

    def decode(self, genetic_code: GeneticCode) -> AminoFragment:
        aminos: List[bytes] = []
        steps = []
        # npath = AminoAcidPath()

        amino_acids = genetic_code.amino_acids()
        amino_abc = genetic_code.amino_alphabet
        # aa_alphabet = Alphabet(b"".join(amino_acids))

        i = 0
        for step in self._path:
            if isinstance(step.state, MuteState):
                mstate = MuteState(step.state.name, amino_abc)
                # npath.append_amino_acid_step(mstate, 0)
                steps.append((mstate, 0))
            else:
                assert isinstance(step.state, CodonState)

                aa = genetic_code.amino_acid(self._codons[i])
                aminos.append(aa)

                name = step.state.name
                nstate = NormalState(name, amino_abc, [log(1.0)])
                # npath.append_amino_acid_step(nstate, 1)
                steps.append((nstate, 1))
                i += 1

        return AminoFragment(b"".join(aminos), AminoPath(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
