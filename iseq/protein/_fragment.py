from math import log
from typing import List

from imm import MuteState, Path, Sequence, SequenceABC
from nmm import BaseAlphabet, Codon, CodonProb, CodonState, FrameState

from iseq.codon import CodonFragment
from iseq.fragment import Fragment
from iseq.typing import CodonPath, CodonStep

from .typing import ProteinStep

__all__ = ["ProteinFragment"]


class ProteinFragment(Fragment[BaseAlphabet, FrameState]):
    def __init__(
        self,
        sequence: SequenceABC[BaseAlphabet],
        path: Path[ProteinStep],
        homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    def decode(self) -> CodonFragment:
        codons: List[Codon] = []
        steps: List[CodonStep] = []

        start = 0
        for step in self.path:
            if isinstance(step.state, MuteState):

                mstate = MuteState.create(step.state.name, step.state.alphabet)
                steps.append(CodonStep.create(mstate, 0))

            elif isinstance(step.state, FrameState):

                subseq = self.sequence[start : start + step.seq_len]
                codon = step.state.decode(subseq)[1]
                codons.append(codon)

                name = step.state.name
                codonp = CodonProb.create(self.sequence.alphabet)
                codonp.set_lprob(codon, log(1.0))
                cstate = CodonState.create(name, codonp)
                steps.append(CodonStep.create(cstate, 3))

            else:
                raise RuntimeError("Wrong state type.")

            start += step.seq_len

        FrameSequence = Sequence[BaseAlphabet]
        abc = self.sequence.alphabet
        sequence = FrameSequence.create(b"".join(c.symbols for c in codons), abc)
        return CodonFragment(sequence, CodonPath.create(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
