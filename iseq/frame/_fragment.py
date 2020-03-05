from math import log
from typing import List

from nmm.alphabet import BaseAlphabet
from nmm.codon import Codon
from nmm.path import Path
from nmm.sequence import Sequence, SequenceABC
from nmm.state import CodonState, FrameState, MuteState

from .._codon import CodonFragment, CodonPath, CodonStep
from .._fragment import Fragment
from ._typing import FrameStep

__all__ = ["FrameFragment"]


class FrameFragment(Fragment[BaseAlphabet, FrameState]):
    def __init__(
        self,
        sequence: SequenceABC[BaseAlphabet],
        path: Path[FrameStep],
        homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    def decode(self) -> CodonFragment:
        codons: List[Codon] = []
        steps: List[CodonStep] = []

        start = 0
        for step in self.path:
            if isinstance(step.state, MuteState):

                mstate = MuteState(step.state.name, step.state.alphabet)
                steps.append(CodonStep.create(mstate, 0))

            elif isinstance(step.state, FrameState):

                subseq = self.sequence[start : start + step.seq_len]
                codon = step.state.decode(subseq)[1]
                codons.append(codon)

                name = step.state.name
                cstate = CodonState(name, self.sequence.alphabet, {codon: log(1.0)})
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
