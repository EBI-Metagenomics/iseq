from typing import Iterator, Tuple, List
from math import log

from nmm import (
    SequenceABC,
    MuteState,
    FrameState,
    CodonState,
    Codon,
    Interval,
    CBaseAlphabet,
)
from ..fragment import Fragment
from ..codon import CodonFragment, CodonPath
from .path import FramePath
from .step import FrameStep


class FrameFragment(Fragment):
    def __init__(
        self,
        alphabet: CBaseAlphabet,
        sequence: SequenceABC,
        path: FramePath,
        homologous: bool,
    ):
        super().__init__(homologous)
        self._alphabet = alphabet
        self._sequence = sequence
        self._path = path

    @property
    def sequence(self) -> SequenceABC:
        return self._sequence

    def items(self) -> Iterator[Tuple[bytes, FrameStep]]:
        start = end = 0
        for step in self._path:
            end += step.seq_len
            yield (self._sequence.symbols[start:end], step)
            start = end

    def decode(self) -> CodonFragment:
        codons: List[Codon] = []
        steps = []
        # npath = CodonPath()

        start: int = 0
        seq = self.sequence
        for step in self._path:
            if isinstance(step.state, MuteState):
                mstate = MuteState(step.state.name, step.state.alphabet)
                # npath.append_codon_step(mstate, 0)
                steps.append((mstate, 0))
            else:
                assert isinstance(step.state, FrameState)

                subseq = seq.slice(Interval(start, start + step.seq_len))
                # codon = step.state.decode(seq[start : start + step.seq_len])[0]
                codon = Codon(b"XXX", self._alphabet)
                step.state.decode(subseq, codon)
                codons.append(codon)

                name = step.state.name
                abc = step.state.alphabet
                cstate = CodonState(name, abc, {codon: log(1.0)})
                # npath.append_codon_step(cstate, 3)
                steps.append((cstate, 3))

            start += step.seq_len

        return CodonFragment(codons, CodonPath(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
