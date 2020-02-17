from typing import Iterator, Tuple

from nmm import SequenceABC
from ..fragment import Fragment
from .path import FramePath
from .step import FrameStep


class FrameFragment(Fragment):
    def __init__(
        self, sequence: SequenceABC, path: FramePath, homologous: bool,
    ):
        super().__init__(homologous)
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

    # def decode(self) -> CodonFragment:
    #     nseq: List[Codon] = []
    #     npath = CodonPath()

    #     start: int = 0
    #     seq = self.sequence
    #     for step in self._path.steps():
    #         if isinstance(step.state, MuteState):
    #             mstate = MuteState(step.state.name, step.state.alphabet)
    #             npath.append_codon_step(mstate, 0)
    #         else:
    #             assert isinstance(step.state, FrameState)

    #             codon = step.state.decode(seq[start : start + step.seq_len])[0]
    #             nseq.append(codon)

    #             cstate = CodonState(step.state.name, step.state.alphabet, {codon: LOG1})
    #             npath.append_codon_step(cstate, 3)

    #         start += step.seq_len

    #     return CodonFragment(nseq, npath, self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
