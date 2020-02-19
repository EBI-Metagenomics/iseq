from typing import Iterator, Tuple

from nmm.sequence import SequenceABC

from ..fragment import Fragment
from .path import StandardPath
from .step import StandardStep


class StandardFragment(Fragment):
    """
    Fragment of the standard profile.

    Parameters
    ----------
    sequence : `SequenceABC`
        Sequence.
    path : `StandardPath`
        Path of the standard profile.
    homologous : `bool`
        Fragment homology.
    """

    def __init__(
        self, sequence: SequenceABC, path: StandardPath, homologous: bool,
    ):
        super().__init__(homologous)
        self._sequence = sequence
        self._path = path

    @property
    def sequence(self) -> SequenceABC:
        return self._sequence

    def items(self) -> Iterator[Tuple[bytes, StandardStep]]:
        start = end = 0
        for step in self._path:
            end += step.seq_len
            yield (self._sequence[start:end], step)
            start = end

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
