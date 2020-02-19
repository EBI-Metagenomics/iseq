from abc import ABC, abstractmethod
from typing import Iterator, Tuple

from nmm.path import CStep
from nmm.sequence import SequenceABC


class Fragment(ABC):
    """
    Fragment of a sequence.

    Fragment is path with homology information.

    Parameters
    ----------
    homologous : `bool`
        Fragment homology.
    """

    def __init__(self, homologous: bool):
        self._homologous = homologous

    @property
    @abstractmethod
    def sequence(self) -> SequenceABC:
        raise NotImplementedError()

    @abstractmethod
    def items(self) -> Iterator[Tuple[bytes, CStep]]:
        raise NotImplementedError()

    @property
    def homologous(self) -> bool:
        return self._homologous

    def __str__(self) -> str:
        return f"[{self.sequence}]"
