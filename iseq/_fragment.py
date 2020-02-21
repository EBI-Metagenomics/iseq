from typing import TypeVar

from nmm.alphabet import Alphabet
from nmm.fragment import Fragment as FragmentBase
from nmm.path import Path, Step
from nmm.sequence import SequenceABC
from nmm.state import State

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)
TState = TypeVar("TState", bound=State)


__all__ = ["Fragment"]


class Fragment(FragmentBase[TAlphabet, TState]):
    """
    Fragment of a sequence.

    Fragment is a sequence path with homology information.

    Parameters
    ----------
    sequence
        Sequence.
    path
        Path.
    homologous
        Fragment homology.
    """

    def __init__(
        self,
        sequence: SequenceABC[TAlphabet],
        path: Path[Step[TState]],
        homologous: bool,
    ):
        super().__init__(sequence, path)
        self._homologous = homologous

    @property
    def homologous(self) -> bool:
        return self._homologous

    @property
    def sequence(self) -> SequenceABC[TAlphabet]:
        return super().sequence

    @property
    def path(self) -> Path[Step[TState]]:
        return super().path

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
