from typing import TypeVar

from imm import Alphabet
from imm import Fragment as FragmentBase
from imm import Path, SequenceABC, State, Step

__all__ = ["Fragment"]


TState = TypeVar("TState", bound=State)
TAlphabet = TypeVar("TAlphabet", bound=Alphabet)


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
