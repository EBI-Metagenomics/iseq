from typing import TypeVar

from nmm.alphabet import CAlphabet
from nmm.state import CState
from nmm.path import CStep, CPath
from nmm.sequence import SequenceABC
from nmm.fragment import Fragment as FragmentBase

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)
TState = TypeVar("TState", bound=CState)


class Fragment(FragmentBase[TAlphabet, TState]):
    """
    Fragment of a sequence.

    Fragment is path with homology information.

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
        path: CPath[CStep[TState]],
        homologous: bool,
    ):
        super().__init__(sequence, path)
        self._homologous = homologous

    @property
    def homologous(self) -> bool:
        return self._homologous

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
