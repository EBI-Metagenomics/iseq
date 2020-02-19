from typing import Sequence, Tuple, Union

from nmm import create_imm_path
from nmm.path import CPath
from nmm.state import CodonState, MuteState

from .step import CodonStep


class CodonPath(CPath):
    """
    Path for codon profile.

    Parameters
    ----------
    steps : `Sequence[Tuple[Union[MuteState, CodonState], int]]`
        Steps.
    """

    def __init__(self, steps: Sequence[Tuple[Union[MuteState, CodonState], int]]):
        self.__steps = [CodonStep(step[0], step[1]) for step in steps]
        super().__init__(create_imm_path(self.__steps), self.__steps)

    def __getitem__(self, i) -> CodonStep:
        return self.__steps[i]

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
