from typing import Sequence, Tuple, Union

from nmm import CPath, FrameState, MuteState, create_imm_path

from .step import FrameStep


class FramePath(CPath):
    """
    Path for frame profile.

    Parameters
    ----------
    steps : `Sequence[Tuple[Union[MuteState, FrameState], int]]`
        Steps.
    """

    def __init__(self, steps: Sequence[Tuple[Union[MuteState, FrameState], int]]):
        self.__steps = [FrameStep(step[0], step[1]) for step in steps]
        super().__init__(create_imm_path(self.__steps), self.__steps)

    def __getitem__(self, i) -> FrameStep:
        return self.__steps[i]

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
