from typing import Union

from nmm import CStep, FrameState, MuteState, create_imm_step


class FrameStep(CStep):
    """
    Path step for the frame profile.

    Parameters
    ----------
    state : `Union[MuteState, FrameState]`
        State.
    seq_len : `int`
        Sequence length.
    """

    def __init__(self, state: Union[MuteState, FrameState], seq_len: int):
        super().__init__(create_imm_step(state, seq_len), state)
        self._state = state

    @property
    def state(self) -> Union[MuteState, FrameState]:
        return self._state

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
