from typing import Union

from nmm import NormalState, CStep, MuteState, create_imm_step


class AminoStep(CStep):
    """
    Path step for the codon profile.

    Parameters
    ----------
    state : `Union[MuteState, NormalState]`
        State.
    seq_len : `int`
        Sequence length.
    """

    def __init__(self, state: Union[MuteState, NormalState], seq_len: int):
        super().__init__(create_imm_step(state, seq_len), state)
        self._state = state

    @property
    def state(self) -> Union[MuteState, NormalState]:
        return self._state

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
