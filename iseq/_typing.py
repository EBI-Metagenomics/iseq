from typing import TypeVar, Union

from imm import Alphabet, MuteState, NormalState, Path, Results, State, Step
from nmm import AminoAlphabet, CodonState

from ._fragment import Fragment

AminoStep = Step[Union[NormalState, MuteState]]
AminoPath = Path[AminoStep]
AminoFragment = Fragment[AminoAlphabet, NormalState]

CodonStep = Step[Union[CodonState, MuteState]]
CodonPath = Path[CodonStep]

TState = TypeVar("TState", bound=State)
TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

MutableState = Union[TState, MuteState]
MutableStep = Step[Union[TState, MuteState]]
MutablePath = Path[Step[Union[TState, MuteState]]]
MutableResults = Results[Union[TState, MuteState]]


__all__ = [
    "MutableState",
    "MutableStep",
    "MutablePath",
]
