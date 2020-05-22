from typing import TypeVar, Union

from imm import Alphabet, MuteState, Path, Results, State, Step

TState = TypeVar("TState", bound=State)
TAlphabet = TypeVar("TAlphabet", bound=Alphabet)


MutableState = Union[TState, MuteState]
MutableStep = Step[Union[TState, MuteState]]
MutablePath = Path[Step[Union[TState, MuteState]]]
MutableResults = Results[Union[TState, MuteState]]


__all__ = [
    "TState",
    "TAlphabet",
    "MutableState",
    "MutableStep",
    "MutablePath",
]
