from typing import TypeVar, Union

from nmm.alphabet import Alphabet
from nmm.path import Path, Step
from nmm.state import MuteState, State
from nmm.result import Results

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
