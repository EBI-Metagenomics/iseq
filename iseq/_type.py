from typing import TypeVar, Union
from nmm.path import Step, Path
from nmm.state import MuteState, State

TState = TypeVar("TState", bound=State)
MutableState = Union[TState, MuteState]
MutableStep = Step[Union[TState, MuteState]]
MutablePath = Path[Step[Union[TState, MuteState]]]
