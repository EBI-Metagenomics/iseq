from typing import TypeVar, Union

from imm import MuteState, NormalState, Path, Results, State, Step
from nmm import AminoAlphabet, CodonState, DNAAlphabet, IUPACAminoAlphabet, RNAAlphabet

from .fragment import Fragment

AminoStep = Step[Union[NormalState, MuteState]]
AminoPath = Path[AminoStep]
AminoFragment = Fragment[AminoAlphabet, NormalState]

CodonStep = Step[Union[CodonState, MuteState]]
CodonPath = Path[CodonStep]

TState = TypeVar("TState", bound=State)

MutableState = Union[TState, MuteState]
MutableStep = Step[Union[TState, MuteState]]
MutablePath = Path[Step[Union[TState, MuteState]]]
MutableResults = Results[Union[TState, MuteState]]

HMMERAlphabet = Union[RNAAlphabet, DNAAlphabet, IUPACAminoAlphabet]


__all__ = [
    "AminoFragment",
    "AminoPath",
    "AminoStep",
    "CodonPath",
    "CodonStep",
    "HMMERAlphabet",
    "MutablePath",
    "MutableResults",
    "MutableState",
    "MutableStep",
]
