from typing import TypeVar, Union

from nmm.alphabet import Alphabet, AminoAlphabet
from nmm.path import Path, Step
from nmm.state import MuteState, NormalState

from .fragment import Fragment

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

AminoStep = Step[Union[NormalState, MuteState]]
AminoPath = Path[AminoStep]
AminoFragment = Fragment[AminoAlphabet, NormalState]

__all__ = ["AminoStep", "AminoPath", "AminoFragment"]
