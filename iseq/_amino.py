from typing import Union

from nmm.alphabet import AminoAlphabet
from nmm.path import Path, Step
from nmm.state import MuteState, NormalState

from ._fragment import Fragment


AminoStep = Step[Union[NormalState, MuteState]]
AminoPath = Path[AminoStep]
AminoFragment = Fragment[AminoAlphabet, NormalState]

__all__ = ["AminoStep", "AminoPath", "AminoFragment"]
