from typing import Union

from imm import MuteState, NormalState, Path, Step
from nmm import AminoAlphabet

from ._fragment import Fragment

AminoStep = Step[Union[NormalState, MuteState]]
AminoPath = Path[AminoStep]
AminoFragment = Fragment[AminoAlphabet, NormalState]

__all__ = ["AminoStep", "AminoPath", "AminoFragment"]
