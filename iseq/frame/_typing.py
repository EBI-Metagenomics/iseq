from typing import Union

from imm import MuteState, Step
from nmm import BaseAlphabet, FrameState

from .._model import AltModel, Node, NullModel, SpecialNode
from .._result import SearchResults

FrameStep = Step[Union[FrameState, MuteState]]
FrameSearchResults = SearchResults[BaseAlphabet, FrameState]

FrameNode = Node[FrameState]
FrameSpecialNode = SpecialNode[FrameState]
FrameNullModel = NullModel[FrameState]
FrameAltModel = AltModel[FrameState]

__all__ = [
    "FrameStep",
    "FrameSearchResults",
    "FrameNode",
    "FrameSpecialNode",
    "FrameNullModel",
    "FrameAltModel",
]
