from typing import Union

from imm import MuteState, Step
from nmm import BaseAlphabet, FrameState

from iseq.model import AltModel, Node, NullModel, SpecialNode
from iseq.result import SearchResults

ProteinStep = Step[Union[FrameState, MuteState]]
ProteinSearchResults = SearchResults[BaseAlphabet, FrameState]

ProteinNode = Node[FrameState]
ProteinSpecialNode = SpecialNode[FrameState]
ProteinNullModel = NullModel[FrameState]
ProteinAltModel = AltModel[FrameState]

__all__ = [
    "ProteinStep",
    "ProteinSearchResults",
    "ProteinNode",
    "ProteinSpecialNode",
    "ProteinNullModel",
    "ProteinAltModel",
]
