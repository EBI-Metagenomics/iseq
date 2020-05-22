from imm import NormalState, Step

from .._fragment import Fragment
from .._model import AltModel, Node, NullModel, SpecialNode
from .._result import SearchResults
from .._typing import MutableStep, TAlphabet

StandardFragment = Fragment[TAlphabet, NormalState]
StandardStep = Step[MutableStep[NormalState]]
StandardSearchResults = SearchResults[TAlphabet, NormalState]

StandardNode = Node[NormalState]
StandardSpecialNode = SpecialNode[NormalState]
StandardNullModel = NullModel[NormalState]
StandardAltModel = AltModel[NormalState]

__all__ = [
    "StandardAltModel",
    "StandardFragment",
    "StandardNode",
    "StandardNullModel",
    "StandardSearchResults",
    "StandardSpecialNode",
    "StandardStep",
]
