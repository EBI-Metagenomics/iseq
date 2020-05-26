from typing import TypeVar

from imm import Alphabet, NormalState, Step

from .._fragment import Fragment
from .._model import AltModel, Node, NullModel, SpecialNode
from .._result import SearchResults
from .._typing import MutableStep

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

StandardFragment = Fragment[TAlphabet, NormalState]
StandardSearchResults = SearchResults[TAlphabet, NormalState]
StandardStep = Step[MutableStep[NormalState]]

StandardAltModel = AltModel[NormalState]
StandardNode = Node[NormalState]
StandardNullModel = NullModel[NormalState]
StandardSpecialNode = SpecialNode[NormalState]

__all__ = [
    "StandardAltModel",
    "StandardFragment",
    "StandardNode",
    "StandardNullModel",
    "StandardSearchResults",
    "StandardSpecialNode",
    "StandardStep",
]
