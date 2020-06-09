from typing import TypeVar

from imm import Alphabet, NormalState, Step

from iseq.fragment import Fragment
from iseq.model import AltModel, Node, NullModel, SpecialNode
from iseq.result import SearchResults
from iseq.typing import MutableStep

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

HMMER3AltModel = AltModel[NormalState]
HMMER3Fragment = Fragment[TAlphabet, NormalState]
HMMER3Node = Node[NormalState]
HMMER3NullModel = NullModel[NormalState]
HMMER3SearchResults = SearchResults[TAlphabet, NormalState]
HMMER3SpecialNode = SpecialNode[NormalState]
HMMER3Step = Step[MutableStep[NormalState]]

__all__ = [
    "HMMER3AltModel",
    "HMMER3Fragment",
    "HMMER3Node",
    "HMMER3NullModel",
    "HMMER3SearchResults",
    "HMMER3SpecialNode",
    "HMMER3Step",
]
