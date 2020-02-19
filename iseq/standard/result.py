from typing import TypeVar

from nmm.alphabet import CAlphabet

from ..result import SearchResult
from .state import UState

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)

StandardSearchResult = SearchResult[TAlphabet, UState]
