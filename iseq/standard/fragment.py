from typing import TypeVar

from nmm.alphabet import CAlphabet

from ..fragment import Fragment
from .state import UState

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)

StandardFragment = Fragment[TAlphabet, UState]
