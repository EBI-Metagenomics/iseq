from importlib import import_module as _import_module

from . import gff, hmmer3, interval, protein
from ._cli import cli
from ._testit import test

try:
    __version__ = getattr(_import_module("iseq._version"), "version", "x.x.x")
except ModuleNotFoundError:
    __version__ = "x.x.x"

__all__ = [
    "__version__",
    "cli",
    "gff",
    "hmmer3",
    "interval",
    "protein",
    "test",
]
