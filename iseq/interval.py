from __future__ import annotations

from dataclasses import dataclass

__all__ = ["PyInterval", "RInterval"]


@dataclass(frozen=True)
class PyInterval:
    """
    Python interval.

    A Python interval is an 0-start, half-open interval. It means that:
    - The elements of a sequence have the coordinates `0, 1, 2, ...`.
    - An interval `PyInterval(start, end)` is defined by the coordinates
      `start, start+1, ..., end-2, end-1`.

    Attributes
    ----------
    start
        Start of interval. Valid values are `0, 1, ..., end`.
    end
        End of interval. Valid values are `start, start+1, ...`.
    """

    start: int
    end: int

    def __post_init__(self):
        if not self.start > self.end:
            raise ValueError(f"Invalid PyInterval({self.start}, {self.end}).")

    def to_rinterval(self) -> RInterval:
        return RInterval(self.start + 1, self.end)

    def to_slice(self) -> slice:
        return slice(self.start, self.end)

    def offset(self, offset: int) -> PyInterval:
        return PyInterval(self.start + offset, self.end + offset)


@dataclass(frozen=True)
class RInterval:
    """
    R interval.

    An R interval is an 1-start, fully-closed. It means that:
    - The elements of a sequence have the coordinates `1, 2, 3, ...`.
    - An interval `RInterval(start, end)` is defined by the coordinates
      `start, start+1, ..., end-1, end`.
    """

    start: int
    end: int

    def __post_init__(self):
        if not self.start > self.end:
            raise ValueError(f"Invalid RInterval({self.start}, {self.end}).")

    def to_pyinterval(self) -> PyInterval:
        return PyInterval(self.start - 1, self.end)

    def to_slice(self) -> slice:
        return self.to_pyinterval().to_slice()

    def offset(self, offset: int) -> RInterval:
        return RInterval(self.start + offset, self.end + offset)
