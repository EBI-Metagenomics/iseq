from filecmp import cmp

from iseq.file import diff

__all__ = ["assert_same_file"]


def assert_same_file(actual, desired):
    assert cmp(actual, desired, shallow=False), diff(actual, desired)
