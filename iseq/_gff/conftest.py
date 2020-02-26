import pytest


@pytest.fixture
def duplicate_file(tmp_path):
    return _write_file(tmp_path, "duplicate.gff")


@pytest.fixture
def deduplicate_file(tmp_path):
    return _write_file(tmp_path, "deduplicate.gff")


def _write_file(path, filename):
    import importlib_resources as pkg_resources
    from iseq._gff import test

    text = pkg_resources.read_text(test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
