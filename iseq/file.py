import os
import stat
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Mapping

__all__ = [
    "file_hash",
    "make_sure_dir_exist",
    "brotli_decompress",
    "tmp_cwd",
    "diff",
    "fetch_file",
    "make_executable",
    "assert_file_hash",
]


def file_hash(filepath: Path) -> str:
    import hashlib

    BLOCK_SIZE = 65536

    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        fb = f.read(BLOCK_SIZE)
        while len(fb) > 0:
            sha256.update(fb)
            fb = f.read(BLOCK_SIZE)

    return sha256.hexdigest()


def make_sure_dir_exist(dirpath: Path):
    dirpath.mkdir(parents=True, exist_ok=True)


def brotli_decompress(filepath: Path):
    """
    Decompress a brotli file.
    """
    from brotli import decompress

    if filepath.suffix != ".br":
        raise ValueError("File suffix must be `.br`.")

    output_filepath = filepath.parents[0] / filepath.stem
    if output_filepath.exists():
        raise RuntimeError(f"`{output_filepath}` already exists.")

    with open(filepath, "rb") as input:
        with open(output_filepath, "wb") as output:
            output.write(decompress(input.read()))

    return output_filepath


@contextmanager
def tmp_cwd():
    """
    Create and enter a temporary directory.
    The previous working directory is saved and switched back when
    leaving the context. The temporary directory is also recursively
    removed at the context ending.
    """
    oldpwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:

        os.chdir(tmpdir)
        try:
            yield
        finally:
            os.chdir(oldpwd)


def diff(filepath_a, filepath_b):
    """
    Difference between two files.
    """
    from difflib import ndiff

    with open(filepath_a, "r") as file_a:
        a = file_a.readlines()

    with open(filepath_b, "r") as file_b:
        b = file_b.readlines()

    return "".join([line.expandtabs(1) for line in ndiff(a, b)])


def make_executable(filepath):
    st = os.stat(filepath)
    os.chmod(filepath, st.st_mode | stat.S_IEXEC)


def fetch_file(
    filename: str, cache_subdir: str, url_base: str, filemap: Mapping[str, str]
) -> Path:
    import requests

    from .environment import ISEQ_CACHE_HOME

    filepath = ISEQ_CACHE_HOME / cache_subdir / filename

    if filename + ".br" in filemap:
        zipped = fetch_file(filename + ".br", cache_subdir, url_base, filemap)
        cleanup_invalid_filepath(filepath, filemap[filepath.name])
        if not filepath.exists():
            brotli_decompress(zipped)
    else:
        if filename not in filemap:
            raise ValueError(f"Unknown filename {filename}.")

        cleanup_invalid_filepath(filepath, filemap[filepath.name])

    if not filepath.exists():
        r = requests.get(f"{url_base}/{filename}")
        r.raise_for_status()
        with open(filepath, "wb") as f:
            f.write(r.content)

        assert_file_hash(filepath, filemap[filename])

    return filepath


def assert_file_hash(filepath: Path, filehash: str):
    if file_hash(filepath) != filehash:
        msg = (
            f"Hash mismatch:\n"
            f"  ACTUAL : {file_hash(filepath)}\n"
            f"  DESIRED: {filehash}"
        )
        raise RuntimeError(msg)


def cleanup_invalid_filepath(filepath: Path, filehash: str):
    if filepath.exists() and file_hash(filepath) != filehash:
        filepath.unlink()
