import os
import tempfile
from contextlib import contextmanager
from pathlib import Path

__all__ = ["file_hash", "make_sure_dir_exist", "brotli_decompress", "tmp_cwd"]


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
