import os
import shutil
import tempfile
import warnings
from contextlib import contextmanager
from pathlib import Path
from subprocess import check_call
from urllib.parse import urlparse
from urllib.request import urlretrieve


def download(url: str) -> Path:
    """
    Download a file to the current directory.

    Parameters
    ----------
    url
        URL path to the file.

    Returns
    -------
    pathlib.Path
        Path to the local file.
    """
    dest = os.getcwd()
    _filename = os.path.basename(urlparse(url).path)
    filepath = os.path.join(dest, _filename)
    local_filename = urlretrieve(url, filepath)[0]
    return Path(local_filename)


class tempdir:
    def __init__(self):
        self._oldpath = os.getcwd()
        self._dirpath = tempfile.mkdtemp()

    def __enter__(self):
        os.chdir(self._dirpath)

    def __exit__(self, *_):
        os.chdir(self._oldpath)
        try:
            shutil.rmtree(self._dirpath)
        except PermissionError as e:
            warnings.warn(str(e) + "\n. I will ignore it and proceed.")


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


def brotli_decompress(filepath: Path):
    """
    Decompress a brotli file.
    """
    cmd = shutil.which("brotli")
    if cmd is None:
        raise RuntimeError("Could not find the `brotli` command-line tool.")

    if filepath.suffix != ".br":
        raise ValueError("File suffix must be `.br`.")

    output_filepath = filepath.parents[0] / filepath.stem
    if output_filepath.exists():
        raise RuntimeError(f"`{output_filepath}` already exists.")

    check_call([cmd, "-d", "-k", str(filepath), "-o", str(output_filepath)])

    return output_filepath


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


def make_sure_dir_exist(dirpath: Path):
    dirpath.mkdir(parents=True, exist_ok=True)


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
