import os
import shutil
import tempfile
import warnings
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
    """
    Create and enter a temporary directory.

    The previous working directory is saved and switched back when
    leaving the context. The temporary directory is also recursively
    removed at the context ending.
    """

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
