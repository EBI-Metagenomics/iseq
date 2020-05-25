import os
from pathlib import Path

from appdirs import user_cache_dir

from ._file import make_sure_dir_exist

__all__ = ["ISEQ_CACHE_HOME"]

ISEQ_CACHE_HOME = Path(
    os.environ.get(
        "ISEQ_CACHE_HOME", default=Path(user_cache_dir("iseq", "EBI-Metagenomics")),
    )
)

make_sure_dir_exist(ISEQ_CACHE_HOME)
make_sure_dir_exist(ISEQ_CACHE_HOME / "test_data")
