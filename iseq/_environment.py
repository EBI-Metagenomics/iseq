import os
from pathlib import Path

from appdirs import user_cache_dir

from ._misc import make_sure_dir_exist

ISEQ_CACHE_HOME = Path(
    os.environ.get(
        "ISEQ_CACHE_HOME", default=Path(user_cache_dir("iseq", "EBI-Metagenomics")),
    )
)

__all__ = ["ISEQ_CACHE_HOME"]

make_sure_dir_exist(ISEQ_CACHE_HOME)
