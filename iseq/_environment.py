from xdg import XDG_CACHE_HOME

from ._misc import make_sure_dir_exist

ISEQ_CACHE_HOME = XDG_CACHE_HOME / "iseq"

__all__ = ["ISEQ_CACHE_HOME"]

make_sure_dir_exist(ISEQ_CACHE_HOME)
