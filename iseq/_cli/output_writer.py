from pathlib import Path
from typing import IO, Optional, Union

from iseq.gff import GFFItem, GFFWriter


class OutputWriter:
    def __init__(self, file: Union[str, Path, IO[str]], window: int):
        self._gff = GFFWriter(file)
        self._profile_name = "-"
        self._profile_acc = "-"
        self._window = window
        self._item_idx = 1

    def set_profile(self, name: str = "-", acc: str = "-"):
        self._profile_name = name
        self._profile_acc = acc

    def write_item(self, seqid: str, start: int, end: int, att: Optional[dict] = None):
        if att is None:
            att = dict()

        item_id = f"item{self._item_idx}"
        atts = f"ID={item_id};Profile_name={self._profile_name}"
        atts += f";Profile_acc={self._profile_acc}"
        atts += f";Window={self._window}"
        for k in sorted(att.keys()):
            atts += f";{k}={att[k]}"

        item = GFFItem(seqid, "nmm", ".", start + 1, end, 0.0, "+", ".", atts)
        self._gff.write_item(item)
        self._item_idx += 1
        return item_id

    def close(self):
        """
        Close the associated stream.
        """
        self._gff.close()
