from pathlib import Path
from typing import IO, Optional, Union

from iseq.gff import GFFItem, GFFWriter
from iseq.profile import ProfileID

__all__ = ["OutputWriter"]


class OutputWriter:
    def __init__(self, file: Union[str, Path, IO[str]]):
        self._gff = GFFWriter(file)
        self._item_idx = 1

    def write_item(
        self,
        seqid: str,
        profid: ProfileID,
        start: int,
        end: int,
        window_length: int,
        att: Optional[dict] = None,
    ):
        if att is None:
            att = dict()

        item_id = f"item{self._item_idx}"
        atts = f"ID={item_id};Profile_name={profid.name}"
        atts += f";Profile_acc={profid.acc}"
        atts += f";Window={window_length}"
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

    def __exit__(self, exception_type, exception_value, traceback):
        del exception_type
        del exception_value
        del traceback
        self.close()
