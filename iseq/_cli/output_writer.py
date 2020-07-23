from pathlib import Path
from typing import IO, Optional, Union

from iseq.gff import GFFItem, GFFWriter
from iseq.profile import ProfileID

__all__ = ["OutputWriter"]


class OutputWriter:
    def __init__(self, file: Union[str, Path, IO[str]], item_prefix="item"):
        self._gff = GFFWriter(file)
        self._item_idx = 1
        self._item_prefix = item_prefix

    def write_item(
        self,
        seqid: str,
        seq_alphabet: str,
        profid: ProfileID,
        prof_alphabet: str,
        start: int,
        end: int,
        window_length: int,
        att: Optional[dict] = None,
    ):
        if att is None:
            att = dict()

        item_id = f"{self._item_prefix}{self._item_idx}"
        atts = f"ID={item_id}"
        atts += f";Target_alph={seq_alphabet}"
        atts += f";Profile_name={profid.name}"
        atts += f";Profile_alph={prof_alphabet}"
        atts += f";Profile_acc={profid.acc}"
        atts += f";Window={window_length}"
        for k in sorted(att.keys()):
            atts += f";{k}={att[k]}"

        item = GFFItem(seqid, "iseq", ".", start + 1, end, "0.0", "+", ".", atts)
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
