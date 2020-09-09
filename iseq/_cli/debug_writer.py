from typing import IO

from iseq.result import DebugRow


class DebugWriter:
    def __init__(self, file: IO[str]):

        self._file = file
        self._file.write(
            "defline\twindow\tstart\tstop\talt_viterbi_score\tnull_viterbi_score\n"
        )

    def write_row(
        self,
        defline: str,
        debug_row: DebugRow,
    ):
        """
        Write item.
        """
        wnum = debug_row.window_number
        start = debug_row.window.start + 1
        stop = debug_row.window.stop
        alt_vit = debug_row.alt_viterbi_score
        null_vit = debug_row.null_viterbi_score
        self._file.write(f"{defline}\t{wnum}\t{start}\t{stop}\t{alt_vit}\t{null_vit}\n")

    def close(self):
        """
        Close the associated stream.
        """
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        del exception_type
        del exception_value
        del traceback
        self.close()
