from typing import IO


class DebugWriter:
    def __init__(self, file: IO[str]):

        self._file = file
        self._file.write(
            "defline\twindow\tstart\tstop\talt_viterbi_score\tnull_viterbi_score\n"
        )

    def write_row(
        self,
        defline: str,
        window: int,
        start: int,
        stop: int,
        alt_viterbi_score: float,
        null_viterbi_score: float,
    ):
        """
        Write item.
        """
        self._file.write(
            f"{defline}\t{window}\t{start}\t{stop}\t{alt_viterbi_score}\t{null_viterbi_score}\n"
        )

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
