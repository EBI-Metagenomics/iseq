from typing import NamedTuple, Generator

from ._misc import decomment

TBLScore = NamedTuple("TBLScore", [("e_value", str), ("score", str), ("bias", str)])


TBLRow = NamedTuple(
    "TBLRow",
    [
        ("target_name", str),
        ("target_accession", str),
        ("query_name", str),
        ("query_accession", str),
        ("full_sequence", TBLScore),
        ("best_1_domain", TBLScore),
        ("exp", str),
        ("reg", str),
        ("clu", str),
        ("ov", str),
        ("env", str),
        ("dom", str),
        ("rep", str),
        ("inc", str),
        ("description", str),
    ],
)


def tblout_reader(file) -> Generator[TBLRow, None, None]:
    import csv

    reader = csv.reader(decomment(file), delimiter=" ", skipinitialspace=True)
    for row in reader:
        full_sequence = TBLScore(e_value=row[4], score=row[5], bias=row[6])
        best_1_domain = TBLScore(e_value=row[7], score=row[8], bias=row[9])
        yield TBLRow(
            target_name=row[0],
            target_accession=row[1],
            query_name=row[2],
            query_accession=row[3],
            full_sequence=full_sequence,
            best_1_domain=best_1_domain,
            exp=row[10],
            reg=row[11],
            clu=row[12],
            ov=row[13],
            env=row[14],
            dom=row[15],
            rep=row[16],
            inc=row[17],
            description=row[18],
        )
