from typing import Dict, Generator, NamedTuple

__all__ = ["TBLScore", "TBLRow", "tblout_reader", "TBLData", "TBLIndex"]

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

TBLIndex = NamedTuple(
    "TBLIndex",
    [
        ("target_name", str),
        ("target_accession", str),
        ("query_name", str),
        ("query_accession", str),
    ],
)


class TBLData:
    def __init__(self, file):
        import csv

        self.rows: Dict[TBLIndex, TBLRow] = {}

        reader = csv.reader(decomment(file), delimiter=" ", skipinitialspace=True)
        for line in reader:
            full_sequence = TBLScore(e_value=line[4], score=line[5], bias=line[6])
            best_1_domain = TBLScore(e_value=line[7], score=line[8], bias=line[9])
            index = TBLIndex(line[0], line[1], line[2], line[3])
            row = TBLRow(
                target_name=line[0],
                target_accession=line[1],
                query_name=line[2],
                query_accession=line[3],
                full_sequence=full_sequence,
                best_1_domain=best_1_domain,
                exp=line[10],
                reg=line[11],
                clu=line[12],
                ov=line[13],
                env=line[14],
                dom=line[15],
                rep=line[16],
                inc=line[17],
                description=line[18],
            )

            if index in row:
                raise RuntimeError("Duplicate tblout index.")
            self.rows[index] = row


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


def decomment(rows):
    for row in rows:
        if row.startswith("#"):
            continue
        yield row
