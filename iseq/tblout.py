from typing import Dict, NamedTuple, Optional

__all__ = ["TBLScore", "TBLRow", "TBLData", "TBLIndex"]

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

    def find(
        self,
        target_name: str = "-",
        target_acc: str = "-",
        query_name: str = "-",
        query_acc: str = "-",
    ) -> Optional[TBLRow]:
        i = TBLIndex(target_name, target_acc, query_name, query_acc)
        return self.rows.get(i, None)


def decomment(rows):
    for row in rows:
        if row.startswith("#"):
            continue
        yield row
