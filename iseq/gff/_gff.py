"""
GFF3 File Format
----------------

The first line of a GFF3 file must be a comment that identifies the version, e.g.

```
##gff-version 3
```

Fields must be tab-separated. Also, all but the final field in each feature line must contain a
value; "empty" columns should be denoted with a '.'.

- seqid: name of the chromosome or scaffold;
- source: name of the program that generated this feature, or the data source (database or project
  name);
- type: type of feature. Must be a term or accession from the SOFA sequence ontology;
- start: Start position of the feature, with sequence numbering starting at 1;
- end: End position of the feature, with sequence numbering starting at 1;
- score: A floating point value;
- strand: defined as + (forward) or - (reverse);
- phase: One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base
  of a codon, '1' that the second base is the first base of a codon, and so on;
- attributes: A semicolon-separated list of tag-value pairs, providing additional information about
  each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF
  documentation for more details;

Example:

```
##gff-version 3
ctg123 . mRNA            1300  9000  .  +  .  ID=mrna0001;Name=sonichedgehog
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mrna0001
ctg123 . exon            1050  1500  .  +  .  ID=exon00002;Parent=mrna0001
ctg123 . exon            3000  3902  .  +  .  ID=exon00003;Parent=mrna0001
ctg123 . exon            5000  5500  .  +  .  ID=exon00004;Parent=mrna0001
ctg123 . exon            7000  9000  .  +  .  ID=exon00005;Parent=mrna0001
```
"""
from __future__ import annotations

import pathlib
from typing import IO, NamedTuple, Union, List
from collections import OrderedDict

__all__ = ["read", "GFF", "GFFItem", "GFFWriter"]

GFFItem = NamedTuple(
    "GFFItem",
    [
        ("seqid", str),
        ("source", str),
        ("type", str),
        ("start", int),
        ("end", int),
        ("score", Union[float, str]),
        ("strand", str),
        ("phase", Union[int, str]),
        ("attributes", str),
    ],
)

_columns = OrderedDict(
    [
        ("seqid", str),
        ("source", str),
        ("type", str),
        ("start", int),
        ("end", int),
        ("score", str),
        ("strand", str),
        ("phase", str),
        ("attributes", str),
    ]
)


def read(file: Union[str, pathlib.Path, IO[str]]) -> GFF:
    from pandas import read_csv

    if isinstance(file, IO):
        close_file = False
    else:
        close_file = True

    if isinstance(file, str):
        file = pathlib.Path(file)

    if isinstance(file, pathlib.Path):
        file = open(file, "r")

    start = file.tell()
    header = file.readline().rstrip()

    names = list(_columns.keys())

    df = read_csv(file, sep="\t", names=names, dtype=_columns)
    gff = GFF(header)
    for _, row in df.iterrows():
        gff.append(GFFItem(**{k: row[k] for k in _columns.keys()}))

    file.seek(start)

    if close_file:
        file.close()

    return gff


class GFF:
    def __init__(self, header: str):
        self._header = header
        self._items: List[GFFItem] = []

    def append(self, item: GFFItem):
        self._items.append(item)

    def _to_dataframe(self):
        from pandas import DataFrame

        df = DataFrame(self._items)
        for name, dtype in _columns.items():
            df[name] = df[name].astype(dtype)
        return df

    def _from_dataframe(self, df):
        self._items = []

        for _, row in df.iterrows():
            self.append(GFFItem(**{k: row[k] for k in _columns.keys()}))

    def deduplicate(self):
        from pandas import concat

        df = self._to_dataframe()
        df = _explode_attributes(df)

        df_stack = []
        for _, df0 in df.groupby(["seqid", "att_Profile"]):
            df_stack.append(_dedup_sequences(df0))
        df = concat(df_stack).sort_index()

        df = _implode_attributes(df)
        self._from_dataframe(df)

    def write_file(self, file: Union[str, pathlib.Path, IO[str]]):
        gff_writer = GFFWriter(file)
        for item in self.items():
            gff_writer.write_item(item)
        gff_writer.close()

    def items(self) -> List[GFFItem]:
        """
        Get the list of all items.
        """
        return self._items

    def __str__(self) -> str:
        return str(self._to_dataframe())


class GFFWriter:
    def __init__(self, file: Union[str, pathlib.Path, IO[str]]):

        if isinstance(file, str):
            file = pathlib.Path(file)

        if isinstance(file, pathlib.Path):
            file = open(file, "w")

        self._file = file
        self._file.write("##gff-version 3\n")

    def write_item(self, item: GFFItem):
        cols = [
            item.seqid,
            item.source,
            item.type,
            str(item.start),
            str(item.end),
            str(item.score),
            item.strand,
            str(item.phase),
            item.attributes,
        ]
        self._file.write("\t".join(cols))
        self._file.write("\n")

    def close(self):
        """
        Close the associated stream.
        """
        self._file.close()


def _explode_attributes(df):
    for i, att in enumerate(df["attributes"]):
        for item in att.split(";"):
            name, value = item.split("=")
            df.at[i, "att_" + name] = value
    return df


def _implode_attributes(df):
    for col in df.columns:
        if col.startswith("att_"):
            del df[col]
    return df


def _dedup_sequences(df):
    df = df.sort_values("start")
    iloc_keep = []
    i = 0
    j = 0
    while i < df.shape[0]:
        start = df.iloc[i]["start"]
        max_end = df.iloc[i]["end"]
        max_end_iloc = i
        j = i + 1
        while j < df.shape[0] and df.iloc[j]["start"] == start:
            end = df.iloc[j]["end"]
            if end > max_end:
                max_end = end
                max_end_iloc = j
            j += 1
        iloc_keep.append(max_end_iloc)
        i = j
        while i < df.shape[0] and df.iloc[i]["end"] <= max_end:
            i += 1
    return df.iloc[iloc_keep]
