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

import dataclasses
import pathlib
from dataclasses import dataclass
from typing import IO, Any, Dict, Iterable, Iterator, List, Optional, Tuple, Type, Union

import pandas as pd

__all__ = ["read", "GFF", "GFFItem", "GFFWriter"]


@dataclass
class GFFItem:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: str

    def attributes_astuple(self):
        attrs = []
        for item in self.attributes.split(";"):
            name, value = item.partition("=")[::2]
            attrs.append((name, value))
        return tuple(attrs)

    def get_attribute(self, name: str):
        return dict(self.attributes_astuple())[name]

    def set_attribute(self, name: str, value: str):
        attrs = self.attributes_astuple()
        new_attrs = []
        found = False
        for n, v in attrs:
            if n == name:
                new_attrs.append((n, value))
                found = True
            else:
                new_attrs.append((n, v))
        if not found:
            raise ValueError(f"Attribute {name} not found.")

        self._set_attributes(new_attrs)

    def _set_attributes(self, attrs: List[Tuple[str, str]]):
        self.attributes = ";".join(f"{n}={v}" for n, v in attrs)

    def __iter__(self):
        return iter(dataclasses.astuple(self))

    @classmethod
    def field_names(cls) -> List[str]:
        return [f.name for f in dataclasses.fields(cls)]

    @classmethod
    def field_types(cls) -> List[Type[Any]]:
        return [f.type for f in dataclasses.fields(cls)]

    def copy(self) -> GFFItem:
        from copy import copy

        return copy(self)


def read(file: Union[str, pathlib.Path, IO[str]]) -> GFF:
    if isinstance(file, IO):
        close_file = False
    else:
        close_file = True

    if isinstance(file, str):
        file = pathlib.Path(file)

    if isinstance(file, pathlib.Path):
        file = open(file, "r")

    start = file.tell()

    gff = GFF.read_csv(file)

    if close_file:
        file.close()
    else:
        file.seek(start)

    return gff


class GFF:
    def __init__(self, header: str):
        self._header = header
        columns = GFFItem.field_names()
        self._df = pd.DataFrame(columns=columns)
        dtype = _column_types()
        for col, typ in dtype.items():
            self._df[col] = self._df[col].astype(typ)
        self._ravel = False

    @classmethod
    def read_csv(cls: Type[GFF], file: IO[str]) -> GFF:
        header = file.readline().rstrip()

        names = GFFItem.field_names()

        dtype = _column_types()
        df = pd.read_csv(
            file,
            sep="\t",
            names=names,
            dtype=dtype,
            engine="c",
        )

        gff = cls(header)
        gff._df = df

        return gff

    @property
    def dataframe(self) -> pd.DataFrame:
        return self._df

    def extend(self, items: Iterable[GFFItem]):
        reset = self._df.shape[0] != 0
        self._df = pd.concat((self._df, pd.DataFrame.from_dict(items)))
        if reset:
            self._df.reset_index(drop=True, inplace=True)

    def append(self, item: GFFItem):
        if self._ravel:
            raise RuntimeError("Please, unravel it first.")
        self._df = self._df.append(dataclasses.asdict(item), ignore_index=True)

    def ravel(self):
        attrs = []
        for i, att in enumerate(self._df["attributes"]):
            data = {}
            for item in att.split(";"):
                name, value = item.partition("=")[::2]
                data[name] = value
            attrs.append(data)
        atts_df = pd.DataFrame.from_records(attrs)
        atts_df.rename(columns={c: f"att_{c}" for c in atts_df.columns}, inplace=True)
        self._df = pd.concat((self._df, atts_df), axis=1)

        self._ravel = True

    def unravel(self):
        df = self._df
        cols = set(df.columns.tolist()) - set(GFFItem.field_names())
        df.drop(columns=list(cols), inplace=True)
        self._ravel = False

    @property
    def header(self) -> str:
        return self._header

    def write_file(self, file: Union[str, pathlib.Path, IO[str]]):
        gff_writer = GFFWriter(file, self._header)
        for item in self.items:
            gff_writer.write_item(item)
        gff_writer.close()

    @property
    def items(self) -> Iterator[GFFItem]:
        """
        Get the list of all items.
        """
        ncols = len(GFFItem.field_names())
        for row in self._df.itertuples(False):
            yield GFFItem(*row[:ncols])

    def __str__(self) -> str:
        return str(self._df)


class GFFWriter:
    def __init__(
        self, file: Union[str, pathlib.Path, IO[str]], header: Optional[str] = None
    ):

        if isinstance(file, str):
            file = pathlib.Path(file)

        if isinstance(file, pathlib.Path):
            file = open(file, "w")

        self._file = file
        if header is None:
            self._file.write("##gff-version 3\n")
        else:
            self._file.write(f"{header}\n")

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

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        del exception_type
        del exception_value
        del traceback
        self.close()


def _column_types() -> Dict[str, Union[str, type]]:
    names = GFFItem.field_names()
    types = GFFItem.field_types()
    dtype: Dict[str, Union[str, type]] = dict(zip(names, types))
    dtype["strand"] = "category"
    dtype["phase"] = "category"
    return dtype
