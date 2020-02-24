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
import pathlib
from typing import IO, NamedTuple, Union

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


class GFFWriter:
    def __init__(self, file: Union[str, pathlib.Path, IO[str]]):
        if isinstance(file, str):
            file = pathlib.Path(file)

        if isinstance(file, pathlib.Path):
            file = open(file, "r")

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
