import os
from io import StringIO

from numpy.testing import assert_equal

from iseq.example import example_filepath
from iseq.gff import read as read_gff
from iseq.testing import assert_same_file


def test_gff_read():
    items = read_gff(example_filepath("duplicate.gff")).items()
    assert_equal(len(items), 14)
    assert_equal(items[3].seqid, "GALNBKIG_00914_ont_01_plus_strand")
    assert_equal(items[6].end, 474)


def test_gff_deduplicate(tmp_path):
    os.chdir(tmp_path)

    gff = read_gff(example_filepath("duplicate.gff"))
    gff.deduplicate()

    gff.write_file("output.gff")

    dedup_file = example_filepath("deduplicate.gff")
    assert_same_file("output.gff", dedup_file)


def test_gff_read_empty():
    file = StringIO("##gff-version 3")
    gff = read_gff(file)
    df = gff.to_dataframe()
    assert_equal(len(df), 0)
