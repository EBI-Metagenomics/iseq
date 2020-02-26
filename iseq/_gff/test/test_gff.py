import filecmp

from numpy.testing import assert_equal

from iseq._gff import read_gff


def test_gff_read_gff(duplicate_file):
    gff = read_gff(duplicate_file)
    items = gff.items()
    assert_equal(len(items), 14)
    assert_equal(items[3].seqid, "GALNBKIG_00914_ont_01_plus_strand")
    assert_equal(items[6].end, 474)


def test_gff_deduplicate(tmpdir, duplicate_file, deduplicate_file):
    tmpdir.chdir()

    gff = read_gff(duplicate_file)
    gff.deduplicate()

    gff.write_file("output.gff")

    assert_equal(filecmp.cmp(deduplicate_file, "output.gff", shallow=False), True)
