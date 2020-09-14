from io import StringIO

from numpy import dtype
from numpy.testing import assert_equal

from iseq.example import example_filepath
from iseq.gff import GFF, GFFItem
from iseq.gff import read as read_gff


def test_gff_read():
    gff = read_gff(example_filepath("duplicate.gff"))
    assert_equal(len(gff.items()), 14)
    assert_equal(gff.items()[3].seqid, "GALNBKIG_00914_ont_01_plus_strand")
    assert_equal(gff.items()[6].end, 474)
    df = gff.to_dataframe()
    assert_equal(df.shape, (14, 13))
    assert_equal(
        df.columns.tolist(),
        [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
            "att_ID",
            "att_Profile",
            "att_Epsilon",
            "att_Window",
        ],
    )
    item = GFFItem(
        "AE014075.1:534-908|dna",
        "iseq",
        ".",
        1,
        306,
        "0.0",
        "+",
        ".",
        "ID=chunk_1_item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_alph=dna;Profile_acc=PF01797.17;Window=0;Epsilon=0.01;E-value=4e-28",
    )
    gff.append(item)
    df = gff.to_dataframe()
    assert_equal(df.shape, (15, 18))
    assert_equal(
        df.columns.tolist(),
        [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
            "att_ID",
            "att_Profile",
            "att_Epsilon",
            "att_Window",
            "att_Target_alph",
            "att_Profile_name",
            "att_Profile_alph",
            "att_Profile_acc",
            "att_E-value",
        ],
    )
    assert_equal(df["att_E-value"].dtype, dtype("float64"))


def test_gff_read_empty():
    file = StringIO("##gff-version 3")
    gff = read_gff(file)
    df = gff.to_dataframe()
    assert_equal(len(df), 0)


def test_gff_append():
    gff = GFF("##gff-version 3")
    item0 = GFFItem(
        "AE014075.1:190-252|dna",
        "iseq",
        ".",
        1,
        63,
        "0.0",
        "+",
        ".",
        "ID=chunk_1_item1;Target_alph=dna;Profile_name=Leader_Thr;Profile_alph=dna;Profile_acc=PF08254.12;Window=0;Epsilon=0.01;E-value=7e-13",
    )
    item1 = GFFItem(
        "AE014075.1:534-908|dna",
        "iseq",
        ".",
        1,
        306,
        "0.0",
        "+",
        ".",
        "ID=chunk_1_item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_alph=dna;Profile_acc=PF01797.17;Window=0;Epsilon=0.01;E-value=4e-28",
    )
    gff.append(item0)
    gff.append(item1)
    assert gff.items()[0] == item0
    assert gff.items()[1] == item1
    gff = gff.filter(max_e_value=7e-14)
    assert gff.items()[0] == item1
    assert gff.header == "##gff-version 3"


def test_gff_concat():
    gff = GFF("##gff-version 3")
    item0 = GFFItem(
        "AE014075.1:190-252|dna",
        "iseq",
        ".",
        1,
        63,
        "0.0",
        "+",
        ".",
        "ID=chunk_1_item1;Target_alph=dna;Profile_name=Leader_Thr;Profile_alph=dna;Profile_acc=PF08254.12;Window=0;Epsilon=0.01;E-value=7e-13",
    )
    item1 = GFFItem(
        "AE014075.1:534-908|dna",
        "iseq",
        ".",
        1,
        306,
        "0.0",
        "+",
        ".",
        "ID=chunk_1_item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_alph=dna;Profile_acc=PF01797.17;Window=0;Epsilon=0.01;E-value=4e-28",
    )
    gff.extend([item0, item1])
    assert gff.items()[0] == item0
    assert gff.items()[1] == item1
    gff = gff.filter(max_e_value=7e-14)
    assert gff.items()[0] == item1
    assert gff.header == "##gff-version 3"
