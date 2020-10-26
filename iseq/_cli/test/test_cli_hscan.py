import os
from io import StringIO

import pandas as pd
from assertpy import assert_that, contents_of
from click.testing import CliRunner
from imm.testing import assert_allclose

from iseq import cli
from iseq.example import example_filepath


def test_cli_hscan(tmp_path):
    os.chdir(tmp_path)

    hmm_filepath = example_filepath("2OG-FeII_Oxy_3-nt.hmm")
    fasta_filepath = example_filepath("2OG-FeII_Oxy_3-nt_unilocal.fasta")

    invoke = CliRunner().invoke
    r = invoke(
        cli,
        [
            "hscan",
            str(hmm_filepath),
            str(fasta_filepath),
            "--output",
            "output.gff",
            "--odebug",
            "debug.tsv",
        ],
    )
    assert r.exit_code == 0, r.output
    df = pd.read_csv("debug.tsv", sep="\t", header=0)

    deflines = [
        "2OG-FeII_Oxy_3-seed0",
        "2OG-FeII_Oxy_3-seed1",
        "2OG-FeII_Oxy_3-seed2",
        "2OG-FeII_Oxy_3-seed3",
        "2OG-FeII_Oxy_3-seed4",
        "2OG-FeII_Oxy_3-seed5",
        "2OG-FeII_Oxy_3-seed6",
        "2OG-FeII_Oxy_3-seed7",
        "2OG-FeII_Oxy_3-seed8",
        "2OG-FeII_Oxy_3-seed9",
    ]
    assert all(a == b for a, b in zip(deflines, df["defline"]))
    assert all(a == 0 for a in df["window"])
    assert all(a == 1 for a in df["start"])

    stops = [331, 205, 438, 354, 784, 570, 727, 1186, 757, 322]
    assert all(a == b for a, b in zip(stops, df["stop"]))

    alt_viterbi_score = [
        -7.503633801392808,
        0.4848468637802794,
        -14.883360097987726,
        -14.771241935678885,
        -11.125960638712687,
        -10.816347267999015,
        -13.826119382329795,
        -15.451590982775484,
        -16.470234922310645,
        -4.954388607990875,
    ]
    null_viterbi_score = [
        -0.995475868010347,
        -0.9927026888392518,
        -0.9965796793428908,
        -0.9957693485995353,
        -0.9980880892564644,
        -0.9973709827987554,
        -0.9979383014590296,
        -0.9987358366164624,
        -0.9980199469219536,
        -0.9953496347175764,
    ]

    assert_allclose(df["alt_viterbi_score"], alt_viterbi_score)
    assert_allclose(df["null_viterbi_score"], null_viterbi_score)


_desired_hscan_window200 = """\
defline	window	start	stop	alt_viterbi_score	null_viterbi_score
2OG-FeII_Oxy_3-seed0	0	1	200	-15.527363477178557	-0.6003021143456335
2OG-FeII_Oxy_3-seed0	1	101	300	-7.850984227099249	-0.6003021143456335
2OG-FeII_Oxy_3-seed0	2	201	331	-7.850984227099249	-0.39215716012528823
2OG-FeII_Oxy_3-seed1	0	1	200	-0.5829157613527629	-0.9683717405833878
2OG-FeII_Oxy_3-seed1	1	101	205	-0.935912455911347	-0.5060837237219715
2OG-FeII_Oxy_3-seed2	0	1	200	-15.245591445669424	-0.453820037046305
2OG-FeII_Oxy_3-seed2	1	101	300	-15.200755825320652	-0.453820037046305
2OG-FeII_Oxy_3-seed2	2	201	400	-15.200755825320652	-0.453820037046305
2OG-FeII_Oxy_3-seed2	3	301	438	-16.746364853317925	-0.31242886972534567
2OG-FeII_Oxy_3-seed3	0	1	200	-15.846003451920733	-0.5613543919867068
2OG-FeII_Oxy_3-seed3	1	101	300	-15.104559443585694	-0.5613543919867068
2OG-FeII_Oxy_3-seed3	2	201	354	-15.104559443585694	-0.4315940802711866
2OG-FeII_Oxy_3-seed4	0	1	200	-18.790566488218154	-0.25366478896811806
2OG-FeII_Oxy_3-seed4	1	101	300	-17.251109265820347	-0.25366478896811806
2OG-FeII_Oxy_3-seed4	2	201	400	-17.251109265820347	-0.25366478896811806
2OG-FeII_Oxy_3-seed4	3	301	500	-18.377272876190982	-0.25366478896811806
2OG-FeII_Oxy_3-seed4	4	401	600	-16.978126950374214	-0.25366478896811806
2OG-FeII_Oxy_3-seed4	5	501	700	-12.79724579302745	-0.25366478896811806
2OG-FeII_Oxy_3-seed4	6	601	784	-14.34160029541395	-0.23326963005610857
2OG-FeII_Oxy_3-seed5	0	1	200	-15.58289781361393	-0.3488169166554522
2OG-FeII_Oxy_3-seed5	1	101	300	-15.58289781361393	-0.3488169166554522
2OG-FeII_Oxy_3-seed5	2	201	400	-16.239911646894587	-0.3488169166554522
2OG-FeII_Oxy_3-seed5	3	301	500	-10.923952174733826	-0.3488169166554522
2OG-FeII_Oxy_3-seed5	4	401	570	-10.923952174733826	-0.29623145183302224
2OG-FeII_Oxy_3-seed6	0	1	200	-17.25817778190801	-0.2735395619701748
2OG-FeII_Oxy_3-seed6	1	101	300	-14.983776198390302	-0.2735395619701748
2OG-FeII_Oxy_3-seed6	2	201	400	-13.935243642869697	-0.2735395619701748
2OG-FeII_Oxy_3-seed6	3	301	500	-13.935243642869697	-0.2735395619701748
2OG-FeII_Oxy_3-seed6	4	401	600	-18.109624170704166	-0.2735395619701748
2OG-FeII_Oxy_3-seed6	5	501	700	-17.370952892948864	-0.2735395619701748
2OG-FeII_Oxy_3-seed6	6	601	727	-16.05873236510965	-0.17319590355900516
2OG-FeII_Oxy_3-seed7	0	1	200	-18.38550781244612	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	1	101	300	-17.11206036796655	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	2	201	400	-17.11206036796655	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	3	301	500	-17.153271562051046	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	4	401	600	-18.34984525692567	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	5	501	700	-18.34984525692567	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	6	601	800	-18.450393451326224	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	7	701	900	-18.38428059004831	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	8	801	1000	-15.569063133896439	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	9	901	1100	-15.569063133896439	-0.1677201953474059
2OG-FeII_Oxy_3-seed7	10	1001	1186	-18.657518340285357	-0.1559207846194477
2OG-FeII_Oxy_3-seed8	0	1	200	-17.43605249682422	-0.26270630877972057
2OG-FeII_Oxy_3-seed8	1	101	300	-17.43605249682422	-0.26270630877972057
2OG-FeII_Oxy_3-seed8	2	201	400	-17.61166557445805	-0.26270630877972057
2OG-FeII_Oxy_3-seed8	3	301	500	-16.727798885282404	-0.26270630877972057
2OG-FeII_Oxy_3-seed8	4	401	600	-16.59086399665851	-0.26270630877972057
2OG-FeII_Oxy_3-seed8	5	501	700	-17.415938135704337	-0.26270630877972057
2OG-FeII_Oxy_3-seed8	6	601	757	-17.48162505234466	-0.20594062396802215
2OG-FeII_Oxy_3-seed9	0	1	200	-5.410597716846649	-0.6170547579713324
2OG-FeII_Oxy_3-seed9	1	101	300	-14.25461877252513	-0.6170547579713324
2OG-FeII_Oxy_3-seed9	2	201	322	-16.55367854007932	-0.37519409906799606
"""


def test_cli_hscan_window200(tmp_path):
    os.chdir(tmp_path)

    hmm_filepath = example_filepath("2OG-FeII_Oxy_3-nt.hmm")
    fasta_filepath = example_filepath("2OG-FeII_Oxy_3-nt_unilocal.fasta")

    invoke = CliRunner().invoke
    r = invoke(
        cli,
        [
            "hscan",
            str(hmm_filepath),
            str(fasta_filepath),
            "--output",
            "output.gff",
            "--odebug",
            "debug.tsv",
            "--window",
            200,
            "--hmmer3-compat",
        ],
    )
    assert r.exit_code == 0, r.output
    actual = pd.read_csv("debug.tsv", sep="\t", header=0)

    desired = pd.read_csv(StringIO(_desired_hscan_window200), sep="\t", header=0)

    assert all(a == b for a, b in zip(actual["defline"], desired["defline"]))
    assert all(a == b for a, b in zip(actual["window"], desired["window"]))
    assert all(a == b for a, b in zip(actual["start"], desired["start"]))
    assert all(a == b for a, b in zip(actual["stop"], desired["stop"]))
    assert_allclose(actual["alt_viterbi_score"], desired["alt_viterbi_score"])


_output_window_auto = """\
##gff-version 3
2OG-FeII_Oxy_3-seed0	iseq	.	252	288	0.0	+	.	ID=item1;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed1	iseq	.	91	162	0.0	+	.	ID=item2;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed2	iseq	.	228	272	0.0	+	.	ID=item3;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed3	iseq	.	240	277	0.0	+	.	ID=item4;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed4	iseq	.	591	716	0.0	+	.	ID=item5;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed5	iseq	.	423	441	0.0	+	.	ID=item6;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed6	iseq	.	312	336	0.0	+	.	ID=item7;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed7	iseq	.	383	448	0.0	+	.	ID=item8;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed7	iseq	.	945	989	0.0	+	.	ID=item9;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed8	iseq	.	495	523	0.0	+	.	ID=item10;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
2OG-FeII_Oxy_3-seed9	iseq	.	60	132	0.0	+	.	ID=item11;Target_alph=dna;Profile_name=2OG-FeII_Oxy_3;Profile_alph=dna;Profile_acc=-;Window=630
"""


def test_cli_hscan_window_auto(tmp_path):
    os.chdir(tmp_path)

    hmm_filepath = example_filepath("2OG-FeII_Oxy_3-nt.hmm")
    fasta_filepath = example_filepath("2OG-FeII_Oxy_3-nt_unilocal.fasta")

    invoke = CliRunner().invoke
    r = invoke(
        cli,
        [
            "hscan",
            str(hmm_filepath),
            str(fasta_filepath),
            "--output",
            "output.gff",
            "--window",
            -1,
            "--quiet",
        ],
    )
    assert r.exit_code == 0, r.output

    with open("desired.gff", "w") as file:
        file.write(_output_window_auto)

    assert_that(contents_of("output.gff")).is_equal_to(contents_of("desired.gff"))


def test_cli_hscan_dna_vs_rna(tmp_path):
    os.chdir(tmp_path)

    hmm_filepath = example_filepath("2OG-FeII_Oxy_3-nt.hmm")
    fasta_dna_filepath = example_filepath("2OG-FeII_Oxy_3-nt_unilocal.fasta")
    fasta_rna_filepath = example_filepath("2OG-FeII_Oxy_3-rna_unilocal.fasta")

    invoke = CliRunner().invoke
    r = invoke(
        cli,
        [
            "hscan",
            str(hmm_filepath),
            str(fasta_dna_filepath),
            "--output",
            "output_dna.gff",
            "--window",
            200,
        ],
    )
    assert r.exit_code == 0, r.output

    r = invoke(
        cli,
        [
            "hscan",
            str(hmm_filepath),
            str(fasta_rna_filepath),
            "--output",
            "output_rna.gff",
            "--window",
            200,
        ],
    )
    assert r.exit_code == 0, r.output

    assert_that(contents_of("output_dna.gff")).is_equal_to(
        contents_of("output_rna.gff")
    )
