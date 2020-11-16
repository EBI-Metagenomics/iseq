import os
from math import isfinite
from pathlib import Path

import pytest
from fasta_reader import read_fasta
from hmmer_reader import open_hmmer
from imm import Sequence
from imm.testing import assert_allclose

from iseq.example import example_filepath
from iseq.hmmer3 import create_profile
from iseq.hmmer_model import HMMERModel


@pytest.mark.slow
def test_hmmer3_pfam_viterbi_scores_compat(tmp_path):
    os.chdir(tmp_path)
    db_filepath = example_filepath("Pfam-A.hmm")
    target_filepath = example_filepath("A0ALD9.fasta")
    iseq_scores = loadtxt(example_filepath("Pfam-A_iseq_viterbi_scores.txt"))

    with read_fasta(target_filepath) as fasta:
        target = list(fasta)[0]

    actual_scores = []
    for hmmprof in open_hmmer(db_filepath):
        prof = create_profile(HMMERModel(hmmprof), hmmer3_compat=True)
        seq = Sequence.create(target.sequence.encode(), prof.alphabet)
        search_results = prof.search(seq)
        score = search_results.results[0].alt_viterbi_score
        actual_scores.append(score)

    iseq_scores = loadtxt(example_filepath("Pfam-A_iseq_viterbi_scores.txt"))
    assert_allclose(actual_scores, iseq_scores)

    hmmer3_scores = loadtxt(example_filepath("Pfam-A_hmmer3.3_viterbi_scores.txt"))
    ok = [i for i, score in enumerate(hmmer3_scores) if isfinite(score)]

    actual_scores = [actual_scores[i] for i in ok]
    hmmer3_scores = [hmmer3_scores[i] for i in ok]

    assert_allclose(actual_scores, hmmer3_scores, 3e-2)


def test_hmmer3_viterbi_dna_scores_compat():
    hmmfile = example_filepath("2OG-FeII_Oxy_3-nt.hmm")
    hmmprof = open_hmmer(hmmfile).read_model()
    prof = create_profile(HMMERModel(hmmprof), hmmer3_compat=True)
    for align in ["local", "unilocal", "glocal", "uniglocal"]:

        fastafile = f"2OG-FeII_Oxy_3-nt_{align}.fasta"
        hmmer3_vitfile = f"hmmer3_2OG-FeII_Oxy_3-nt_{align}.fasta.viterbi"
        iseq_vitfile = f"iseq_2OG-FeII_Oxy_3-nt_{align}.fasta.viterbi"

        hmmer3_scores = loadtxt(example_filepath(hmmer3_vitfile))
        iseq_scores = loadtxt(example_filepath(iseq_vitfile))

        # HMMER3 viterbi filter has very low accuracy (2 bytes of integer arithmetic)
        # while we use 8 bytes of floating point arithmetic. Therefore we have
        # to allow for relatively high viterbi score differences.
        ok = [i for i, s in enumerate(hmmer3_scores) if isfinite(s) and abs(s) > 1]

        hmmer3_scores = [hmmer3_scores[i] for i in ok]
        iseq_scores = [iseq_scores[i] for i in ok]

        assert_allclose(iseq_scores, hmmer3_scores, rtol=7e-2)

        actual_scores = []
        with read_fasta(example_filepath(fastafile)) as fasta:
            for target in fasta:
                seq = Sequence.create(target.sequence.encode(), prof.alphabet)
                search_results = prof.search(seq)
                score = search_results.results[0].alt_viterbi_score
                actual_scores.append(score)

        iseq_scores = loadtxt(example_filepath(iseq_vitfile))
        assert_allclose(actual_scores, iseq_scores)


def test_hmmer3_viterbi_amino_scores_compat():
    hmmfile = example_filepath("2OG-FeII_Oxy_3.hmm")
    hmmprof = open_hmmer(hmmfile).read_model()
    prof = create_profile(HMMERModel(hmmprof), hmmer3_compat=True)
    for align in ["local", "unilocal", "glocal", "uniglocal"]:

        fastafile = f"2OG-FeII_Oxy_3_{align}.fasta"
        hmmer3_vitfile = f"hmmer3_2OG-FeII_Oxy_3_{align}.fasta.viterbi"
        iseq_vitfile = f"iseq_2OG-FeII_Oxy_3_{align}.fasta.viterbi"

        hmmer3_scores = loadtxt(example_filepath(hmmer3_vitfile))
        iseq_scores = loadtxt(example_filepath(iseq_vitfile))

        # HMMER3 viterbi filter has very low accuracy (2 bytes of integer arithmetic)
        # while we use 8 bytes of floating point arithmetic. Therefore we have
        # to allow for relatively high viterbi score differences.
        ok = [i for i, s in enumerate(hmmer3_scores) if isfinite(s)]

        hmmer3_scores = [hmmer3_scores[i] for i in ok]
        iseq_scores = [iseq_scores[i] for i in ok]

        assert_allclose(iseq_scores, hmmer3_scores, rtol=1e-2)

        actual_scores = []
        with read_fasta(example_filepath(fastafile)) as fasta:
            for target in fasta:
                seq = Sequence.create(target.sequence.encode(), prof.alphabet)
                search_results = prof.search(seq)
                score = search_results.results[0].alt_viterbi_score
                actual_scores.append(score)

        iseq_scores = loadtxt(example_filepath(iseq_vitfile))
        assert_allclose(actual_scores, iseq_scores)


def loadtxt(filepath: Path):
    arr = []
    with open(filepath, "r") as file:
        for line in file:
            arr.append(float(line))
    return arr
