import os
from math import isfinite
from pathlib import Path

import pytest
from fasta_reader import open_fasta
from hmmer_reader import open_hmmer

from imm import Sequence
from imm.testing import assert_allclose
from iseq import file_example
from iseq.standard import create_profile

from .._hmmdata import HMMData


@pytest.mark.slow
def test_hmmer3_viterbi_scores_compat(tmp_path):
    os.chdir(tmp_path)
    db_filepath = file_example("Pfam-A.hmm")
    target_filepath = file_example("A0ALD9.fasta")
    iseq_scores = loadtxt(file_example("Pfam-A_iseq_viterbi_scores.txt"))

    with open_fasta(target_filepath) as fasta:
        target = list(fasta)[0]

    actual_scores = []
    for hmmprof in open_hmmer(db_filepath):
        prof = create_profile(HMMData(hmmprof), hmmer3_compat=True)
        seq = Sequence.create(target.sequence.encode(), prof.alphabet)
        search_results = prof.search(seq, 0)
        score = search_results.results[0].viterbi_score
        actual_scores.append(score)

    iseq_scores = loadtxt(file_example("Pfam-A_iseq_viterbi_scores.txt"))
    assert_allclose(actual_scores, iseq_scores)

    hmmer3_scores = loadtxt(file_example("Pfam-A_hmmer3.3_viterbi_scores.txt"))
    ok = [i for i, score in enumerate(hmmer3_scores) if isfinite(score)]

    actual_scores = [actual_scores[i] for i in ok]
    hmmer3_scores = [hmmer3_scores[i] for i in ok]

    assert_allclose(actual_scores, hmmer3_scores, 3e-2)


def loadtxt(filepath: Path):
    arr = []
    with open(filepath, "r") as file:
        for line in file:
            arr.append(float(line))
    return arr
