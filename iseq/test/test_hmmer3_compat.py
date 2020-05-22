import pytest
from numpy import asarray, isfinite, loadtxt
from numpy.testing import assert_allclose
from tqdm import tqdm

from fasta_reader import open_fasta
from hmmer_reader import open_hmmer
from imm import Sequence
from iseq._misc import tmp_cwd
from iseq._test_data import get_filepath
from iseq.standard import create_hmmer3_profile

from .._hmmdata import HMMData


@pytest.mark.slow
def test_hmmer3_viterbi_scores_compat():
    with tmp_cwd():
        profiles_filepath = get_filepath("Pfam-A.hmm")
        target_filepath = get_filepath("A0ALD9.fasta")
        scores_iseq = loadtxt(get_filepath("Pfam-A_iseq_viterbi_scores.txt"))

        with open_fasta(target_filepath) as fasta:
            target = list(fasta)[0]

        breakpoint()
        actual_scores = []
        # i = 0
        for hmmprof in tqdm(open_hmmer(profiles_filepath)):
            prof = create_hmmer3_profile(HMMData(hmmprof), hmmer3_compat=True)
            seq = Sequence.create(target.sequence.encode(), prof.alphabet)
            search_results = prof.search(seq, 0)
            actual_score = search_results.results[0].viterbi_score
            actual_scores.append(actual_score)
            # if i == 100:
            #     break
            # i += 1

        import iseq._model
        import pandas as pd

        elapsed = iseq._model._ELAPSED
        df = pd.DataFrame(elapsed, columns=("M", "elapsed"))
        df.to_csv("/Users/horta/viterbi_perf.csv")

        actual_scores = asarray(actual_scores)
        scores_hmmer3 = loadtxt(get_filepath("Pfam-A_hmmer3.3_viterbi_scores.txt"))
        ok = isfinite(scores_hmmer3)
        assert_allclose(actual_scores[ok], scores_hmmer3[ok], 3e-2)

        scores_iseq = loadtxt(get_filepath("Pfam-A_iseq_viterbi_scores.txt"))
        assert_allclose(actual_scores[ok], scores_iseq[ok])
