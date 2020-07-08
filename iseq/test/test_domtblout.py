import pytest

from iseq.domtblout import DomTBLData
from iseq.example import example_filepath


def test_domtblout():
    with open(example_filepath("domtblout.txt"), "r") as file:
        tbl = DomTBLData(file)
        it = iter(tbl)
        row = next(it)

        assert row.target.name == "Leader_Thr"
        assert row.target.accession == "PF08254.12"
        assert row.target.length == 22

        assert row.query.name == "AE014075.1:190-252|amino|11"
        assert row.query.accession == "-"
        assert row.query.length == 21

        assert row.full_sequence.e_value == "5.3e-10"
        assert row.full_sequence.score == "38.8"
        assert row.full_sequence.bias == "17.5"

        assert row.domain.id == 1
        assert row.domain.size == 1
        assert row.domain.c_value == "3e-14"
        assert row.domain.i_value == "5.5e-10"
        assert row.domain.score == "38.8"
        assert row.domain.bias == "17.5"

        assert row.hmm_coord.start == 1
        assert row.hmm_coord.stop == 22

        assert row.ali_coord.start == 1
        assert row.ali_coord.stop == 21

        assert row.env_coord.start == 1
        assert row.env_coord.stop == 21

        assert row.acc == "0.99"
        assert row.description == "Threonine leader peptide"

        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        next(it)
        row = next(it)

        assert row.acc == "0.89"
        assert row.description == "Probable molybdopterin binding domain"

        with pytest.raises(StopIteration):
            next(it)
