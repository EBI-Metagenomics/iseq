from iseq.example import example_filepath
from iseq.tblout import tblout_reader


def test_tblout():
    with open(example_filepath("tblout.txt"), "r") as file:
        reader = tblout_reader(file)

        row = next(reader)
        assert row.target_name == "item2"
        assert row.full_sequence.e_value == "1.2e-07"
        assert row.best_1_domain.e_value == "1.2e-07"

        row = next(reader)
        assert row.target_name == "item3"
        assert row.full_sequence.e_value == "1.2e-07"
        assert row.best_1_domain.e_value == "1.2e-07"
