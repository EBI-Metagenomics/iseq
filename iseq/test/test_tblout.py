from iseq.example import example_filepath
from iseq.tblout import TBLData


def test_tblout():
    with open(example_filepath("tblout.txt"), "r") as file:
        reader = TBLData(file)

        row = reader.find("item2", query_name="Octapeptide", query_acc="PF03373.14")
        assert row.target_name == "item2"
        assert row.full_sequence.e_value == "1.2e-07"
        assert row.best_1_domain.e_value == "1.2e-07"

        row = reader.find("item3", query_name="Octapeptide", query_acc="PF03373.14")
        assert row.target_name == "item3"
        assert row.full_sequence.e_value == "1.2e-07"
        assert row.best_1_domain.e_value == "1.2e-07"
