import pytest


@pytest.fixture
def PF03373(tmp_path):
    return _write_file(tmp_path, "PF03373.hmm")


@pytest.fixture
def GALNBKIG_cut(tmp_path):
    return {
        "fasta": _write_file(tmp_path, "GALNBKIG_cut.fasta"),
        "gff": _write_file(tmp_path, "PF03373_GALNBKIG_cut.gff"),
        "amino.fasta": _write_file(tmp_path, "PF03373_GALNBKIG_cut.amino.fasta"),
        "codon.fasta": _write_file(tmp_path, "PF03373_GALNBKIG_cut.codon.fasta"),
    }


@pytest.fixture
def GALNBKIG_cut_amino_fasta(tmp_path):
    return _write_file(tmp_path, "GALNBKIG_cut.amino.fasta")


@pytest.fixture
def tblout(tmp_path):
    return _write_file(tmp_path, "tblout.txt")


@pytest.fixture
def amino1(tmp_path):
    return _write_file(tmp_path, "amino1.fasta")


@pytest.fixture
def database1(tmp_path):
    return _write_file(tmp_path, "database1.hmm")


@pytest.fixture
def output1(tmp_path):
    return _write_file(tmp_path, "output1.gff")


@pytest.fixture
def output1_evalue(tmp_path):
    return _write_file(tmp_path, "output1_evalue.gff")


def _write_file(path, filename):
    import importlib_resources as pkg_resources
    import iseq._data

    text = pkg_resources.read_text(iseq._data, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
