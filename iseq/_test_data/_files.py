from pathlib import Path
from .._misc import brotli_decompress, download, tmp_cwd, file_hash
from .._environment import ISEQ_CACHE_HOME

_filenames = {
    "A0ALD9_dna_huge.fasta": "d20c21ff64e4a52bad11b8c0902f5d83195c6a25810d1308e34364438d640f69",
    "A0ALD9.fasta": "03d238022670c2cdd0e99bb1f6f9c8e7134f995e5da431ef0a3cadf2836bcfdc",
    "PF00113.hmm": "3701257b5462049d46a8af0ef4655be3ba1507e41c7e5ceeb1e4251a0639cc1c",
    "PF07476.11.hmm": "fe8ea5da2ffecd4ca6614165fefd29de18b147f8de0a9339361fef56579acdef",
    "PF15449.6.hmm": "7dcc7b807306d2575700d60625ffe47d97224356610fc0ef9621b23b25fad18a",
    "Pfam-A.hmm": "deaf15657dc19b0e6c3d7370c4c94a505b109dde103c327862213679958b5dd3",
    "Pfam-A_hmmer3.3_viterbi_scores.txt": "a9087a7bdd77dcfe8a467a949f3132a8b33a6768dafef419d51394552b67e33c",
    "PF00113_A0ALD9_dna_huge_output1776.gff": "4276b590cbf94df33d3c75ff091a910019a6dae6fdcb32f084a490bfb2a30901",
    "PF00113_A0ALD9_dna_huge_codon1776.fasta": "53804506d7152139bfb4fb9eada69986fcdf7ded7df991dd47279a4d1a42206a",
    "PF00113_A0ALD9_dna_huge_amino1776.fasta": "f817d721ab19c588effdaf6986973cb35063a424bb509017c7d075df81fb5a5f",
    "Pfam-A_iseq_viterbi_scores.txt": "",
}


def get_filepath(filename: str):

    if filename not in _filenames:
        raise ValueError(f"Unknown filename {filename}.")

    test_data_folder = Path(ISEQ_CACHE_HOME / "test_data")
    filepath = test_data_folder / filename

    if filepath.exists() and file_hash(filepath) != _filenames[filename]:
        filepath.unlink()

    if not filepath.exists():
        with tmp_cwd():
            url_base = "https://rest.s3for.me/iseq"
            compressed_filepath = download(f"{url_base}/{filename}.br")
            brotli_decompress(compressed_filepath).rename(filepath)

    if file_hash(filepath) != _filenames[filename]:
        msg = (
            f"Hash mismatch:\n"
            f"  ACTUAL : {file_hash(filepath)}\n"
            f"  DESIRED: {_filenames[filename]}"
        )
        raise RuntimeError(msg)

    return filepath
