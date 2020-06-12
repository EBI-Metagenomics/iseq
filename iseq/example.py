from pathlib import Path

from .environment import ISEQ_CACHE_HOME
from .file import brotli_decompress, file_hash

__all__ = ["example_filepath"]

_filenames = {
    "A0ALD9.fasta": "03d238022670c2cdd0e99bb1f6f9c8e7134f995e5da431ef0a3cadf2836bcfdc",
    "GALNBKIG_cut.fasta": "a0105af1bc1d57fc4c577b67fc1207a5aa2cb23cf163c39d026285091f993017",
    "PF03373.hmm": "0e725a276ae927ccef1432b5c78a4307a8c8770455e3ac1808f65a73681fbb77",
    "PF03373_GALNBKIG_cut.amino.fasta": "095d7ab79e35b54273f50a83d0c1d3da092c5d223c4e7e2fc230c0b780e37d2e",
    "PF03373_GALNBKIG_cut.codon.fasta": "8f38982e6a349d558e01640d69649b6504176ea93a174fe94fbaa5fc025762c0",
    "PF03373_GALNBKIG_cut.gff": "4ada39ff8b7d26d06795813023217a3f6798cb3a763a4968330ecaa6eac7154c",
    "PF07476.11.hmm": "fe8ea5da2ffecd4ca6614165fefd29de18b147f8de0a9339361fef56579acdef",
    "PF07476.11.hmm.br": "08d64e5ce06e3b7710b90200eca59f92bdb9507dce0cd3cdedccc65df9d57327",
    "PF15449.6.hmm": "7dcc7b807306d2575700d60625ffe47d97224356610fc0ef9621b23b25fad18a",
    "PF15449.6.hmm.br": "d212c46b45aa9ede405aa40cc818eb3ad197082d723305b03e34fe69f4f8c5d7",
    "Pfam-A.hmm": "deaf15657dc19b0e6c3d7370c4c94a505b109dde103c327862213679958b5dd3",
    "Pfam-A.hmm.br": "64e2626f1704d18a78a00c58eab6c3a417b97787db101d0ce2cc62a9a860ecf8",
    "Pfam-A_hmmer3.3_viterbi_scores.txt": "a9087a7bdd77dcfe8a467a949f3132a8b33a6768dafef419d51394552b67e33c",
    "Pfam-A_hmmer3.3_viterbi_scores.txt.br": "143ad8ef5430e27791e73b2b7044a419fe34f6472201dc337b8e99c507d72c8c",
    "Pfam-A_iseq_viterbi_scores.txt": "af4ecc5a6f04614b00a1697f3849b1f4e677fee172cdf565e0d6bc3bf107609a",
    "Pfam-A_iseq_viterbi_scores.txt.br": "7224a2f7778f0e9dd265a2b564b1255f28dbea9f43a3b330503a041f2e3a6c40",
    "amino1.fasta": "5e0ec53f3a49bae4980c62a124a2df315845767d79b3596a594df833ce723c02",
    "database1.hmm": "c079de36d5f6c3d2e59a8c708a4524f2c0bc7409f5d43ab3ac8fda2045bdabdf",
    "deduplicate.gff": "f5c8b1fcf83b27f563c1ae7d5cd2df0d2639ab4624a50cc51e876283df96090b",
    "duplicate.gff": "522f425659a9299914d2143d30c9ecd4778f7902c75092d1c7ed67b07ebb6938",
    "ecori.hmm": "615843744493c57e93f616bc13fc46b81935becb738766008f907065761d3d19",
    "large_rna_seq.fasta": "d5470b074541a4c7e54f294e9a459bb279e097f1476feeda07ad1c1ce5eac12e",
    "large_rna_seq_amino0.fasta": "aeb020d25ad221062a9475ea1642add361667cd9a57c4455f43c9a77a0d328e1",
    "large_rna_seq_amino48.fasta": "aeb020d25ad221062a9475ea1642add361667cd9a57c4455f43c9a77a0d328e1",
    "large_rna_seq_codon0.fasta": "1e7c5a4c7faeeb3a12b373ddd97a94f5aff7b33360a450dd991ec5d795ef7fae",
    "large_rna_seq_codon48.fasta": "1e7c5a4c7faeeb3a12b373ddd97a94f5aff7b33360a450dd991ec5d795ef7fae",
    "large_rna_seq_output0.gff": "8e4a3f1a7fff2241e07e81c069ef5c27086b875ffc39cb3e92f6c9e6c44f62ab",
    "large_rna_seq_output48.gff": "73f3196445d2db76125ffccadcb5efe91bd719c24e979029c6dd812be48641a7",
    "output1.gff": "053c938f94823054c52729a814aeb015960977be50d2f75b0fcc9df808176e61",
    "output1_evalue.gff": "3837202c76e8a1e5687ad700a81aba19e936f7cbccab94dcc669657cd5d8979e",
    "problematic1.fasta": "8053144931f5d2fbce5cdf0eef4283ada5e2f715bd7bb75075dad31f6114316e",
    "problematic1.hmm": "0fe62f7306510f1f9226e0e26148567a0e1ac4732da7cc89e74d6e208adf955b",
    "tblout.txt": "b1d95803d13017e98e8431c3a926920979470b936706177c5cfbd788585833ff",
    "PF00113.hmm": "3701257b5462049d46a8af0ef4655be3ba1507e41c7e5ceeb1e4251a0639cc1c",
    "PF00113.hmm.br": "991b6badfc261c80e296c77cc73b7c9982e877bf27bfa4c7aa9303e5fd2a95ff",
    "A0ALD9_dna_huge.fasta": "d20c21ff64e4a52bad11b8c0902f5d83195c6a25810d1308e34364438d640f69",
    "A0ALD9_dna_huge.fasta.br": "0b99da9a5b85d8c9118d133a46e07bb45e5ba60b7bc4bec9036028c340b8fa22",
    "PF00113_A0ALD9_dna_huge_output1776.gff": "9a2ac951e91ef6689f24e331e75ce2ba622d3aa89f4f0e618da776909ffef836",
    "PF00113_A0ALD9_dna_huge_output1776.gff.br": "53bd9a2a0dfcb8c41907e92d5c3226df3d6e2f064b8db07cb143e347f0d3f4a6",
    "PF00113_A0ALD9_dna_huge_codon1776.fasta": "53804506d7152139bfb4fb9eada69986fcdf7ded7df991dd47279a4d1a42206a",
    "PF00113_A0ALD9_dna_huge_codon1776.fasta.br": "1939b4e25fb09aaab29f1dcce6027a3994e211b651957461c4e04f644b245721",
    "PF00113_A0ALD9_dna_huge_amino1776.fasta": "f817d721ab19c588effdaf6986973cb35063a424bb509017c7d075df81fb5a5f",
    "PF00113_A0ALD9_dna_huge_amino1776.fasta.br": "88e96aec425ca07a8dd6764b7ada9b606dda8a5f52b107f3faf7387c6ca531da",
}


def example_filepath(filename: str) -> Path:
    import requests

    url = "https://iseq-py.s3.amazonaws.com"
    test_data_folder = ISEQ_CACHE_HOME / "test_data"
    filepath = test_data_folder / filename

    if filename + ".br" in _filenames:
        zipped = example_filepath(filename + ".br")
        cleanup_invalid_filepath(filepath)
        if not filepath.exists():
            brotli_decompress(zipped)
    else:
        if filename not in _filenames:
            raise ValueError(f"Unknown filename {filename}.")

        cleanup_invalid_filepath(filepath)

    if not filepath.exists():
        r = requests.get(f"{url}/{filename}")
        r.raise_for_status()
        with open(filepath, "wb") as f:
            f.write(r.content)

        if file_hash(filepath) != _filenames[filename]:
            msg = (
                f"Hash mismatch:\n"
                f"  ACTUAL : {file_hash(filepath)}\n"
                f"  DESIRED: {_filenames[filename]}"
            )
            raise RuntimeError(msg)

    return filepath


def cleanup_invalid_filepath(filepath: Path):
    if filepath.exists() and file_hash(filepath) != _filenames[filepath.name]:
        filepath.unlink()
