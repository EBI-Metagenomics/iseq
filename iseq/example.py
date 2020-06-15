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
    "PF03373_GALNBKIG_cut.gff": "e7b7322f4d4039f76fd6f476d31b13f3468a2a77b3e83c0b36f66b7d37ad2a18",
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
    "large_rna_seq_output0.gff": "6fbc15524a39ae3748d39bc63c525c9992fda67a02e47c7fc0680daed1bd2fa0",
    "large_rna_seq_output48.gff": "6509d2951f9a3040600b83f326625e3f67922c6a4c3331f5f38338af3872aa49",
    "output1.gff": "053c938f94823054c52729a814aeb015960977be50d2f75b0fcc9df808176e61",
    "output1_evalue.gff": "3837202c76e8a1e5687ad700a81aba19e936f7cbccab94dcc669657cd5d8979e",
    "problematic1.fasta": "8053144931f5d2fbce5cdf0eef4283ada5e2f715bd7bb75075dad31f6114316e",
    "problematic1.hmm": "0fe62f7306510f1f9226e0e26148567a0e1ac4732da7cc89e74d6e208adf955b",
    "tblout.txt": "b1d95803d13017e98e8431c3a926920979470b936706177c5cfbd788585833ff",
    "PF00113.hmm": "3701257b5462049d46a8af0ef4655be3ba1507e41c7e5ceeb1e4251a0639cc1c",
    "PF00113.hmm.br": "991b6badfc261c80e296c77cc73b7c9982e877bf27bfa4c7aa9303e5fd2a95ff",
    "A0ALD9_dna_huge.fasta": "d20c21ff64e4a52bad11b8c0902f5d83195c6a25810d1308e34364438d640f69",
    "A0ALD9_dna_huge.fasta.br": "0b99da9a5b85d8c9118d133a46e07bb45e5ba60b7bc4bec9036028c340b8fa22",
    "PF00113_A0ALD9_dna_huge_output1776.gff": "a9556e86178d17e87f7d94cad53fd3d369f363ca5101b169ba4052dc7104774a",
    "PF00113_A0ALD9_dna_huge_output1776.gff.br": "ea168bf89d52241138b738f0d33710c01686d7ed5be52040381352bbcbb58422",
    "PF00113_A0ALD9_dna_huge_codon1776.fasta": "53804506d7152139bfb4fb9eada69986fcdf7ded7df991dd47279a4d1a42206a",
    "PF00113_A0ALD9_dna_huge_codon1776.fasta.br": "1939b4e25fb09aaab29f1dcce6027a3994e211b651957461c4e04f644b245721",
    "PF00113_A0ALD9_dna_huge_amino1776.fasta": "f817d721ab19c588effdaf6986973cb35063a424bb509017c7d075df81fb5a5f",
    "PF00113_A0ALD9_dna_huge_amino1776.fasta.br": "88e96aec425ca07a8dd6764b7ada9b606dda8a5f52b107f3faf7387c6ca531da",
    "2OG-FeII_Oxy_3-nt.hmm": "c15c9afaf79b2f2502142bb88edc53b9ed7e3e3c29f857f682bbde1c6978ec4b",
    "2OG-FeII_Oxy_3-nt_glocal.fasta": "d7865a392e79bd107db85402839596744c7cbc623a267e332deb742771b7eacc",
    "2OG-FeII_Oxy_3-nt_local.fasta": "d9cdd3cf7e8b6c2acd919ded45b3ee33f2467bafb96d583f2afe74152b8cc0a3",
    "2OG-FeII_Oxy_3-nt_uniglocal.fasta": "1f217a358185a55e4c929ee6ad5de614bcbf16c69d85fdf8418475ba4363acef",
    "2OG-FeII_Oxy_3-nt_unilocal.fasta": "0408df372cb444495dfe1e30d2ba4691a28fdb768be9881ec3fa361516ef056e",
    "2OG-FeII_Oxy_3-sample1_local.fasta": "d0335318c93c7787e9336cc2e104294dfd065c737274ac65b88061d8f54992e2",
    "2OG-FeII_Oxy_3.hmm": "1aba1f7c07ae71bcc8c82fd214f37dc57750f0cd02d96ad7f9e835427dda7689",
    "2OG-FeII_Oxy_3_glocal.fasta": "7bf4069eec53a63225699fafb5aa35750b1f64284a791cdf36775d3f27de7097",
    "2OG-FeII_Oxy_3_local.fasta": "3755eb981aaa1821c3ed68ffd9dd0d9c5651b0d19fdfe817cebd1fea68bbeeef",
    "2OG-FeII_Oxy_3_uniglocal.fasta": "0b8294a682642fed2b0afb83f8fb7b6c3789fec622599f321df0170a4c6e93bc",
    "2OG-FeII_Oxy_3_unilocal.fasta": "a9aaba5b50367938cede6ec54a4d21b97e52c17da1a027cc6e7fb044b21e2820",
    "hmmer3_2OG-FeII_Oxy_3-nt_glocal.fasta.viterbi": "3bb7052efa88f1a83721d4355437c21d560003ce6871bb19209b38bf597fa45e",
    "hmmer3_2OG-FeII_Oxy_3-nt_local.fasta.viterbi": "a0f5285abd9fc607f15c3dde6444c93aa647b71bc15641ee661789992614d747",
    "hmmer3_2OG-FeII_Oxy_3-nt_uniglocal.fasta.viterbi": "2fb774f60bf13ec86fee7186e439db8e1c40174fca3fed12b356c77b6c91ad8f",
    "hmmer3_2OG-FeII_Oxy_3-nt_unilocal.fasta.viterbi": "e41a948662d0f655b24f2c7619fd17921fdd83da8675d1acd522260402610650",
    "hmmer3_2OG-FeII_Oxy_3_glocal.fasta.viterbi": "528d2a9e6eaf635612b529a4adf04e6f4834c3d22385792806fd58d118067a26",
    "hmmer3_2OG-FeII_Oxy_3_local.fasta.viterbi": "01910ec5e4f1e48ee615838f93ee6b8147ca46141909c368bde2788dad851476",
    "hmmer3_2OG-FeII_Oxy_3_uniglocal.fasta.viterbi": "5837154cd694e2405d186251dd00ab53e0768d26c74a988774035093eadc5b8e",
    "hmmer3_2OG-FeII_Oxy_3_unilocal.fasta.viterbi": "40be45e0d60b64c3e122dd4c198281c69082618f682656732ff73d6db48c5368",
    "iseq_2OG-FeII_Oxy_3-nt_glocal.fasta.viterbi": "f9e5eefbcfb1ed6836ec7e69e342be8ad93fea728ae05f215c3919191c438511",
    "iseq_2OG-FeII_Oxy_3-nt_local.fasta.viterbi": "c757d19c66dcf744a56836c8be80e6a345bd8e70b42eb905c981cfc61d379fc8",
    "iseq_2OG-FeII_Oxy_3-nt_uniglocal.fasta.viterbi": "9625f15805cc6c37c4191d45ec79b14f86687de839b802fe02af1ed0581176aa",
    "iseq_2OG-FeII_Oxy_3-nt_unilocal.fasta.viterbi": "37dafb8ad9acbd8bd0e3555a05c89683b5f6c77051f01e6be7b8e5f5fa117204",
    "iseq_2OG-FeII_Oxy_3_glocal.fasta.viterbi": "1179695c1044e075db8bf5c3a74534e721c42704a48ad144b0467a2a911bb8df",
    "iseq_2OG-FeII_Oxy_3_local.fasta.viterbi": "377048eb0617bf54f1412c915ae9ae2f8b69a09979016fe1b4579917792b8fde",
    "iseq_2OG-FeII_Oxy_3_uniglocal.fasta.viterbi": "83e2de881ae0d5a1ba66cf790d7b59314676001b9a992e3a2486fefd4b39b67f",
    "iseq_2OG-FeII_Oxy_3_unilocal.fasta.viterbi": "8595c542009117c33bd881ad9d77299709ff4bf6b7db2d68883674d9344151a1",
    "2OG-FeII_Oxy_3-rna_unilocal.fasta": "8ef0ce57699412bacc6f697a15fba63b970e09f7989758fd6ab50a4c13bc99b3",
    "PF03373_rna_most_likely.fasta": "e10ce706a7ad7b8b437e4f04c0d6eaff8eb716249f80b6ebd40d6f722844ef73",
    "PF03373_rna_oamino.fasta": "e945d4af434477a8b7ccda1bfcce886c0964fc81e729b2f79f90275c7cbfe891",
    "PF03373_rna_output.gff": "ea0635f3dcbc94f924097a40b654c72b697b6567759a3b31fe311e2dec5a4f55",
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
