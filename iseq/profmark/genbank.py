from pathlib import Path

import nmm
from Bio import Entrez, GenBank, SeqIO
from Bio.Alphabet import IUPAC, DNAAlphabet, RNAAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from iseq.codon_table import CodonTable
from iseq.file import assert_file_hash, cleanup_invalid_filepath
from iseq.gencode import GeneticCode

__all__ = ["get_accession", "download", "extract_cds"]


def get_accession(filepath: Path) -> str:
    with open(filepath, "r") as file:
        rec = next(GenBank.parse(file))
        return rec.version


def download(accession: str, rettype: str, output: Path):
    """
    Parameters
    ----------
    accession
        Accession number.
    rettype
        Accepted values are ``"gb"`` and ``"fasta"``.
    """

    cleanup_invalid_filepath(output, accession_hash(accession)[rettype])
    if not output.exists():
        download_efetch(output, accession, rettype)
    assert_file_hash(output, accession_hash(accession)[rettype])


def extract_cds(gb_filepath: Path, amino_filepath: Path, nucl_filepath: Path):

    accession = get_accession(gb_filepath)

    cleanup_invalid_filepath(nucl_filepath, accession_hash(accession)["nucl"])
    cleanup_invalid_filepath(amino_filepath, accession_hash(accession)["amino"])

    if nucl_filepath.exists() and amino_filepath.exists():
        return

    nucl_output = open(nucl_filepath, "w")
    amino_output = open(amino_filepath, "w")

    rec: SeqRecord = next(SeqIO.parse(gb_filepath, "genbank"))

    nucl_name = rec.annotations["molecule_type"].lower()
    assert nucl_name in ["dna", "rna"]

    starts = set()
    # TODO: remove the subset
    jj = 0
    for feature in tqdm(rec.features, desc="Features"):
        if feature.type != "CDS":
            continue

        if "protein_id" not in feature.qualifiers:
            continue

        if feature.strand != 1:
            continue

        nucl_rec: SeqRecord = feature.extract(rec)
        if is_alphabet_ambiguous(nucl_rec.seq):
            continue

        assert len(feature.qualifiers["translation"]) == 1
        if is_extended_protein(feature.qualifiers["translation"][0]):
            continue

        assert accession == nucl_rec.id

        # TODO: remove the subset
        if jj == 10:
            break
        jj += 1

        amino_rec: SeqRecord = SeqRecord(
            Seq(feature.qualifiers["translation"][0], IUPAC.protein),
            id=nucl_rec.id,
            name=nucl_rec.name,
            description=nucl_rec.description,
        )

        assert len(feature.qualifiers["transl_table"]) == 1
        transl_table = int(feature.qualifiers["transl_table"][0])
        nucl_seq, amino_seq = remove_stop_codon(
            nucl_rec.seq, amino_rec.seq, transl_table
        )
        nucl_rec.seq = nucl_seq
        amino_rec.seq = amino_seq

        assert len(feature.qualifiers["codon_start"]) == 1
        assert feature.qualifiers["codon_start"][0] == "1"

        start = int(feature.location.start) + 1
        end = start + len(nucl_rec.seq) - 1

        if start in starts:
            raise ValueError("Feature with same start position.")
        starts.add(start)

        nucl_rec.id = f"{accession}:{start}-{end}|{nucl_name}"
        amino_rec.id = f"{accession}:{start}-{end}|amino|{transl_table}"

        nucl_output.write(nucl_rec.format("fasta"))
        amino_output.write(amino_rec.format("fasta"))

    nucl_output.close()
    amino_output.close()

    # TODO: uncomment it
    # assert_file_hash(nucl_filepath, accession_hash(accession)["nucl"])
    # assert_file_hash(amino_filepath, accession_hash(accession)["amino"])


def download_efetch(filepath: Path, accession: str, rettype: str):
    Entrez.email = "horta@ebi.ac.uk"
    efetch = Entrez.efetch

    with efetch(db="nuccore", id=accession, rettype=rettype, retmode="text") as handle:
        with open(filepath, "w") as file:
            file.write(handle.read())


def accession_hash(accession: str):
    accessions = {
        "AE014075.1": {
            "amino": "d14d16c7455caeb1411101bf54967643adfe7619a12b8816e3c755613f337f38",
            "nucl": "96308d9d2da397993f596f28559935ac0c1f98ab817fc92627e2217826ba1cef",
            "gb": "c26ac7e4208de948a0d3644423b48d9eaab3f4a3bf0e3ce1b64af6ded8e2d7af",
            "fasta": "d7e170b7758aa367b45f99f1c87e31977c1f1309b687a9e36d4402b3ccf65d82",
        }
    }

    return accessions[accession]


def is_alphabet_ambiguous(seq):

    if isinstance(seq.alphabet, DNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_dna.letters))
        if remains > 0:
            return True

    elif isinstance(seq.alphabet, RNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_rna.letters))
        if remains > 0:
            return True

    else:
        raise ValueError("Unkown alphabet.")

    return False


def get_nucl_alphabet(seq):
    if isinstance(seq.alphabet, DNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_dna.letters))
        assert remains == 0
        return "dna"

    if isinstance(seq.alphabet, RNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_rna.letters))
        assert remains == 0
        return "rna"

    raise ValueError("Unkown alphabet.")


def is_extended_protein(seq: str):
    remains = len(set(seq) - set(IUPAC.protein.letters))
    return remains > 0


def remove_stop_codon(nucl_seq: Seq, amino_seq: Seq, trans_table_num: int):
    abc_name = get_nucl_alphabet(nucl_seq)
    if abc_name == "dna":
        base_abc = nmm.DNAAlphabet()
    else:
        base_abc = nmm.RNAAlphabet()

    amino_abc = nmm.IUPACAminoAlphabet()
    codon_table = CodonTable(base_abc, amino_abc, GeneticCode(id=trans_table_num))
    nucl_str = str(nucl_seq)

    aminos = []
    for start in range(0, len(nucl_str), 3):
        stop = min(start + 3, len(nucl_str))
        codon = nmm.Codon.create(nucl_str[start:stop].encode(), base_abc)
        if start == 0:
            if codon in codon_table.start_codons:
                aminos.append("M")
                continue
        aminos.append(codon_table.decode(codon).decode())

    amino_str = "".join(aminos)
    assert amino_str[-1] == "*"
    amino_str = amino_str[:-1]
    assert "*" not in amino_str
    assert str(amino_seq) == amino_str
    assert str(amino_seq) == amino_str
    assert (len(str(nucl_seq)) % 3) == 0
    return nucl_seq[:-3], amino_seq
