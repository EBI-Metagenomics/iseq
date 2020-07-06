from typing import Dict, List, Set, Union

from nmm import AminoAlphabet, DNAAlphabet, RNAAlphabet, Codon
from .gencode import GeneticCode, get_genetic_code

__all__ = ["CodonTable"]


class CodonTable:
    """
    Codon table.

    Parameters
    ----------
    base_abc
        Base alphabet.
    amino_abc
        Amino acid alphabet.
    name
        NCBI `translation table name`_. Defaults to `"Standard"`.

    .. _translation table name: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """

    def __init__(
        self,
        base_abc: Union[DNAAlphabet, RNAAlphabet],
        amino_abc: AminoAlphabet,
        gencode: GeneticCode = get_genetic_code("Standard"),
    ):

        self._base_alphabet = base_abc
        self._amino_alphabet = amino_abc

        table = translation_table(gencode)

        self._map: Dict[bytes, List[Codon]] = {aa: [] for aa in table}

        for aa, triplets in table.items():
            gcode = self._map[aa]
            for triplet in triplets:
                if isinstance(base_abc, RNAAlphabet):
                    triplet = triplet.replace(b"T", b"U")
                gcode.append(Codon.create(triplet, base_abc))

        self._amino_acid: Dict[Codon, bytes] = {}
        for aa, codons in self._map.items():
            for codon in codons:
                self._amino_acid[codon] = aa

    def codons(self, amino_acid: bytes) -> List[Codon]:
        amino_acid = amino_acid.upper()
        return self._map.get(amino_acid, [])

    def codons_prob(self, amino_acid: bytes) -> Dict[Codon, float]:
        codons = self.codons(amino_acid)
        n = len(codons)
        if n == 0:
            return {}
        return {codon: 1 / n for codon in codons}

    def amino_acid(self, codon: Codon) -> bytes:
        return self._amino_acid[codon]

    def amino_acids(self) -> Set[bytes]:
        return set(self._map.keys())

    @property
    def base_alphabet(self) -> Union[DNAAlphabet, RNAAlphabet]:
        return self._base_alphabet

    @property
    def amino_alphabet(self) -> AminoAlphabet:
        return self._amino_alphabet


def translation_table(gencode: GeneticCode):
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PendingDeprecationWarning)
        from Bio.Data.CodonTable import unambiguous_dna_by_id, unambiguous_dna_by_name

    if gencode.id not in unambiguous_dna_by_id:
        names = str(list(unambiguous_dna_by_name.keys()))
        msg = f"Unknown translation table {gencode.name}. Possible names are: {names}."
        raise ValueError(msg)

    table = unambiguous_dna_by_id[gencode.id]

    btable: Dict[bytes, List[bytes]] = {}

    for codon, aa in table.forward_table.items():
        baa = aa.encode()
        if baa not in btable:
            btable[baa] = []
        btable[baa].append(codon.encode())

    return btable


# aa_3code = {
#     "Ala": "A",
#     "Arg": "R",
#     "Asn": "N",
#     "Asp": "D",
#     "Asx": "B",
#     "Cys": "C",
#     "End": "*",
#     "Gln": "Q",
#     "Glu": "E",
#     "Glx": "Z",
#     "Gly": "G",
#     "His": "H",
#     "Ile": "I",
#     "Leu": "L",
#     "Lys": "K",
#     "Met": "M",
#     "Phe": "F",
#     "Pro": "P",
#     "Ser": "S",
#     "Thr": "T",
#     "Trp": "W",
#     "Tyr": "Y",
#     "Val": "V",
# }


def get_codon_probs(name: str = "uniform"):
    import iseq
    from io import StringIO
    import re
    import requests
    import nmm
    from pandas import read_csv
    import importlib_resources as pkg_resources

    buffer = pkg_resources.open_binary(iseq, "species.table.2007.09.12.tsv.xz")

    df = read_csv(
        buffer,
        sep="\t",
        header=None,
        names=["Name", "Taxonomy ID"],
        dtype={"Name": str, "Taxonomy ID": str},
        skiprows=[0],
        compression="xz",
    )

    row = df[df["Name"] == name]

    if len(row) == 0:
        msg = "Unknown specie name. Please, check http://www.kazusa.or.jp/codon/"
        msg += " for the recognized names."
        raise ValueError(msg)

    taxid = row["Taxonomy ID"].values[0]

    url = f"http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={taxid}&aa=1&style=GCG"
    response = requests.get(url)
    response.raise_for_status()

    match = re.search(r"<PRE>(.+?)</PRE>", response.text, re.MULTILINE | re.DOTALL)
    if match is None:
        raise RuntimeError(f"Could not find table at {url}.")

    txt = match.group(1)
    txt = txt.strip().replace("\n \n", "\n").replace("   ..", "")
    txt = re.sub(" +", " ", txt).replace(" ", "\t")

    df = read_csv(StringIO(txt), header=0, sep="\t", dtype={"Number": int})
    df["Prob"] = df["Number"] / df["Number"].sum()
    df = df.drop(columns=["Number", "/1000", "Fraction"])

    # with open("/Users/horta/tmp.tsv", "w") as file:
    #     file.write(txt)

    pass
    # return dict(zip(df["Name"].values, df["Taxonomy ID"].values))
