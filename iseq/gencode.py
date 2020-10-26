# NCBI genetic code table version 4.6
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

__all__ = ["GeneticCode"]

genetic_codes: List[Tuple[str, str, int]] = []
names: Dict[str, int] = {}
alt_names: Dict[str, int] = {}
ids: Dict[int, int] = {}


@dataclass
class GeneticCode:
    name: str
    alt_name: str
    id: int

    def __init__(
        self,
        name: Optional[str] = None,
        alt_name: Optional[str] = None,
        id: Optional[int] = None,
    ):
        n = sum([name is None, alt_name is None, id is None])

        if n != 2:
            raise ValueError("You must use one, and only one, parameter.")

        if name is not None:
            if name in names:
                self.name, self.alt_name, self.id = genetic_codes[names[name]]
                return
            raise ValueError(f"Unknown name {name}.")

        if alt_name is not None:
            if alt_name in alt_names:
                self.name, self.alt_name, self.id = genetic_codes[alt_names[alt_name]]
                return
            raise ValueError(f"Unknown alternative name {alt_name}.")

        assert id is not None
        self.name, self.alt_name, self.id = genetic_codes[ids[id]]


def register_ncbi_genetic_code(name: str, alt_name: str, id: int):
    names[name] = len(genetic_codes)
    if alt_name != "":
        alt_names[alt_name] = len(genetic_codes)
    ids[id] = len(genetic_codes)
    genetic_codes.append((name, alt_name, id))


register_ncbi_genetic_code(
    "Standard",
    "SGC0",
    1,
)

register_ncbi_genetic_code(
    "Vertebrate Mitochondrial",
    "SGC1",
    2,
)

register_ncbi_genetic_code(
    "Yeast Mitochondrial",
    "SGC2",
    3,
)

register_ncbi_genetic_code(
    "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate "
    "Mitochondrial; Mycoplasma; Spiroplasma",
    "SGC3",
    4,
)

register_ncbi_genetic_code(
    "Invertebrate Mitochondrial",
    "SGC4",
    5,
)

register_ncbi_genetic_code(
    "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
    "SGC5",
    6,
)

register_ncbi_genetic_code(
    "Echinoderm Mitochondrial; Flatworm Mitochondrial",
    "SGC8",
    9,
)

register_ncbi_genetic_code(
    "Euplotid Nuclear",
    "SGC9",
    10,
)

register_ncbi_genetic_code(
    "Bacterial, Archaeal and Plant Plastid",
    "",
    11,
)

register_ncbi_genetic_code(
    "Alternative Yeast Nuclear",
    "",
    12,
)

register_ncbi_genetic_code(
    "Ascidian Mitochondrial",
    "",
    13,
)

register_ncbi_genetic_code(
    "Alternative Flatworm Mitochondrial",
    "",
    14,
)

register_ncbi_genetic_code(
    "Blepharisma Macronuclear",
    "",
    15,
)

register_ncbi_genetic_code(
    "Chlorophycean Mitochondrial",
    "",
    16,
)

register_ncbi_genetic_code(
    "Trematode Mitochondrial",
    "",
    21,
)

register_ncbi_genetic_code(
    "Scenedesmus obliquus Mitochondrial",
    "",
    22,
)

register_ncbi_genetic_code(
    "Thraustochytrium Mitochondrial",
    "",
    23,
)

register_ncbi_genetic_code(
    "Rhabdopleuridae Mitochondrial",
    "",
    24,
)

register_ncbi_genetic_code(
    "Candidate Division SR1 and Gracilibacteria",
    "",
    25,
)

register_ncbi_genetic_code(
    "Pachysolen tannophilus Nuclear",
    "",
    26,
)

register_ncbi_genetic_code(
    "Karyorelict Nuclear",
    "",
    27,
)

register_ncbi_genetic_code(
    "Condylostoma Nuclear",
    "",
    28,
)

register_ncbi_genetic_code(
    "Mesodinium Nuclear",
    "",
    29,
)

register_ncbi_genetic_code(
    "Peritrich Nuclear",
    "",
    30,
)

register_ncbi_genetic_code(
    "Blastocrithidia Nuclear",
    "",
    31,
)

register_ncbi_genetic_code(
    "Balanophoraceae Plastid",
    "",
    32,
)

register_ncbi_genetic_code(
    "Cephalodiscidae Mitochondrial",
    "",
    33,
)
