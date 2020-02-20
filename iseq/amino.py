from typing import TypeVar, Union

# from typing import List, Sequence, Tuple, TypeVar, Union

from nmm.alphabet import Alphabet, AminoAlphabet
from nmm.path import Path, Step
from nmm.sequence import SequenceABC
from nmm.state import MuteState, NormalState

# from hmmer_reader import HMMERProfile

from .fragment import Fragment

# from .model import AltModel, Node, NullModel, SpecialNode, Transitions
# from .profile import Profile
# from .result import SearchResult

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)

# AminoFragment = Fragment[TAlphabet, Union[NormalState, MuteState]]
AminoStep = Step[Union[NormalState, MuteState]]
AminoPath = Path[AminoStep]
# AminoSearchResult = SearchResult[TAlphabet, Union[NormalState, MuteState]]

# AminoNode = Node[NormalState]
# AminoSpecialNode = SpecialNode[NormalState]
# AminoNullModel = NullModel[NormalState]
# AminoAltModel = AltModel[NormalState]


class AminoFragment(Fragment[AminoAlphabet, NormalState]):
    def __init__(
        self,
        sequence: SequenceABC[AminoAlphabet],
        path: Path[Step[NormalState]],
        homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    # @property
    # def sequence(self) -> SequenceABC:
    #     from nmm import Sequence as Seq

    #     seq = Seq(self._aminos, self._path[0].state.alphabet)
    #     return seq

    # def items(self) -> Iterator[Tuple[bytes, AminoStep]]:
    #     # start = end = 0
    #     idx = 0
    #     for step in self._path:
    #         if isinstance(step.state, MuteState):
    #             # end += step.seq_len
    #             yield (b"", step)
    #         else:
    #             assert isinstance(step.state, NormalState)
    #             yield (self._aminos[idx : idx + 1], step)
    #             idx += 1
    #         # start = end

    # def decode(self, genetic_code: GeneticCode) -> AminoAcidFragment:
    #     nseq: List[bytes] = []
    #     npath = AminoAcidPath()

    #     amino_acids = genetic_code.amino_acids()
    #     aa_alphabet = Alphabet(b"".join(amino_acids))

    #     start: int = 0
    #     i = 0
    #     for step in self._path.steps():
    #         if isinstance(step.state, MuteState):
    #             mstate = MuteState(step.state.name, aa_alphabet)
    #             npath.append_amino_acid_step(mstate, 0)
    #         else:
    #             assert isinstance(step.state, TableState)
    #             aa = genetic_code.amino_acid(self._codons[i])
    #             nseq.append(aa)

    #             nstate = NormalState(step.state.name, aa_alphabet, {aa: LOG1})
    #             npath.append_amino_acid_step(nstate, 1)
    #             i += 1

    #         start += step.seq_len

    #     return AminoAcidFragment(b"".join(nseq), npath, self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
