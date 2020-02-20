from typing import Union, TypeVar
from math import log

from nmm.alphabet import Alphabet, BaseAlphabet
from nmm.state import CodonState, MuteState, NormalState
from nmm.sequence import SequenceABC, Sequence
from nmm.path import Step, Path
from nmm.codon import Codon
from nmm.prob import AminoTable, lprob_zero
from nmm import GeneticCode

from .fragment import Fragment
from .amino import AminoFragment, AminoStep, AminoPath

CodonStep = Step[Union[CodonState, MuteState]]
CodonPath = Path[CodonStep]


class CodonFragment(Fragment[BaseAlphabet, CodonState]):
    def __init__(
        self,
        sequence: SequenceABC[BaseAlphabet],
        path: Path[Step[CodonState]],
        homologous: bool,
    ):
        super().__init__(sequence, path, homologous)

    # @property
    # def sequence(self) -> SequenceABC:
    #     from nmm import Sequence as Seq

    #     codon_seq = b"".join(c.symbols for c in self._codons)
    #     seq = Seq(codon_seq, self._path[0].state.alphabet)
    #     return seq

    # def items(self) -> Iterator[Tuple[bytes, CodonStep]]:
    #     # start = end = 0
    #     idx = 0
    #     for step in self._path:
    #         # end += step.seq_len
    #         if isinstance(step.state, MuteState):
    #             yield (b"", step)
    #         else:
    #             assert isinstance(step.state, CodonState)
    #             yield (self._codons[idx].symbols, step)
    #             idx += 1
    #         # start = end

    def decode(self, genetic_code: GeneticCode) -> AminoFragment:
        aminos: List[bytes] = []
        steps = []
        # npath = AminoAcidPath()

        amino_acids = genetic_code.amino_acids()
        amino_abc = genetic_code.amino_alphabet
        base_abc = genetic_code.base_alphabet
        # aa_alphabet = Alphabet(b"".join(amino_acids))

        seq = self.sequence
        start = 0
        for step in self._path:
            if isinstance(step.state, MuteState):

                mstate = MuteState(step.state.name, amino_abc)
                steps.append(AminoStep.create(mstate, 0))

            elif isinstance(step.state, CodonState):

                subseq = seq[start : start + step.seq_len]
                aa = genetic_code.amino_acid(Codon.create(bytes(subseq), base_abc))
                aminos.append(aa)

                lprobs = [lprob_zero()] * amino_abc.length
                lprobs[amino_abc.symbol_idx(aa)] = log(1.0)
                name = step.state.name
                nstate = NormalState(name, amino_abc, lprobs)
                steps.append(AminoStep.create(nstate, 1))

            else:
                raise RuntimeError("Wrong state type.")

            start += step.seq_len

        sequence = Sequence.create(b"".join(aminos), amino_abc)
        return AminoFragment(sequence, AminoPath.create(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
