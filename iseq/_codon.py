from math import log
from typing import List, Tuple, Union

from imm import MuteState, NormalState, Path, Sequence, SequenceABC, Step, lprob_zero
from nmm import BaseAlphabet, Codon, CodonState, GeneticCode

from ._amino import AminoFragment, AminoPath, AminoStep
from ._fragment import Fragment

CodonStep = Step[Union[CodonState, MuteState]]
CodonPath = Path[CodonStep]

__all__ = ["CodonFragment", "CodonStep", "CodonPath"]


class CodonFragment(Fragment[BaseAlphabet, CodonState]):
    def decode(self, genetic_code: GeneticCode) -> AminoFragment:
        aminos: List[bytes] = []
        steps: List[AminoStep] = []

        amino_abc = genetic_code.amino_alphabet

        start = 0
        for step in self.path:
            if isinstance(step.state, MuteState):

                mstate = MuteState.create(step.state.name, amino_abc)
                steps.append(AminoStep.create(mstate, 0))

            elif isinstance(step.state, CodonState):

                subseq = self.sequence[start : start + step.seq_len]
                name = step.state.name

                amino, new_step = _create_step(genetic_code, subseq, name)
                aminos.append(amino)
                steps.append(new_step)

            else:
                raise RuntimeError("Wrong state type.")

            start += step.seq_len

        sequence = Sequence[BaseAlphabet].create(b"".join(aminos), amino_abc)
        return AminoFragment(sequence, AminoPath.create(steps), self.homologous)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


def _create_step(
    gcode: GeneticCode, seq: SequenceABC[BaseAlphabet], state_name: bytes
) -> Tuple[bytes, AminoStep]:
    base_abc = seq.alphabet
    amino = gcode.amino_acid(Codon.create(bytes(seq), base_abc))

    amino_abc = gcode.amino_alphabet
    lprobs = [lprob_zero()] * amino_abc.length
    lprobs[amino_abc.symbol_idx(amino)] = log(1.0)
    state = NormalState.create(state_name, amino_abc, lprobs)

    return amino, AminoStep.create(state, 1)
