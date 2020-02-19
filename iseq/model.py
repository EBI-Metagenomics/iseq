from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Sequence, Tuple, Generic, TypeVar, List, Union, Dict

from nmm import HMM, CData
from nmm.path import Path, Step
from nmm.prob import LPROB_ZERO
from nmm.result import CResults
from nmm.sequence import CSequence
from nmm.state import CState, MuteState

TEmissionState = TypeVar("TEmissionState", bound=CState)


@dataclass
class Transitions:
    MM: float = LPROB_ZERO
    MI: float = LPROB_ZERO
    MD: float = LPROB_ZERO
    IM: float = LPROB_ZERO
    II: float = LPROB_ZERO
    DM: float = LPROB_ZERO
    DD: float = LPROB_ZERO

    def normalize(self):
        from numpy import logaddexp

        m_norm: float = logaddexp(logaddexp(self.MM, self.MI), self.MD)
        self.MM -= m_norm
        self.MI -= m_norm
        self.MD -= m_norm

        i_norm: float = logaddexp(self.IM, self.II)
        self.IM -= i_norm
        self.II -= i_norm

        d_norm: float = logaddexp(self.DM, self.DD)
        self.DM -= d_norm
        self.DD -= d_norm


@dataclass
class SpecialTransitions:
    NN: float = 0.0
    NB: float = 0.0
    EC: float = 0.0
    CC: float = 0.0
    CT: float = 0.0
    EJ: float = 0.0
    JJ: float = 0.0
    JB: float = 0.0
    RR: float = 0.0
    BM: float = 0.0
    ME: float = 0.0


class Node(Generic[TEmissionState]):
    def __init__(self, M: TEmissionState, I: TEmissionState, D: MuteState):
        self._M = M
        self._I = I
        self._D = D

    @property
    def M(self) -> TEmissionState:
        return self._M

    @property
    def I(self) -> TEmissionState:
        return self._I

    @property
    def D(self) -> MuteState:
        return self._D

    def states(self) -> List[Union[TEmissionState, MuteState]]:
        return [self._M, self._I, self._D]


class SpecialNode(Generic[TEmissionState]):
    def __init__(
        self,
        S: MuteState,
        N: TEmissionState,
        B: MuteState,
        E: MuteState,
        J: TEmissionState,
        C: TEmissionState,
        T: MuteState,
    ):
        self._S = S
        self._N = N
        self._B = B
        self._E = E
        self._J = J
        self._C = C
        self._T = T

    @property
    def S(self) -> MuteState:
        return self._S

    @property
    def N(self) -> TEmissionState:
        return self._N

    @property
    def B(self) -> MuteState:
        return self._B

    @property
    def E(self) -> MuteState:
        return self._E

    @property
    def J(self) -> TEmissionState:
        return self._J

    @property
    def C(self) -> TEmissionState:
        return self._C

    @property
    def T(self) -> MuteState:
        return self._T

    def states(self) -> List[Union[TEmissionState, MuteState]]:
        return [self._S, self._N, self._B, self._E, self._J, self._C, self._T]


class NullModel(Generic[TEmissionState]):
    def __init__(self, state: TEmissionState):
        self._hmm = HMM(state.alphabet)
        self._hmm.add_state(state, 0.0)
        self._state = state

    @property
    def state(self) -> TEmissionState:
        return self._state

    def set_transition(self, lprob: float):
        self._hmm.set_transition(self.state, self.state, lprob)

    def likelihood(self, sequence: CSequence):
        path = Path([Step(self.state, 1) for i in range(sequence.length)])
        return self._hmm.likelihood(sequence, path)


class AltModel(Generic[TEmissionState]):
    def __init__(
        self,
        special_node: SpecialNode,
        nodes_trans: Sequence[Tuple[Node, Transitions]],
    ):
        self._special_node = special_node
        self._core_nodes = [nt[0] for nt in nodes_trans]
        self._states: Dict[CData, Union[TEmissionState, MuteState]] = {}

        for node in self._core_nodes:
            for state in node.states():
                self._states[state.imm_state] = state

        for state in special_node.states():
            self._states[state.imm_state] = state

        hmm = HMM(special_node.S.alphabet)
        hmm.add_state(special_node.S, 0.0)
        hmm.add_state(special_node.N)
        hmm.add_state(special_node.B)
        hmm.add_state(special_node.E)
        hmm.add_state(special_node.J)
        hmm.add_state(special_node.C)
        hmm.add_state(special_node.T)

        self._special_transitions = SpecialTransitions()

        if len(nodes_trans) > 0:
            node, trans = nodes_trans[0]
            hmm.add_state(node.M)
            hmm.add_state(node.I)
            hmm.add_state(node.D)
            prev = node

            for node, trans in nodes_trans[1:]:
                hmm.add_state(node.M)
                hmm.add_state(node.I)
                hmm.add_state(node.D)

                hmm.set_transition(prev.M, node.M, trans.MM)
                hmm.set_transition(prev.M, prev.I, trans.MI)
                hmm.set_transition(prev.M, node.D, trans.MD)
                hmm.set_transition(prev.I, node.M, trans.IM)
                hmm.set_transition(prev.I, prev.I, trans.II)
                hmm.set_transition(prev.D, node.M, trans.DM)
                hmm.set_transition(prev.D, node.D, trans.DD)
                prev = node

        self._hmm = hmm

    def set_transition(self, a: CState, b: CState, lprob: float):
        self._hmm.set_transition(a, b, lprob)

    def core_nodes(self) -> Sequence[Node]:
        return self._core_nodes

    @property
    def special_node(self) -> SpecialNode:
        return self._special_node

    @property
    def special_transitions(self) -> SpecialTransitions:
        return self._special_transitions

    @property
    def length(self) -> int:
        return len(self._core_nodes)

    def viterbi(
        self, seq: CSequence, window_length: int = 0
    ) -> Tuple[float, Path[Step[Union[TEmissionState, MuteState]]]]:

        results = self._hmm.viterbi(seq, self.special_node.T, window_length)
        # TODO: implement multiple windows
        assert len(results) == 1

        path = results[0].path
        score = results[0].loglikelihood

        steps = [
            Step(self._states[step.state.imm_state], step.seq_len) for step in path
        ]
        new_path = Path(steps)

        return (score, new_path)
