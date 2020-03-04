from dataclasses import dataclass
from math import log
from typing import Dict, Generic, List, Tuple, Union

from nmm import HMM, CData
from nmm.path import Path, Step
from nmm.prob import lprob_zero, lprob_is_zero
from nmm.sequence import Sequence
from nmm.state import MuteState, State

from ._typing import MutableState, TState, MutableResults

__all__ = ["AltModel", "Node", "NullModel", "SpecialNode", "Transitions"]


@dataclass
class Transitions:
    MM: float = lprob_zero()
    MI: float = lprob_zero()
    MD: float = lprob_zero()
    IM: float = lprob_zero()
    II: float = lprob_zero()
    DM: float = lprob_zero()
    DD: float = lprob_zero()

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


class Node(Generic[TState]):
    def __init__(self, M: TState, I: TState, D: MuteState):
        self._M = M
        self._I = I
        self._D = D

    @property
    def M(self) -> TState:
        return self._M

    @property
    def I(self) -> TState:
        return self._I

    @property
    def D(self) -> MuteState:
        return self._D

    def states(self) -> List[MutableState[TState]]:
        return [self._M, self._I, self._D]


class SpecialNode(Generic[TState]):
    def __init__(
        self,
        S: MuteState,
        N: TState,
        B: MuteState,
        E: MuteState,
        J: TState,
        C: TState,
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
    def N(self) -> TState:
        return self._N

    @property
    def B(self) -> MuteState:
        return self._B

    @property
    def E(self) -> MuteState:
        return self._E

    @property
    def J(self) -> TState:
        return self._J

    @property
    def C(self) -> TState:
        return self._C

    @property
    def T(self) -> MuteState:
        return self._T

    def states(self) -> List[MutableState[TState]]:
        return [self._S, self._N, self._B, self._E, self._J, self._C, self._T]


class NullModel(Generic[TState]):
    def __init__(
        self, state: TState, special_trans: SpecialTransitions,
    ):
        self._hmm = HMM(state.alphabet)
        self._hmm.add_state(state, 0.0)
        self._state = state
        self._special_transitions = special_trans

    @property
    def state(self) -> TState:
        return self._state

    def set_transition(self, lprob: float):
        self._hmm.set_transition(self.state, self.state, lprob)

    def likelihood(self, sequence: Sequence):
        steps = [Step.create(self.state, 1) for i in range(len(sequence))]
        path = Path.create(steps)
        return self._hmm.likelihood(sequence, path)

    def update_special_transitions(self):
        self.set_transition(self._special_transitions.RR)


class AltModel(Generic[TState]):
    def __init__(
        self,
        special_node: SpecialNode,
        core_nodes: List[Node],
        core_trans: List[Transitions],
        special_trans: SpecialTransitions,
    ):
        self._special_node = special_node
        self._core_nodes = core_nodes
        self._states: Dict[CData, MutableState[TState]] = {}
        self._special_transitions = special_trans
        self._log_occ, self._logZ = _calculate_occupancy(core_trans)

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
        self._hmm = hmm

        for node in core_nodes:
            hmm.add_state(node.M)
            hmm.add_state(node.I)
            hmm.add_state(node.D)

        M1 = core_nodes[0].M
        hmm.set_transition(special_node.B, M1, core_trans[0].MM)

        for i, trans in enumerate(core_trans[1:-1]):
            prev = core_nodes[i]
            next = core_nodes[i + 1]
            hmm.set_transition(prev.M, prev.I, trans.MI)
            hmm.set_transition(prev.I, prev.I, trans.II)
            hmm.set_transition(prev.M, next.M, trans.MM)
            hmm.set_transition(prev.I, next.M, trans.IM)
            hmm.set_transition(prev.M, next.D, trans.MD)
            hmm.set_transition(prev.D, next.D, trans.DD)
            hmm.set_transition(prev.D, next.M, trans.DM)

        Mm = core_nodes[-1].M
        hmm.set_transition(Mm, special_node.E, core_trans[-1].MM)

        # self.view(core_model_only=True)

    def view_emissions(self):
        abc = self._core_nodes[0].M.alphabet
        print("V", self._core_nodes[0].M.lprob(Sequence.create(b"V", abc)))
        print("E", self._core_nodes[1].M.lprob(Sequence.create(b"E", abc)))
        print("V", self._core_nodes[2].M.lprob(Sequence.create(b"V", abc)))
        print("E", self._core_nodes[3].I.lprob(Sequence.create(b"E", abc)))
        print("V", self._core_nodes[4].M.lprob(Sequence.create(b"V", abc)))
        print("Y", self._core_nodes[5].M.lprob(Sequence.create(b"Y", abc)))

    def view(self, core_model_only=False):
        from graphviz import Digraph

        dot = Digraph(comment="Profile")
        dot.attr(rankdir="LR")

        dot.node(self._special_node.B.name.decode())

        if not core_model_only:
            states = [self._special_node.S, self._special_node.N]
            for state in states:
                dot.node(state.name.decode())

        rank = 0
        for core in self._core_nodes:
            with dot.subgraph() as cluster:
                cluster.attr(rank=str(rank))
                cluster.node(core.M.name.decode(), shape="square")
                cluster.node(core.D.name.decode(), shape="circle")
                cluster.node(core.I.name.decode(), shape="diamond")
            rank += 1

        dot.node(self._special_node.E.name.decode())

        if not core_model_only:
            states = [self._special_node.C, self._special_node.T]
            for state in states:
                dot.node(state.name.decode())

            state = self._special_node.J
            dot.node(state.name.decode())

        for s0 in self._states.values():
            for s1 in self._states.values():
                t = self._hmm.transition(s0, s1)
                if lprob_is_zero(t):
                    continue
                label = f"{t:.4f}"
                special = self._is_special_state(s0) or self._is_special_state(s1)
                core_trans = self._is_core_state(s0) and self._is_core_state(s1)
                if core_model_only and not core_trans:
                    continue

                name0 = s0.name.decode()
                name1 = s1.name.decode()
                if special:
                    color = "lightgrey"
                else:
                    color = "black"
                dot.edge(name0, name1, label=label, color=color)

        dot.view()

    def _is_special_state(self, state: MutableState[TState]) -> bool:
        return state in self._special_node.states()

    def _is_core_state(self, state: MutableState[TState]) -> bool:
        B = self._special_node.B
        E = self._special_node.E
        return state not in self._special_node.states() or state in [B, E]

    def set_transition(self, a: State, b: State, lprob: float):
        self._hmm.set_transition(a, b, lprob)

    def core_nodes(self) -> List[Node]:
        return self._core_nodes

    @property
    def special_node(self) -> SpecialNode:
        return self._special_node

    @property
    def length(self) -> int:
        return len(self._core_nodes)

    def viterbi(self, seq: Sequence, window_length: int = 0) -> MutableResults[TState]:
        return self._hmm.viterbi(seq, self.special_node.T, window_length)

    def set_fragment_length(self):
        if self.length == 0:
            return

        B = self.special_node.B
        E = self.special_node.E

        # Uniform local alignment fragment length distribution
        t = self._special_transitions
        t.BM = log(2) - log(self.length) - log(self.length + 1)
        t.ME = 0.0
        for node in self.core_nodes():
            self.set_transition(B, node.M, t.BM)
            self.set_transition(node.M, E, t.ME)

        for node in self.core_nodes()[1:]:
            self.set_transition(node.D, E, 0.0)

    def set_entry_transitions(self):
        log_odds = [locc - self._logZ for locc in self._log_occ]

        B = self.special_node.B
        for lodd, node in zip(log_odds, self.core_nodes()):
            self.set_transition(B, node.M, lodd)

    def set_exit_transitions(self):
        if self.length == 0:
            return

        E = self.special_node.E

        for node in self.core_nodes():
            self.set_transition(node.M, E, 0.0)

        for node in self.core_nodes()[1:]:
            self.set_transition(node.D, E, 0.0)

    def update_special_transitions(self, hmmer3=False):
        t = self._special_transitions
        node = self.special_node

        # TODO: THIS IS FOR HMMER3 VITERBI FILTER COMPARISON
        if hmmer3:
            t.NN = 0.0
            t.CC = 0.0
            t.JJ = 0.0

        self.set_transition(node.S, node.B, t.NB)
        self.set_transition(node.S, node.N, t.NN)
        self.set_transition(node.N, node.N, t.NN)
        self.set_transition(node.N, node.B, t.NB)

        self.set_transition(node.E, node.T, t.EC + t.CT)
        self.set_transition(node.E, node.C, t.EC + t.CC)
        self.set_transition(node.C, node.C, t.CC)
        self.set_transition(node.C, node.T, t.CT)

        self.set_transition(node.E, node.B, t.EJ + t.JB)
        self.set_transition(node.E, node.J, t.EJ + t.JJ)
        self.set_transition(node.J, node.J, t.JJ)
        self.set_transition(node.J, node.B, t.JB)

        self.set_transition(node.B, self._core_nodes[1].D, lprob_zero())
        self.set_transition(node.B, self._core_nodes[0].I, lprob_zero())


def _calculate_occupancy(core_trans: List[Transitions],) -> Tuple[List[float], float]:
    from numpy import logaddexp
    from math import exp
    import sys

    log_occ = [logaddexp(core_trans[0].MI, core_trans[0].MM)]
    for trans in core_trans[1:-1]:
        val0 = log_occ[-1] + logaddexp(trans.MM, trans.MI)
        val1 = log1_p(log_occ[-1]) + trans.DM
        log_occ.append(logaddexp(val0, val1))

    print([exp(i) for i in log_occ])

    logZ = lprob_zero()
    for i, locc in enumerate(log_occ):
        logZ = logaddexp(logZ, locc + log(len(log_occ) - i))

    print("Z", exp(logZ))

    return log_occ, logZ


def log1_p(log_p: float):
    """
    Computes log(1 - p) given log(p).
    """
    from scipy.special import logsumexp

    return logsumexp([0.0, log_p], b=[1.0, -1.0])


class MSVModel(Generic[TState]):
    def __init__(
        self,
        special_node: SpecialNode,
        nodes_trans: List[Tuple[Node, Transitions]],
        special_trans: SpecialTransitions,
    ):
        self._special_node = special_node
        self._core_nodes = [nt[0] for nt in nodes_trans]
        self._states: Dict[CData, Union[TState, MuteState]] = {}
        self._special_transitions = special_trans

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

        if len(nodes_trans) > 0:
            node = nodes_trans[0][0]
            hmm.add_state(node.M)
            prev = node

            for node, _ in nodes_trans[1:]:
                hmm.add_state(node.M)

                hmm.set_transition(prev.M, node.M, 0.0)
                prev = node

        self._hmm = hmm

    def set_transition(self, a: State, b: State, lprob: float):
        self._hmm.set_transition(a, b, lprob)

    def core_nodes(self) -> List[Node]:
        return self._core_nodes

    @property
    def special_node(self) -> SpecialNode:
        return self._special_node

    @property
    def length(self) -> int:
        return len(self._core_nodes)

    def viterbi(self, seq: Sequence, window_length: int = 0) -> MutableResults[TState]:
        return self._hmm.viterbi(seq, self.special_node.T, window_length)

    def set_fragment_length(self):
        if self.length == 0:
            return

        B = self.special_node.B
        E = self.special_node.E

        # Uniform local alignment fragment length distribution
        t = self._special_transitions
        t.BM = log(2) - log(self.length) - log(self.length + 1)
        t.ME = 0.0
        for node in self.core_nodes():
            self.set_transition(B, node.M, t.BM)
            self.set_transition(node.M, E, t.ME)

    def update_special_transitions(self):
        t = self._special_transitions
        node = self.special_node

        self.set_transition(node.S, node.B, t.NB)
        self.set_transition(node.S, node.N, t.NN)
        self.set_transition(node.N, node.N, t.NN)
        self.set_transition(node.N, node.B, t.NB)

        self.set_transition(node.E, node.T, t.EC + t.CT)
        self.set_transition(node.E, node.C, t.EC + t.CC)
        self.set_transition(node.C, node.C, t.CC)
        self.set_transition(node.C, node.T, t.CT)

        self.set_transition(node.E, node.B, t.EJ + t.JB)
        self.set_transition(node.E, node.J, t.EJ + t.JJ)
        self.set_transition(node.J, node.J, t.JJ)
        self.set_transition(node.J, node.B, t.JB)
