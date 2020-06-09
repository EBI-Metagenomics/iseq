from typing import List

from imm import HMM

from .model import Node, SpecialNode

__all__ = ["core_nodes", "special_node"]

# def transitions(
#     hmm: HMM, nodes: List[Node], special_node: SpecialNode
# ) -> List[Transitions]:
#     breakpoint()
#     MM = hmm.get_transition(nodes[0].M, nodes[1].M)
#     MI = hmm.get_transition(nodes[0].M, nodes[0].I)
#     MD = hmm.get_transition(nodes[0].M, nodes[1].D)
#     IM = hmm.get_transition(nodes[0].I, nodes[1].M)
#     II = hmm.get_transition(nodes[0].I, nodes[0].I)
#     DM = -0.0
#     # DD = -inf
#     pass


def core_nodes(hmm: HMM) -> List[Node]:
    Ms = GrowingList()
    Is = GrowingList()
    Ds = GrowingList()
    for state in hmm.states().values():
        if state.name[:1] == b"M":
            idx = int(state.name[1:]) - 1
            Ms[idx] = state
        elif state.name[:1] == b"I":
            idx = int(state.name[1:]) - 1
            Is[idx] = state
        elif state.name[:1] == b"D":
            idx = int(state.name[1:]) - 1
            Ds[idx] = state

    return [Node(i[0], i[1], i[2]) for i in zip(Ms, Is, Ds)]


def special_node(hmm: HMM) -> SpecialNode:
    S = hmm.find_state(b"S")
    N = hmm.find_state(b"N")
    B = hmm.find_state(b"B")
    E = hmm.find_state(b"E")
    J = hmm.find_state(b"J")
    C = hmm.find_state(b"C")
    T = hmm.find_state(b"T")
    return SpecialNode(S, N, B, E, J, C, T)


class GrowingList(list):
    def __setitem__(self, index, value):
        if index >= len(self):
            self.extend([None] * (index + 1 - len(self)))
        list.__setitem__(self, index, value)
