from collections import defaultdict as dd
from symbol import Sym, ε, φ
from typing import Set, Tuple

from fsa import FSA
from semiring import Semiring
from state import State


class Transformer:
    @staticmethod
    def _add_trim_arcs(F: FSA, T: FSA, AC: Set[State]):
        for i in AC:
            for a, j, w in F.arcs(i):
                if j in AC:
                    T.add_arc(i, a, j, w)

    @staticmethod
    def trim(F: FSA) -> FSA:
        """trims the machine"""

        # compute accessible and co-accessible arcs
        A, C = F.accessible(), F.coaccessible()
        AC = A.intersection(C)

        # create a new F with only the pruned arcs
        T = F.spawn()
        Transformer._add_trim_arcs(F, T, AC)

        # add initial state
        for q, w in F.I:
            if q in AC:
                T.set_I(q, w)

        # add final state
        for q, w in F.F:
            if q in AC:
                T.set_F(q, w)

        return T

    @staticmethod
    def partition(fsa, partition_symbol: Sym = ε) -> Tuple[FSA, FSA]:
        """Partition FSA into two
        (one with arcs of the partition symbol and one with all others)

        Args:
            fsa (FSA): The input FSA
            partition_symbol (Sym, optional): The symbol based on which to
            partition the input FSA

        Returns:
            Tuple[FSA, FSA]: The FSA with non-partition symbol arcs
                             and the FSA with only the partition symbol arcs
        """

        E = fsa.spawn()
        N = fsa.spawn(keep_init=True, keep_final=True)

        for q in fsa.Q:
            E.add_state(q)
            N.add_state(q)

        for i in fsa.Q:
            for a, j, w in fsa.arcs(i):
                if a == partition_symbol:
                    E.add_arc(i, a, j, w)
                else:
                    N.add_arc(i, a, j, w)

        return N, E

    @staticmethod
    def _phi_transitive_closure(
        fsa: FSA,
        q: State,
        single_hops: dd[State, dd[Sym, dd[State, Semiring]]],
    ) -> None:
        """Determines the one-hop arcs possible from a given state `q` according to
        the original transition function and the failure arcs.

        Args:
            fsa (FSA): The input FSA
            q (State): The current state
            single_hops (dd[State, dd[Sym, dd[State, Semiring]]]):
                The dictionary of one-hop arcs found so far in the backward traversal of
                the failure arcs.
        """

        # Keep the list of all taken transitions, and the failure target for `q`
        valid_transitions, fallback_target = set(), None
        for a, t, w in fsa.arcs(q):
            if a == φ:
                fallback_target = t
            else:
                # If the transition is not φ, add it directly to the expanded FSA
                valid_transitions.add(a)
                single_hops[q][a][t] = w

        # If `q` has no fallback target, we are done
        if fallback_target is None:
            return

        # If there is an outgoing failure arc,
        # add the non-existing outgoing transitions
        # according to the fallback target
        for a in single_hops[fallback_target]:
            if a in (fsa.Sigma - valid_transitions - set([φ])):
                for t, w in single_hops[fallback_target][a].items():
                    single_hops[q][a][t] = w

    @staticmethod
    def _failure_expanded_fsa(
        fsa: FSA,
        single_hops: dd[State, dd[Sym, dd[State, Semiring]]],
    ) -> FSA:
        """Construct the expanded FSA based on the one-hop arcs
        found by traversing the failure arcs.

        Args:
            fsa (FSA): The input FSA
            single_hops (dd[State, dd[Sym, dd[State, Semiring]]]):
                The one-hop transitions possible from each state based on the original
                transition function and the failure arcs

        Returns:
            FSA: The failure-arc expanded FSA
        """

        Aʼ = FSA(fsa.R)
        for q in fsa.Q:
            for a in single_hops[q]:
                for t, w in single_hops[q][a].items():
                    Aʼ.add_arc(q, a, t, w)

        for q, w in fsa.I:
            Aʼ.add_I(q, w)
        for q, w in fsa.F:
            Aʼ.add_F(q, w)

        return Aʼ

    @staticmethod
    def expand_phi_arcs(fsa: FSA) -> FSA:
        """Generates an equivalent FSA without failure (φ) transitions
           by creating the missing transitions to the fallback state.
           We assume all failure transitions have weight 1,
           that there is at most 1 outgoing failure arc per state,
           and that they form an *acyclic* subgraph in the FSA.
        Args:
            fsa (FSA): The input FSA with failure arcs according to the
                specifications above.

        Returns:
            FSA: The FSA with the failure arcs expanded.
        """

        _, Aφ = Transformer.partition(fsa, φ)

        # Hacky: preprocess the φ-only FSA to have initial states
        Iφ = set(Aφ.Q)
        for q in Aφ:
            for _, t, _ in Aφ.arcs(q):
                # Remove the target state from the states with no incoming edges
                Iφ.discard(t)

        for q in Iφ:
            Aφ.add_I(q, fsa.R.one)

        # We assume that E is *acyclic*
        assert Aφ.acyclic

        # Keep the map of the transitions accessible using failure arcs
        # from all the states in the FSA
        hops = dd(lambda: dd(lambda: dd(lambda: fsa.R.zero)))

        for q in Aφ.toposort(rev=True):
            Transformer._phi_transitive_closure(fsa, q, hops)

        # After determining all one-hop transitions considering the failure arcs,
        # create the new expanded FSA
        Aʼ = Transformer._failure_expanded_fsa(fsa, hops)

        return Aʼ
