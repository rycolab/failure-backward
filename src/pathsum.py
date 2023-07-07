from collections import defaultdict
from symbol import Sym, φ
from typing import DefaultDict, Dict, Set, Tuple

from aggregator import Aggregator
from fsa import FSA
from semiring import Semiring
from state import State


class Strategy:
    VITERBI = 1
    FAILURE_MEMORIZATION = 2
    FAILURE_RING = 3
    FAILURE_GENERAL = 4


class FailureBackward:
    def __init__(self, A: FSA):
        self.A = A
        self.R = A.R
        self._β = dict()
        self.Σ = tuple(self.A.Sigma)
        self._β = defaultdict(lambda: self.R.zero)
        self._β_S = defaultdict(lambda: defaultdict(lambda: self.R.zero))

    def β(self, q: State):
        if self._β[q] == self.R.zero:
            self._β[q] = self.β_Σ(q) + self.A.ρ[q]
        return self._β[q]

    def β_S(self, q: State, S: Set[Sym]):
        if self._β_S[q][S] == self.R.zero:
            self._β_S[q][S] = sum(list(self.β_a(q, a) for a in S), self.R.zero)
        return self._β_S[q][S]

    def β_Σ(self, q: State):  # Line 5
        Σ_q = tuple(self.A.out_symbols(q, ignore_phi=True))
        if φ not in list(self.A.out_symbols(q)):
            return self.β_S(q, Σ_q)
        return (
            self.β_S(self.A.qφ(q), self.Σ)
            + self.β_S(q, Σ_q)
            - self.β_S(self.A.qφ(q), Σ_q)
        )  # Line 9

    def β_a(self, q: State, a: Sym):  # Line 10
        if a in set(self.A.out_symbols(q, ignore_phi=True)):
            r = self.A.R.zero
            for qʼ, w in self.A.a_out_arcs(q, a):
                r += w * self._β[qʼ]
            return r
        elif φ in set(self.A.out_symbols(q)):
            return self.β_a(self.A.qφ(q), a)
        else:
            return self.A.R.zero

    def ring_compute(self) -> Dict[State, Semiring]:
        for q in self.A.toposort(rev=True):
            self._β[q] = self.β(q)
        return self._β


class Pathsum:
    def __init__(self, fsa):
        # basic FSA stuff
        self.fsa = fsa
        self.R = fsa.R
        self.N = self.fsa.num_states

    def pathsum(self, strategy: int) -> Semiring:
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"

        if strategy == Strategy.VITERBI:
            return self.viterbi_pathsum()

        elif strategy == Strategy.FAILURE_MEMORIZATION:
            return self.memorization_pathsum()

        elif strategy == Strategy.FAILURE_RING:
            return self.failure_ring_pathsum()

        elif strategy == Strategy.FAILURE_GENERAL:
            return self.general_failure_pathsum()

        else:
            raise NotImplementedError

    def backward(self, strategy: int) -> DefaultDict[State, Semiring]:
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"

        if strategy == Strategy.VITERBI:
            return self.viterbi_backward()

        if strategy == Strategy.FAILURE_MEMORIZATION:
            return self.memorization_backward()

        elif strategy == Strategy.FAILURE_RING:
            return self.failure_ring_backward()

        elif strategy == Strategy.FAILURE_GENERAL:
            return self.general_failure_backward()

        else:
            raise NotImplementedError

    def viterbi_pathsum(self):
        pathsum = self.R.zero
        β = self.viterbi_backward()
        for q in self.fsa.Q:
            pathsum += self.fsa.λ[q] * β[q]
        return pathsum

    def memorization_pathsum(self):
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
        pathsum = self.R.zero
        β = self.memorization_backward()
        for q in self.fsa.Q:
            pathsum += self.fsa.λ[q] * β[q]
        return pathsum

    def failure_ring_pathsum(self):
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
        pathsum = self.R.zero
        β = self.failure_ring_backward()
        for q in self.fsa.Q:
            pathsum += self.fsa.λ[q] * β[q]
        return pathsum

    def general_failure_pathsum(self):
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
        pathsum = self.R.zero
        β = self.general_failure_backward()
        for q in self.fsa.Q:
            pathsum += self.fsa.λ[q] * β[q]
        return pathsum

    def viterbi_backward(self) -> DefaultDict[State, Semiring]:
        """The Viterbi algorithm run backwards"""

        assert self.fsa.acyclic

        # chart
        β = self.R.chart()

        # base case (paths of length 0)
        for q, w in self.fsa.F:
            β[q] = w

        # recursion
        for p in self.fsa.toposort(rev=True):
            for _, q, w in self.fsa.arcs(p):
                β[p] += w * β[q]

        return β

    def memorization_backward(self) -> DefaultDict[State, Semiring]:
        """Svete et al. (2022) Algorithm 3"""
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"
        A = self.fsa

        β_q = defaultdict(lambda: self.R.zero)
        β_qa = defaultdict(lambda: self.R.zero)
        β_φ = defaultdict(lambda: self.R.zero)

        for q in A.toposort(rev=True):  # Line 2
            for a in A.out_symbols(q, ignore_phi=True):  # Line 3
                for qʼ, w in A.a_out_arcs(q, a):
                    β_qa[(q, a)] += w * β_q[qʼ]

            if A.has_fallback_state(q):  # Line 6
                # Line 7
                for b in set(A.Sigma) - set(A.out_symbols(q)) - {φ}:
                    β_qa[(q, b)] = β_qa[(A.qφ(q), b)]
                    β_φ[(q, "Σ - Σ(q)")] += β_qa[(A.qφ(q), b)]

            β_φ[(q, "Σ(q)")] = sum(
                list(β_qa[(q, a)] for a in A.out_symbols(q)), self.R.zero
            )
            β_φ[(q, "Σ")] += β_φ[(q, "Σ(q)")] + β_φ[(q, "Σ - Σ(q)")]  # Line 10

            β_q[q] = β_φ[(q, "Σ")] + A.ρ[q]  # Line 11

        return β_q

    def failure_ring_backward(self) -> Dict[State, Semiring]:
        """Svete et al. (2022) Algorithm 4"""
        assert self.fsa.acyclic, "Viterbi requires an acyclic FSA"

        return FailureBackward(self.fsa).ring_compute()

    def _failure_trees(self) -> Tuple[FSA, Dict[State, State]]:
        from transformer import Transformer

        _, Aφ = Transformer.partition(self.fsa, partition_symbol=φ)

        Iφ, Fφ = set(Aφ.Q), set(Aφ.Q)
        for q in Aφ:
            for _, t, _ in Aφ.arcs(q):
                Iφ.discard(t)
                Fφ.discard(q)

        for q in Iφ:
            Aφ.add_I(q, Aφ.R.one)
        for q in Fφ:
            Aφ.add_F(q, Aφ.R.one)

        Ts = dict()
        Aφ_r = Aφ.reverse()
        for q in Fφ:
            U = set(Aφ_r.dfs(Is={q})[1].keys())
            for u in U:
                Ts[u] = q

        return Aφ, Ts

    # TODO: remove noqa: C901
    def general_failure_backward(self) -> Dict[State, Semiring]:  # noqa: C901
        A = self.fsa

        Aφ, Ts = self._failure_trees()
        intervals = Aφ.reverse().dfs(intervals=True)[1]

        γs = dict()
        qs = dict()

        fb = FailureBackward(A)

        def _visit(γ, q):  # Algorithm 5 Line 1
            for a in A.out_symbols(q, ignore_phi=True):
                γ.set(a, fb.β_a(q, a))

        def _leave(γ, q):  # Algorithm 5 Line 4
            γ.undo(len(list(A.out_symbols(q, ignore_phi=True))))

        def _visit_plus(γ, q, qT):  # Line 14
            if q != qT:
                _visit_plus(γ, A.qφ(q), qT)
                _visit(γ, q)  # NEW; you don't visit if q == qT since it's already
                # visited by definition of qT

        def _descendant(q, qʼ) -> bool:  # Footnote 9
            """Tests whether q is a descendant of qʼ"""
            return (
                intervals[q][0] <= intervals[qʼ][0]
                and intervals[q][1] >= intervals[qʼ][1]
            )

        for q in A.toposort(rev=True):  # Line 2
            T = Ts[q]

            if not A.has_fallback_state(q):  # Line 4
                γs[T] = Aggregator(A.R, A.Sigma)
                _visit(γs[T], q)
                qs[T] = q

            while not _descendant(qs[T], q):  # Line 8
                _leave(γs[T], qs[T])
                qs[T] = A.qφ(qs[T])

            _visit_plus(γs[T], q, qs[T])  # Line 10
            qs[T] = q

            fb._β[q] = γs[T].value() + A.ρ[q]  # Line 12

        return fb._β
