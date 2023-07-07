from math import ceil, log2
from symbol import Sym
from typing import Type

from alphabet import Alphabet
from semiring import Semiring


class FenwickTree:
    """Fenwick tree class from Svete et al. (2022)."""

    def __init__(self, R: Type[Semiring], Sigma: Alphabet):
        """
        Args:
            R (type): The semiring to use.
            Sigma (Alphabet): The alphabet of input symbols.
        """
        self.R = R
        self.N = len(Sigma)

        self.D = int(2 ** ceil(log2(self.N)))

        self.sym2idx = {a: self.D + i for i, a in enumerate(Sigma)}
        self.values = [self.R.zero for _ in range(2 * self.D + 1)]

    def update(self, a: Sym, v: Semiring):
        i = self.sym2idx[a]

        self.values[i] = v
        while i > 0:
            i //= 2
            self.values[i] = self.values[2 * i] + self.values[2 * i + 1]

    def __getitem__(self, a: Sym) -> Semiring:
        return self.values[self.sym2idx[a]]


class Aggregator:
    """Aggregator class from Svete et al. (2022)."""

    def __init__(self, R: Type[Semiring], Sigma: Alphabet):
        """
        Args:
            R (type): The semiring to use.
            Sigma (Alphabet): The alphabet of input symbols.
        """
        self.R = R
        self.Sigma = Sigma
        self.tree = FenwickTree(self.R, self.Sigma)

        self.backlog = list()  # Keeps the list of the most recent updates to the
        # aggregator in the form (a, v_old).
        # This is used to undo the updates when the aggregator is reset.

    def set(self, a: Sym, v: Semiring):
        self.backlog.append((a, self.tree[a]))
        self.tree.update(a, v)

    def get(self, a: Sym) -> Semiring:
        return self.tree[a]

    def value(self) -> Semiring:
        return self.tree.values[1]

    def undo(self, n: int):
        for i in range(n):
            a, v_old = self.backlog.pop()
            self.tree.update(a, v_old)
