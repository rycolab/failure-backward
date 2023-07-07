from collections.abc import Sequence
from symbol import Sym, to_sym, ε
from typing import List, Set, Union


class Alphabet(Sequence):
    def __init__(self, symbols: Union[List[Sym], Set[Sym]]):
        assert len(symbols) > 0, "Alphabet must be non-empty."
        assert ε not in symbols, "ε must not be in the alphabet."
        # We want to have this a set so that the order is deterministic.
        self.symbols = symbols if isinstance(symbols, list) else list(symbols)

    def __str__(self):
        return "{" + ", ".join(str(sym) for sym in self.symbols) + "}"

    def __repr__(self):
        return "{" + ", ".join(str(sym) for sym in self.symbols) + "}"

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return isinstance(other, Alphabet) and self.symbols == other.symbols

    def __len__(self):
        return len(self.symbols)

    def __contains__(self, item):
        return item in self.symbols

    def __iter__(self):
        return iter(self.symbols)

    def __getitem__(self, index):
        return self.symbols[index]

    def __add__(self, other):
        return Alphabet(self.symbols + other.symbols)

    def add(self, symbol: Sym):
        self.symbols.append(symbol)

    def discard(self, symbol: Sym):
        if symbol in self.symbols:
            self.symbols.remove(symbol)

    def union(self, other: "Alphabet"):
        return Alphabet(self.symbols + other.symbols)


def to_alphabet(
    symbols: Union[List[Union[Sym, str]], Set[Union[Sym, str]], str]
) -> Alphabet:
    """Creates an Alphabet from a list of symbols.

    Args:
        symbols (Union[List[Union[Sym, str]], Set[Union[Sym, str]]]): The symbols
            which will form the alphabet.

    Returns:
        Alphabet: The alphabet.
    """

    if isinstance(symbols, Alphabet):
        return symbols
    else:
        return Alphabet(
            [to_sym(sym) if isinstance(sym, str) else sym for sym in symbols]
        )
