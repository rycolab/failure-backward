import random
import string
from fractions import Fraction
from math import sqrt
from symbol import Sym, ε, φ
from typing import List, Optional, Tuple, Type, Union

from alphabet import Alphabet
from fsa import FSA, State
from pathsum import Pathsum
from semiring import (
    Boolean,
    Integer,
    MaxPlus,
    Rational,
    Real,
    Semiring,
    String,
    Tropical,
)


def _random_weight(semiring, **kwargs):  # noqa: C901
    if semiring is String:
        str_len = int(random.random() * 8 + 1)
        return semiring(
            "".join(random.choice(string.ascii_lowercase) for _ in range(str_len))
        )

    elif semiring is Boolean:
        return semiring(True)

    elif semiring is Real:
        tol = 1e-3
        s = kwargs.get("divide_by", 2)
        random_weight = round(random.random() / s, 3)
        while random_weight < sqrt(tol):
            random_weight = round(random.random() / s, 3)
        return semiring(random_weight)

    elif semiring is Rational:
        return semiring(Fraction(f"{random.randint(1, 1)}/{random.randint(10, 15)}"))

    elif semiring is Tropical:
        return semiring(random.randint(0, 50))

    elif semiring is Integer:
        return semiring(random.randint(1, 10))

    elif semiring is MaxPlus:
        return semiring(random.randint(-10, -1))


def _add_arc(
    i: int,
    a: Union[Sym, Tuple[Sym, Sym]],
    j: int,
    used_a: List[Union[Sym, Tuple[Sym, Sym]]],
    A: FSA,
    bias: float = 0.25,
    acyclic: bool = False,
    deterministic: bool = True,
    **kwargs,
) -> bool:
    """Handles adding a state to the random machine.

    Args:
        i (int): _description_
        a (Sym): _description_
        j (int): _description_
        used_a (List[Sym]): _description_
        fsa (FSA): _description_
        bias (float, optional): _description_. Defaults to 0.25.
        acyclic (bool, optional): _description_. Defaults to False.
        deterministic (bool, optional): _description_. Defaults to True.
        fst (bool, optional): _description_. Defaults to False.
        kwargs: Arguments for random weight generation

    Returns:
        Whether the arc was added.
    """

    if (deterministic or a == φ) and a in used_a:
        # always add at most one failure arc
        # or at most of any symbol if machine should be deterministic
        return False

    bias = bias if a != φ else kwargs.get("phi_bias", bias)

    if random.random() < bias:
        w = _random_weight(A.R, **kwargs) if a != φ else A.R.one
        if acyclic or a == φ:
            # Make sure that the failure arcs *always* form an acyclic subgraph
            if i < j:
                assert isinstance(a, Sym)
                A.add_arc(State(i), a, State(j), w)
                used_a.append(a)
                return True
            else:
                return False
        else:
            assert isinstance(a, Sym)
            A.add_arc(State(i), a, State(j), w)
            used_a.append(a)
            return True
    else:
        return False


def _random_machine(
    Σ: Alphabet,
    R: Type[Semiring],
    num_states: int,
    bias: float = 0.25,
    no_eps: bool = False,
    no_phi: bool = True,
    acyclic: bool = False,
    deterministic: bool = True,
    **kwargs,
) -> FSA:
    fsa = FSA(R=R)

    if not no_eps:
        Σ.add(ε)
    else:
        Σ.discard(ε)

    if not no_phi:
        Σ.add(φ)
    else:
        Σ.discard(φ)

    for i in range(num_states):
        used_a = []
        for j in random.choices(range(num_states), k=num_states):
            for a in Σ:
                _add_arc(
                    i=i,
                    a=a,
                    j=j,
                    used_a=used_a,
                    A=fsa,
                    bias=bias,
                    acyclic=acyclic,
                    deterministic=deterministic,
                    **kwargs,
                )
    fsa.set_I(State(0), _random_weight(fsa.R, divide_by=1))
    fsa.set_F(State(fsa.num_states - 1), _random_weight(fsa.R, divide_by=1))

    return fsa


def random_machine(
    Sigma: Alphabet,
    R: Type[Semiring],
    num_states: int,
    bias: float = 0.25,
    no_eps: bool = False,
    no_phi: bool = True,
    eigen: bool = False,
    acyclic: bool = False,
    deterministic: bool = True,
    trimmed: bool = True,
    pushed: bool = False,
    fst: bool = False,
    seed: Optional[int] = None,
    **kwargs,
) -> FSA:
    """
    Creates a random WFSA or WFST.
    It takes a number of parameters that control the properties of the machine.

    Args:
        Sigma (Alphabet): The alphabet of the WFSA.
        R (Type[Semiring]): The semiring of the WFSA.
        num_states (int): The number of states of the WFSA.
        bias (float, optional): The probability of realising an edge between
                                any pair of states (q, q') with a specific symbol.
                                Defaults to 0.25.
        no_eps (bool, optional): When true, the WFSA contains no ε transitions.
                                 Defaults to False.
        no_phi (bool, optional): When true, the WFSA contains no φ transitions.
                               Defaults to True.
        eigen (bool, optional): _description_. Defaults to False.
        acyclic (bool, optional): When true, the WFSA will be acyclic by design.
                                  Defaults to False.
        deterministic (bool, optional): When true, the WFSA will be deterministic.
                                        Defaults to True.
        trimmed (bool, optional): When true, trims the machine to make it smaller.
                                  Defaults to True.
        pushed (bool, optional): When true, pushes the machine to make it locally
                                 normalised.
                                 Defaults to False.
        fst (bool, optional): Whether to create a random _transducer_.
                              Defaults to False.
        seed (int, optional): The seed for the random number generator.
        kwargs: Arguments for random weight generation

    Returns:
        FSA: A random WFSA of WFST.
    """

    random.seed(seed)

    fsa = None
    while fsa is None or not fsa.num_states:
        fsa = _random_machine(
            Sigma,
            R,
            num_states,
            bias=bias,
            no_eps=no_eps,
            no_phi=no_phi,
            acyclic=acyclic,
            deterministic=deterministic,
            fst=fst,
            **kwargs,
        )

        # Trim the machine to make it smaller
        if trimmed:
            fsa = fsa.trim()

        if eigen and R is Real:
            pathsum = Pathsum(fsa)
            if pathsum.max_eval() >= 1.0:
                fsa = None

    if pushed:
        fsa = fsa.push()

    return fsa
