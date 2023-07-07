import pytest

from alphabet import to_alphabet
from pathsum import Pathsum, Strategy
from random_fsa import random_machine
from semiring import Real, Tropical
from transformer import Transformer


@pytest.mark.parametrize("R", [Real, Tropical])
def test_memorization_backward(R):
    Sigma = to_alphabet("abcdefg")

    for _ in range(50):
        A = random_machine(
            Sigma=Sigma,
            R=R,
            num_states=20,
            bias=0.2,
            no_eps=True,
            no_phi=False,
            acyclic=True,
            deterministic=False,
            trimmed=True,
            phi_bias=0.5,
        )

        Ae = Transformer.expand_phi_arcs(A)

        assert Pathsum(Ae).pathsum(Strategy.VITERBI) == Pathsum(A).pathsum(
            Strategy.FAILURE_MEMORIZATION
        )


def test_ring_backward():
    Sigma = to_alphabet("abcdefg")

    for _ in range(50):
        A = random_machine(
            Sigma=Sigma,
            R=Real,
            num_states=20,
            bias=0.2,
            no_eps=True,
            no_phi=False,
            acyclic=True,
            deterministic=False,
            trimmed=True,
            phi_bias=0.5,
        )

        Ae = Transformer.expand_phi_arcs(A)

        assert Pathsum(Ae).pathsum(Strategy.VITERBI) == Pathsum(A).pathsum(
            Strategy.FAILURE_RING
        )


@pytest.mark.parametrize("R", [Real, Tropical])
def test_general_backward(R):
    Sigma = to_alphabet("abcdefg")

    for _ in range(50):
        A = random_machine(
            Sigma=Sigma,
            R=R,
            num_states=20,
            bias=0.2,
            no_eps=True,
            no_phi=False,
            acyclic=True,
            deterministic=False,
            trimmed=True,
            phi_bias=0.5,
        )

        Ae = Transformer.expand_phi_arcs(A)

        assert Pathsum(Ae).pathsum(Strategy.VITERBI) == Pathsum(A).pathsum(
            Strategy.FAILURE_GENERAL
        )
