from collections import defaultdict as dd
from fractions import Fraction
from math import exp, log


# base code from
# https://github.com/timvieira/hypergraphs/blob/master/hypergraphs/semirings/boolean.py
class Semiring:
    zero: "Semiring"
    one: "Semiring"
    idempotent = False

    def __init__(self, score):
        self.score = score

    @classmethod
    def chart(cls, default=None):
        if default is None:
            default = cls.zero
        return dd(lambda: default)

    def __add__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        return self.score == other.score

    def __hash__(self):
        return hash(self.score)


def lcp(str1, str2):
    """computes the longest common prefix"""
    prefix = ""
    for n in range(min(len(str1), len(str2))):
        if str1[n] == str2[n]:
            prefix += str1[n]
        else:
            break
    return prefix


class String(Semiring):
    def __init__(self, score):
        super().__init__(score)

    def star(self):
        return String.one

    def __add__(self, other):
        if other is self.zero:
            return self
        if self is self.zero:
            return other
        return String(lcp(self.score, other.score))

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return String(self.score + other.score)

    def __truediv__(self, other):
        from rayuela.base.misc import lcp

        prefix = lcp(self.score, other.score)
        return String(self.score[len(prefix) :])

    def __eq__(self, other):
        return self.score == other.score

    def __repr__(self):
        return f"{self.score}"

    def __hash__(self):
        return hash(self.score)


# unique "infinity" string
String.zero = String("∞")
# empty string
String.one = String("")
String.idempotent = False
String.cancellative = False


class Boolean(Semiring):
    def __init__(self, score):
        super().__init__(score)

    def star(self):
        return Boolean.one

    def __add__(self, other):
        return Boolean(self.score or other.score)

    def __mul__(self, other):
        if other.score is self.one:
            return self.score
        if self.score is self.one:
            return other.score
        if other.score is self.zero:
            return self.zero
        if self.score is self.zero:
            return self.zero
        return Boolean(other.score and self.score)

    # TODO: is this correct?
    def __invert__(self):
        return Boolean.one

    def __truediv__(self, other):
        return Boolean.one

    def __eq__(self, other):
        return self.score == other.score

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        return f"{self.score}"

    def __str__(self):
        return str(self.score)

    def __hash__(self):
        return hash(self.score)


Boolean.zero = Boolean(False)
Boolean.one = Boolean(True)
Boolean.idempotent = True
# TODO: check
Boolean.cancellative = True


class MaxPlus(Semiring):
    def __init__(self, score):
        super().__init__(score)

    def star(self):
        return self.one

    def __float__(self):
        return float(self.score)

    def __add__(self, other):
        return MaxPlus(max(self.score, other.score))

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return MaxPlus(self.score + other.score)

    def __invert__(self):
        return MaxPlus(-self.score)

    def __truediv__(self, other):
        return MaxPlus(self.score - other.score)

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        return f"MaxPlus({self.score})"


MaxPlus.zero = MaxPlus(float("-inf"))
MaxPlus.one = MaxPlus(0.0)
MaxPlus.idempotent = True
MaxPlus.superior = True
MaxPlus.cancellative = True


class Tropical(Semiring):
    def __init__(self, score):
        self.score = score

    def star(self):
        return self.one

    def __float__(self):
        return float(self.score)

    def __int__(self):
        return int(self.score)

    def __add__(self, other):
        return Tropical(min(self.score, other.score))

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return Tropical(self.score + other.score)

    def __invert__(self):
        return Tropical(-self.score)

    def __truediv__(self, other):
        return Tropical(self.score - other.score)

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        return f"Tropical({self.score})"

    def __str__(self):
        return str(self.score)


Tropical.zero = Tropical(float("inf"))
Tropical.one = Tropical(0.0)
Tropical.idempotent = True
Tropical.superior = True
Tropical.cancellative = True


class Rational(Semiring):
    def __init__(self, score):
        self.score = Fraction(score)

    def star(self):
        return Rational(Fraction("1") / (Fraction("1") - self.score))

    @classmethod
    @property
    def is_field(self):
        return True

    def __float__(self):
        return float(self.score)

    def __add__(self, other):
        return Rational(self.score + other.score)

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return Rational(self.score * other.score)

    def __invert__(self):
        return Rational(1 / self.score)

    def __truediv__(self, other):
        return Rational(self.score / other.score)

    def __eq__(self, other):
        return abs(float(self.score) - float(other.score)) < 1e-6

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        # return f'Real({self.score})'
        return f"{self.score}"

    # TODO: find out why this wasn't inherited
    def __hash__(self):
        return hash(self.score)


Rational.zero = Rational(Fraction("0"))
Rational.one = Rational(Fraction("1"))
Rational.idempotent = False
Rational.cancellative = True


class Real(Semiring):
    def __init__(self, score):
        # TODO: this is hack to deal with the fact
        # that we have to hash weights
        self.score = score

    def star(self):
        return Real(1.0 / (1.0 - self.score))

    @classmethod
    @property
    def is_field(self):
        return True

    def __float__(self):
        return float(self.score)

    def __add__(self, other):
        return Real(self.score + other.score)

    def __sub__(self, other):
        return Real(self.score - other.score)

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return Real(self.score * other.score)

    def __invert__(self):
        return Real(1.0 / self.score)

    def __pow__(self, other):
        return Real(self.score**other)

    def __truediv__(self, other):
        return Real(self.score / other.score)

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        # return f'Real({self.score})'
        return f"{round(self.score, 3)}"

    def __eq__(self, other):
        return abs(float(self.score) - float(other.score)) < 1e-6

    # TODO: find out why this wasn't inherited
    def __hash__(self):
        return hash(self.score)


Real.zero = Real(0.0)
Real.one = Real(1.0)
Real.idempotent = False
Real.cancellative = True


class Log(Semiring):
    def __init__(self, score):
        # TODO: this is hack to deal with the fact
        # that we have to hash weights
        self.score = score

    def star(self):
        return Log(-log(1 / exp(self.score) - 1) - self.score)

    def __float__(self):
        return float(self.score)

    def __add__(self, other):
        # stolen from https://github.com/timvieira/crf/blob/master/crf/basecrf.py
        if self.score > other.score:
            return Log(self.score + log(exp(other.score - self.score) + 1))
            # return Log(self.score + log(sum(exp(other.score-self.score)).sum()))
        return Log(other.score + log(exp(self.score - other.score + 1)))

    def __mul__(self, other):
        return Log(self.score + other.score)

    def __repr__(self):
        # return f'Real({self.score})'
        return f"{round(self.score, 15)}"

    def __eq__(self, other):
        return abs(float(self.score) - float(other.score)) < 1e-6

    # TODO: find out why this wasn't inherited
    def __hash__(self):
        return hash(self.score)


Log.zero = Log(-float("inf"))
Log.one = Log(0.0)
Log.idempotent = False
Log.cancellative = True


class Integer(Semiring):
    def __init__(self, score):
        # TODO: this is hack to deal with the fact
        # that we have to hash weights
        self.score = score

    def __float__(self):
        return float(self.score)

    def __add__(self, other):
        return Integer(self.score + other.score)

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return Integer(self.score * other.score)

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        return f"Integer({self.score})"

    def __eq__(self, other):
        return float(self.score) == float(other.score)

    def __hash__(self):
        return hash(self.score)


Integer.zero = Integer(0)
Integer.one = Integer(1)
Integer.idempotent = False
Integer.cancellative = True
