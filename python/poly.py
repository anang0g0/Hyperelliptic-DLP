from typing import Generic, List
from itertools import dropwhile, zip_longest

PolynomialBaseType = TypeVar("PolynomialBaseType", bound=Field)

class Polynomial(EuclidianRing, Generic[PolynomialBaseType]):
    def __init__(self: "Polynomial[PolynomialBaseType]", *args: PolynomialBaseType) -> None:
        self.coeffs = list(dropwhile(lambda x: x == x.zero(), args))
        if not self.coeffs:
            self.coeffs = [args[0].zero()]

    @staticmethod
    def _mul_poly(a: List[PolynomialBaseType], b: List[PolynomialBaseType]) -> List[PolynomialBaseType]:
        c = [a[0].zero() for i in range(len(a) + len(b) - 1)]
        for i, aa in enumerate(reversed(a)):
            for j, bb in enumerate(reversed(b)):
                c[i + j] += aa * bb
        return list(reversed(c))

    @staticmethod
    def _div_poly(a: List[PolynomialBaseType], b: List[PolynomialBaseType]) -> Tuple[List[PolynomialBaseType], List[PolynomialBaseType]]:
        q = [a[0].zero() for i in range(max(len(a) - len(b) + 1, 0))]
        r = a[:]
        for i in range(len(q)):
            q[i] = r[i] / b[0]
            for k in range(len(b)):
                r[i + k] -= b[k] * q[i]
        return ([a[0].zero()] + q, r)

    def __repr__(self: "Polynomial") -> str:
        s = " " + " + ".join("%s x ^ %d" % (repr(c), len(self.coeffs) - i - 1) for i, c in enumerate(self.coeffs) if c != c.zero()) + " "
        s = s.replace(" + -", " - ")
        s = s.replace(" x ^ 0 ", " ")
        s = s.replace(" x ^ 1 ", " x ")
        s = s.replace(" 1 x ", " x ")
        return s.strip()

    def __floordiv__(self: "Polynomial[PolynomialBaseType]", x: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        q, r = Polynomial._div_poly(self.coeffs, x.coeffs)
        return Polynomial[PolynomialBaseType](*q)

    def __mod__(self: "Polynomial[PolynomialBaseType]", x: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        q, r = Polynomial._div_poly(self.coeffs, x.coeffs)
        return Polynomial[PolynomialBaseType](*r)

    def __mul__(self: "Polynomial[PolynomialBaseType]", x: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        y = Polynomial._mul_poly(self.coeffs, x.coeffs)
        return Polynomial[PolynomialBaseType](*y)

    def __add__(self: "Polynomial[PolynomialBaseType]", x: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        y = list(a + b for a, b in zip_longest(reversed(self.coeffs), reversed(x.coeffs), fillvalue=self.coeffs[0].zero()))
        y.reverse()
        return Polynomial[PolynomialBaseType](*y)

    def __neg__(self: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        return Polynomial[PolynomialBaseType](*(-a for a in self.coeffs))

    def __eq__(self, x):
        if not isinstance(x, Polynomial):
            return NotImplemented
        return self.coeffs == x.coeffs

    def identity(self: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        return Polynomial[PolynomialBaseType](self.coeffs[0].identity())

    def zero(self: "Polynomial[PolynomialBaseType]") -> "Polynomial[PolynomialBaseType]":
        return Polynomial[PolynomialBaseType](self.coeffs[0].zero())

    def euclid_gcd(self: "Polynomial[PolynomialBaseType]", b: "Polynomial[PolynomialBaseType]") -> Tuple["Polynomial[PolynomialBaseType]", "Polynomial[PolynomialBaseType]", "Polynomial[PolynomialBaseType]"]:
        p, q, r = super().euclid_gcd(b)
        r.coeffs = [x / r.coeffs[0] for x in r.coeffs]
        return (p, q, r)

f = Polynomial[Q](Q(Z(1)), Q(Z(1)))  # x + 1
g = Polynomial[Q](Q(Z(1)), -Q(Z(1))) # x - 1
print(f * g) # -> x ^ 2 - 1

x = Polynomial[Q](Q(Z(2)), Q(Z(1)), Q(Z(-3)), Q(Z(2))) # 2 x ^ 3 + x ^ 2 - 3 x + 2
y = Polynomial[Q](Q(Z(1)), Q(Z(0)), Q(Z(1))) # x ^ 2 + 1
print(x // y) # -> 2 x + 1
print(x % y)  # -> -5 x + 1