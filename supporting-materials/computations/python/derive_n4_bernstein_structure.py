"""Derive the exact polynomial structure used by the ``n=4`` proof.

The four-knot reduction has three polynomial regions and one analytic
AM--GM region.  This script starts from the uniform-simplex cap formula,
performs exact polynomial division by the nonnegative boundary factors,
constructs a fixed tetrahedralization of each remaining polytope, and
finds every identically zero Bernstein coefficient over ``QQ(b,c,d,e)``.

The rational-function calculation is implemented with a sparse polynomial
numerator and a factored denominator.  Consequently a coefficient is
classified as a structural zero only after exact cancellation in a
canonical rational calculation; no floating-point tolerance is involved.
The output contains the lower-degree polynomial formulas, tetrahedra, and
zero-index lists consumed by :mod:`n4_gaffke_certificate`.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sp


REGION_TETRAHEDRA = {
    "A": [["r0", "r1", "r2", "r3"]],
    "B": [
        ["r3", "r2", "r1", "p0"],
        ["r3", "q2", "q3", "p0"],
        ["r3", "q2", "r2", "p0"],
    ],
    "C": [
        ["r3", "q2", "r2", "pc"],
        ["q3", "r3", "q4", "pc"],
        ["q3", "r3", "q2", "pc"],
    ],
}

VERTEX_FORMULAS = {
    "r0": ["0", "0", "0"],
    "r1": ["0", "0", "e/(1-d)"],
    "r2": ["0", "e/(1-c-d)", "e/(1-c-d)"],
    "r3": ["e/(1-b-c-d)"] * 3,
    "p0": ["0", "0", "1"],
    "q2": ["0", "(d+e)/(1-c)", "1"],
    "q3": ["(d+e)/(1-b-c)", "(d+e)/(1-b-c)", "1"],
    "pc": ["0", "1", "1"],
    "q4": ["(c+d+e)/(1-b)", "1", "1"],
}


class Polynomial4:
    """Sparse polynomial over ``QQ`` in ``b,c,d,e``."""

    def __init__(self, coefficients=None):
        self.coefficients = {
            index: Fraction(value)
            for index, value in (coefficients or {}).items()
            if value
        }

    @classmethod
    def constant(cls, value):
        return cls({(0, 0, 0, 0): Fraction(value)})

    @classmethod
    def from_sympy(cls, expression, variables):
        coefficients = {}
        for powers, value in sp.Poly(expression, *variables).terms():
            coefficients[powers] = Fraction(int(value.p), int(value.q))
        return cls(coefficients)

    def __bool__(self):
        return bool(self.coefficients)

    def __add__(self, other):
        if not isinstance(other, Polynomial4):
            other = Polynomial4.constant(other)
        coefficients = dict(self.coefficients)
        for index, value in other.coefficients.items():
            updated = coefficients.get(index, Fraction(0)) + value
            if updated:
                coefficients[index] = updated
            else:
                coefficients.pop(index, None)
        return Polynomial4(coefficients)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial4({index: -value
                            for index, value in self.coefficients.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return Polynomial4.constant(other) - self

    def __mul__(self, other):
        if not isinstance(other, Polynomial4):
            other = Polynomial4.constant(other)
        coefficients = {}
        for left_index, left_value in self.coefficients.items():
            for right_index, right_value in other.coefficients.items():
                index = tuple(a + b for a, b in zip(
                    left_index, right_index))
                updated = (coefficients.get(index, Fraction(0))
                           + left_value * right_value)
                if updated:
                    coefficients[index] = updated
                else:
                    coefficients.pop(index, None)
        return Polynomial4(coefficients)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        result = Polynomial4.constant(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent //= 2
        return result


B4 = Polynomial4({(1, 0, 0, 0): 1})
C4 = Polynomial4({(0, 1, 0, 0): 1})
D4 = Polynomial4({(0, 0, 1, 0): 1})
E4 = Polynomial4({(0, 0, 0, 1): 1})
ONE4 = Polynomial4.constant(1)
DENOMINATOR_BASES = (
    ONE4 - D4,
    ONE4 - C4 - D4,
    ONE4 - B4 - C4 - D4,
    ONE4 - C4,
    ONE4 - B4 - C4,
    ONE4 - B4,
)
DENOMINATOR_POWERS = [[Polynomial4.constant(1)] for _ in DENOMINATOR_BASES]


def _denominator_power(index, exponent):
    powers = DENOMINATOR_POWERS[index]
    while len(powers) <= exponent:
        powers.append(powers[-1] * DENOMINATOR_BASES[index])
    return powers[exponent]


class RationalFunction:
    """Element of ``QQ(b,c,d,e)`` with a factored positive denominator."""

    def __init__(self, numerator=0, denominator=None):
        self.numerator = (numerator if isinstance(numerator, Polynomial4)
                          else Polynomial4.constant(numerator))
        self.denominator = tuple(denominator or (0,) * 6)

    @classmethod
    def polynomial(cls, value):
        return cls(value)

    def __bool__(self):
        return bool(self.numerator)

    def __neg__(self):
        return RationalFunction(-self.numerator, self.denominator)

    def __add__(self, other):
        if not isinstance(other, RationalFunction):
            other = RationalFunction(other)
        common = tuple(max(a, b) for a, b in zip(
            self.denominator, other.denominator))
        left = self.numerator
        right = other.numerator
        for index, (target, a, b) in enumerate(zip(
                common, self.denominator, other.denominator)):
            if target > a:
                left *= _denominator_power(index, target - a)
            if target > b:
                right *= _denominator_power(index, target - b)
        return RationalFunction(left + right, common)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return RationalFunction(other) - self

    def __mul__(self, other):
        if isinstance(other, Polynomial3):
            return NotImplemented
        if not isinstance(other, RationalFunction):
            other = RationalFunction(other)
        return RationalFunction(
            self.numerator * other.numerator,
            tuple(a + b for a, b in zip(
                self.denominator, other.denominator)),
        )

    __rmul__ = __mul__

    def __pow__(self, exponent):
        return RationalFunction(
            self.numerator ** exponent,
            tuple(exponent * value for value in self.denominator),
        )


class Polynomial3:
    """Sparse polynomial in tetrahedral coordinates with rational functions."""

    def __init__(self, coefficients=None):
        self.coefficients = {
            index: value
            for index, value in (coefficients or {}).items()
            if value
        }

    @classmethod
    def constant(cls, value):
        if not isinstance(value, RationalFunction):
            value = RationalFunction(value)
        return cls({(0, 0, 0): value})

    def __add__(self, other):
        if not isinstance(other, Polynomial3):
            other = Polynomial3.constant(other)
        coefficients = dict(self.coefficients)
        for index, value in other.coefficients.items():
            updated = coefficients.get(index, RationalFunction()) + value
            if updated:
                coefficients[index] = updated
            else:
                coefficients.pop(index, None)
        return Polynomial3(coefficients)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial3({index: -value
                            for index, value in self.coefficients.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return Polynomial3.constant(other) - self

    def __mul__(self, other):
        if not isinstance(other, Polynomial3):
            other = Polynomial3.constant(other)
        coefficients = {}
        for (i, j, k), left in self.coefficients.items():
            for (p, q, r), right in other.coefficients.items():
                index = (i + p, j + q, k + r)
                updated = (coefficients.get(index, RationalFunction())
                           + left * right)
                if updated:
                    coefficients[index] = updated
                else:
                    coefficients.pop(index, None)
        return Polynomial3(coefficients)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        result = Polynomial3.constant(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent //= 2
        return result


U = Polynomial3({(1, 0, 0): RationalFunction(1)})
V = Polynomial3({(0, 1, 0): RationalFunction(1)})
W = Polynomial3({(0, 0, 1): RationalFunction(1)})


def _vertex(name):
    zero = RationalFunction(0)
    one = RationalFunction(1)
    if name == "r0":
        return zero, zero, zero
    if name == "r1":
        value = RationalFunction(E4, (1, 0, 0, 0, 0, 0))
        return zero, zero, value
    if name == "r2":
        value = RationalFunction(E4, (0, 1, 0, 0, 0, 0))
        return zero, value, value
    if name == "r3":
        value = RationalFunction(E4, (0, 0, 1, 0, 0, 0))
        return value, value, value
    if name == "p0":
        return zero, zero, one
    if name == "q2":
        value = RationalFunction(D4 + E4, (0, 0, 0, 1, 0, 0))
        return zero, value, one
    if name == "q3":
        value = RationalFunction(D4 + E4, (0, 0, 0, 0, 1, 0))
        return value, value, one
    if name == "pc":
        return zero, one, one
    if name == "q4":
        value = RationalFunction(C4 + D4 + E4,
                                 (0, 0, 0, 0, 0, 1))
        return value, one, one
    raise ValueError(f"unknown vertex {name}")


def _compose(polynomial, xyz, bcde, vertices):
    affine = []
    for coordinate in range(3):
        affine.append(
            Polynomial3.constant(vertices[0][coordinate])
            + U * (vertices[1][coordinate] - vertices[0][coordinate])
            + V * (vertices[2][coordinate] - vertices[0][coordinate])
            + W * (vertices[3][coordinate] - vertices[0][coordinate])
        )
    output = Polynomial3()
    for powers, coefficient in sp.Poly(polynomial, *xyz).terms():
        coefficient = RationalFunction.polynomial(
            Polynomial4.from_sympy(coefficient, bcde))
        output += (affine[0] ** powers[0]
                   * affine[1] ** powers[1]
                   * affine[2] ** powers[2]
                   * coefficient)
    return output


def _bernstein_coefficients(polynomial, degree):
    output = {}
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            for k in range(degree + 1 - i - j):
                ell = degree - i - j - k
                value = RationalFunction()
                for (p, q, r), coefficient in polynomial.coefficients.items():
                    if i >= p and j >= q and k >= r:
                        denominator = (
                            math.comb(degree, p + q + r)
                            * math.factorial(p + q + r)
                            // (math.factorial(p)
                                * math.factorial(q)
                                * math.factorial(r))
                        )
                        multiplier = Fraction(
                            math.comb(i, p) * math.comb(j, q)
                            * math.comb(k, r),
                            denominator,
                        )
                        value += coefficient * multiplier
                output[(i, j, k, ell)] = value
    return output


def _symbolic_regions():
    x, y, z, b, c, d, e = sp.symbols("x y z b c d e")
    s = b * x + c * y + d * z + e
    common = (1 - x) * (1 - y) * (1 - z)

    alpha_a = (1 - e) ** 4
    polynomial_a = alpha_a * common - (1 - s) ** 4

    p1 = d + e
    alpha_b = (1 - p1) ** 3 * (1 + 3 * p1)
    denominator_b = common * z * (z - x) * (z - y)
    numerator_b = (
        alpha_b * denominator_b
        - (1 - s) ** 4 * z * (z - x) * (z - y)
        + (z - s) ** 4 * (1 - x) * (1 - y)
    )
    quotient_b, remainder_b = sp.div(
        sp.Poly(numerator_b, z, y, x),
        sp.Poly(1 - z, z, y, x),
    )
    if not remainder_b.is_zero:
        raise AssertionError("region-B factorization failed")

    p2 = c + d + e
    alpha_c = (1 - p2) ** 2 * (1 + 2 * p2 + 3 * p2 ** 2)
    denominator_c = (common * z * (z - x) * (z - y)
                     * y * (y - x))
    numerator_c = (
        alpha_c * denominator_c
        - (1 - s) ** 4 * z * (z - x) * (z - y) * y * (y - x)
        + (z - s) ** 4 * (1 - x) * (1 - y) * y * (y - x)
        - (y - s) ** 4 * (1 - x) * (1 - z) * z * (z - x)
    )
    quotient_c = sp.Poly(numerator_c, z, y, x)
    for factor in (1 - z, 1 - y, z - y):
        quotient_c, remainder = sp.div(
            quotient_c, sp.Poly(factor, z, y, x))
        if not remainder.is_zero:
            raise AssertionError("region-C factorization failed")

    return (
        (x, y, z),
        (b, c, d, e),
        {
            "A": {
                "degree": 4,
                "alpha_identity": "alpha=(1-e)^4",
                "extracted_nonnegative_factors": [],
                "polynomial": polynomial_a,
            },
            "B": {
                "degree": 6,
                "alpha_identity":
                    "alpha=(1-d-e)^3*(1+3*(d+e))",
                "extracted_nonnegative_factors": ["1-z"],
                "polynomial": quotient_b.as_expr(),
            },
            "C": {
                "degree": 6,
                "alpha_identity": (
                    "alpha=(1-c-d-e)^2*"
                    "(1+2*(c+d+e)+3*(c+d+e)^2)"
                ),
                "extracted_nonnegative_factors":
                    ["1-z", "1-y", "z-y"],
                "polynomial": quotient_c.as_expr(),
            },
        },
    )


def derive():
    xyz, bcde, regions = _symbolic_regions()
    output = {
        "schema_version": 1,
        "claim": (
            "Exact polynomial and tetrahedral Bernstein structure for the "
            "n=4 Gaffke-domination certificate."
        ),
        "variables": [str(value) for value in bcde],
        "vertex_formulas": VERTEX_FORMULAS,
        "regions": {},
    }

    for region_name, region in regions.items():
        polynomial = region["polynomial"]
        polynomial_coefficients = [
            {
                "powers": list(powers),
                "expression": str(sp.factor(coefficient)),
            }
            for powers, coefficient in sp.Poly(polynomial, *xyz).terms()
        ]
        tetrahedra = []
        for vertex_names in REGION_TETRAHEDRA[region_name]:
            vertices = [_vertex(name) for name in vertex_names]
            transformed = _compose(polynomial, xyz, bcde, vertices)
            coefficients = _bernstein_coefficients(
                transformed, region["degree"])
            zero_indices = [
                list(index) for index, value in coefficients.items()
                if not value
            ]
            tetrahedra.append({
                "vertices": vertex_names,
                "structural_zero_indices": zero_indices,
                "coefficient_count": len(coefficients),
                "nonzero_coefficient_count": (
                    len(coefficients) - len(zero_indices)
                ),
            })
            print(
                f"region {region_name}, tetrahedron {len(tetrahedra)-1}: "
                f"{len(zero_indices)} structural zeros"
            )
        output["regions"][region_name] = {
            "degree": region["degree"],
            "alpha_identity": region["alpha_identity"],
            "extracted_nonnegative_factors":
                region["extracted_nonnegative_factors"],
            "polynomial_coefficients": polynomial_coefficients,
            "tetrahedra": tetrahedra,
        }
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    args.out.write_text(json.dumps(derive(), indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
