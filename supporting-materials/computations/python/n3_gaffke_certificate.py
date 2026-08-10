"""Exact certificate for the ``n=3`` conventional-level theorem.

For ``alpha`` equal to 0.01, 0.05, or 0.10, this module certifies that
the binomial-factor Stringer bound pointwise dominates the one-sided
Gaffke/Learned--Miller--Thomas bound.  The latter has finite-sample
distribution-free coverage, so the pointwise comparison proves Stringer
coverage at 99%, 95%, and 90% confidence when ``n=3``.

The proof reduces the comparison to the volume of a halfspace cap of the
uniform three-simplex.  On each of the three possible knot regions, the
needed inequality is either an AM--GM consequence or the nonnegativity of a
low-degree polynomial.  Polynomial nonnegativity is certified by expressing
the polynomial in a Bernstein basis over one or two triangles and checking
every coefficient.

All arithmetic in this file is rational interval arithmetic.  The
Clopper--Pearson factors are enclosed by dyadic intervals whose endpoint
signs are checked using integer arithmetic in :mod:`stringer`.  The static
formula file was derived symbolically from the displayed cap-volume formula;
it contains the 42 degree-five Bernstein coefficients used for the middle
region.  No floating-point result is used to decide a sign.

Usage::

    python n3_gaffke_certificate.py
    python n3_gaffke_certificate.py --out n3-gaffke-certificate.json
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple

from stringer import exact_binomial_factor_brackets


ALPHAS = ("0.01", "0.05", "0.10")
FACTOR_BITS = 120
HERE = Path(__file__).resolve().parent
CERTIFICATE_DIR = HERE.parent / "certificates"
FORMULA_PATH = CERTIFICATE_DIR / "n3-gaffke-bernstein-formulas.json"


@dataclass(frozen=True)
class Interval:
    """Closed interval with exact rational endpoints."""

    lower: Fraction
    upper: Fraction

    def __init__(self, lower=0, upper=None):
        if isinstance(lower, Interval):
            object.__setattr__(self, "lower", lower.lower)
            object.__setattr__(self, "upper", lower.upper)
            return
        lo = Fraction(lower)
        hi = lo if upper is None else Fraction(upper)
        if lo > hi:
            raise ValueError("interval endpoints are reversed")
        object.__setattr__(self, "lower", lo)
        object.__setattr__(self, "upper", hi)

    def __add__(self, other):
        if isinstance(other, Polynomial2):
            return NotImplemented
        other = Interval(other)
        return Interval(self.lower + other.lower,
                        self.upper + other.upper)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        if isinstance(other, Polynomial2):
            return NotImplemented
        return self + (-Interval(other))

    def __rsub__(self, other):
        return Interval(other) - self

    def __mul__(self, other):
        if isinstance(other, Polynomial2):
            return NotImplemented
        other = Interval(other)
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lower <= 0 <= self.upper:
            raise ZeroDivisionError("interval contains zero")
        return Interval(1 / self.upper, 1 / self.lower)

    def __truediv__(self, other):
        return self * Interval(other).reciprocal()

    def __rtruediv__(self, other):
        return Interval(other) / self

    def __pow__(self, exponent: int):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("only nonnegative integer powers are supported")
        result = Interval(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result


class Polynomial2:
    """Polynomial in ``U,W`` with rational-interval coefficients."""

    def __init__(self, coefficients=None):
        self.coefficients: Dict[Tuple[int, int], Interval] = {}
        for index, value in (coefficients or {}).items():
            value = Interval(value)
            if value.lower != 0 or value.upper != 0:
                self.coefficients[index] = value

    @classmethod
    def constant(cls, value):
        return cls({(0, 0): Interval(value)})

    def __add__(self, other):
        if not isinstance(other, Polynomial2):
            other = Polynomial2.constant(other)
        coefficients = dict(self.coefficients)
        for index, value in other.coefficients.items():
            coefficients[index] = coefficients.get(
                index, Interval(0)) + value
        return Polynomial2(coefficients)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial2({index: -value
                            for index, value in self.coefficients.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return Polynomial2.constant(other) - self

    def __mul__(self, other):
        if not isinstance(other, Polynomial2):
            other = Polynomial2.constant(other)
        coefficients: Dict[Tuple[int, int], Interval] = {}
        for (i, j), left in self.coefficients.items():
            for (k, ell), right in other.coefficients.items():
                index = (i + k, j + ell)
                coefficients[index] = coefficients.get(
                    index, Interval(0)) + left * right
        return Polynomial2(coefficients)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * Interval(other).reciprocal()

    def __pow__(self, exponent: int):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("only nonnegative integer powers are supported")
        result = Polynomial2.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result


U = Polynomial2({(1, 0): Interval(1)})
W = Polynomial2({(0, 1): Interval(1)})


def bernstein_coefficients(polynomial: Polynomial2, degree: int):
    """Return triangular degree-``degree`` Bernstein coefficients.

    On the standard triangle ``U >= 0, W >= 0, U + W <= 1``,
    ``U**p W**q`` has coefficient

        C(i,p) C(j,q) / (C(degree,p+q) C(p+q,p))

    at the Bernstein basis function indexed by ``(i,j,degree-i-j)``.
    All multipliers are nonnegative rationals, so interval enclosure is
    preserved.
    """
    output = {}
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            value = Interval(0)
            for (p, q), coefficient in polynomial.coefficients.items():
                if i >= p and j >= q:
                    multiplier = Fraction(
                        math.comb(i, p) * math.comb(j, q),
                        math.comb(degree, p + q) * math.comb(p + q, p),
                    )
                    value += coefficient * multiplier
            output[(i, j, degree - i - j)] = value
    return output


def _evaluate_formula(text: str, values: Mapping[str, Interval]) -> Interval:
    """Evaluate a controlled arithmetic expression by exact intervals."""

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return Interval(node.value)
        if isinstance(node, ast.Name) and node.id in values:
            return values[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            if isinstance(node.op, ast.Pow):
                if not (isinstance(node.right, ast.Constant)
                        and isinstance(node.right.value, int)):
                    raise ValueError("formula exponent is not an integer")
                return left ** node.right.value
            right = visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError(f"unsupported formula syntax: {ast.dump(node)}")

    return visit(ast.parse(text, mode="eval"))


def _fraction_record(value: Fraction):
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": f"{float(value):.12e}",
    }


def _interval_record(value: Interval):
    return {
        "lower": _fraction_record(value.lower),
        "upper": _fraction_record(value.upper),
    }


def _minimum_positive(coefficients: Mapping[Tuple[int, int, int], Interval],
                      excluded: Iterable[Tuple[int, int, int]] = ()):
    excluded = set(excluded)
    candidates = [
        (value.lower, index, value)
        for index, value in coefficients.items()
        if index not in excluded
    ]
    lower, index, interval = min(candidates)
    if lower <= 0:
        raise AssertionError(
            f"Bernstein coefficient {index} is not certified positive: "
            f"[{interval.lower}, {interval.upper}]")
    return index, interval


def _load_formulas():
    raw = FORMULA_PATH.read_bytes()
    data = json.loads(raw)
    if data["degree"] != 5 or len(data["triangles"]) != 2:
        raise ValueError("unexpected n=3 formula schema")
    return data, hashlib.sha256(raw).hexdigest()


def certify_level(alpha_text: str, formula_data):
    alpha = Fraction(alpha_text)
    brackets = exact_binomial_factor_brackets(
        3, alpha_text, FACTOR_BITS)
    p = [Interval(lower, upper) for lower, upper in brackets]

    # Coefficients of knots (t_(1), t_(2), t_(3), 1), in ascending order.
    a = 1 - p[2]
    b = p[2] - p[1]
    c = p[1] - p[0]
    d = p[0]
    A = 1 - b - c  # A = a + d > 0.
    for name, weight in zip("abcd", (a, b, c, d)):
        if weight.lower <= 0:
            raise AssertionError(f"weight {name} is not strictly positive")

    # Region C: weighted AM--GM applies if b,c <= p_2/3.
    region_c_margins = {
        "p2_over_3_minus_b": p[2] / 3 - b,
        "p2_over_3_minus_c": p[2] / 3 - c,
    }
    if any(value.lower <= 0 for value in region_c_margins.values()):
        raise AssertionError("region-C AM--GM condition is not certified")

    # Region A.  In X=1-x, Y=1-y coordinates its domain is the triangle
    # (a/A,a/A), (1,1), (1,(a+b)/(1-c)).  Certify
    # alpha*X*Y - (a+bX+cY)^3 >= 0 in the degree-three Bernstein basis.
    x0 = a / A
    y2 = (a + b) / (1 - c)
    X = x0 + U * (1 - x0) + W * (1 - x0)
    Y = x0 + U * (1 - x0) + W * (y2 - x0)
    L = a + b * X + c * Y
    region_a_polynomial = Interval(alpha) * X * Y - L ** 3
    region_a_coefficients = bernstein_coefficients(
        region_a_polynomial, 3)
    # At (U,W)=(1,0), alpha=(a+b+c)^3 makes the coefficient exactly zero.
    region_a_zero = (3, 0, 0)
    region_a_zero_interval = region_a_coefficients[region_a_zero]
    if not (region_a_zero_interval.lower <= 0
            <= region_a_zero_interval.upper):
        raise AssertionError(
            "region-A structural-zero enclosure does not contain zero")
    region_a_min_index, region_a_min = _minimum_positive(
        region_a_coefficients, (region_a_zero,))

    # Region B.  The formula file gives the degree-five Bernstein
    # coefficients of A^3 Q on the two triangles in (rho,tau), where
    # Q = y^2(1-x)(y-x)^2 * partial V / partial y.
    formula_values = {"b": b, "c": c, "d": d}
    region_b = []
    for triangle_number, formulas in enumerate(formula_data["triangles"]):
        coefficient_intervals = {}
        structural_zeros = []
        for item in formulas:
            index = tuple(item["index"])
            expression = item["expression"]
            value = _evaluate_formula(expression, formula_values)
            coefficient_intervals[index] = value
            if expression == "0":
                structural_zeros.append(index)
                if value.lower != 0 or value.upper != 0:
                    raise AssertionError("structural zero did not evaluate to zero")
            elif value.lower <= 0:
                raise AssertionError(
                    f"region-B triangle {triangle_number} coefficient "
                    f"{index} is not positive: {value}")
        minimum_index, minimum_interval = _minimum_positive(
            coefficient_intervals, structural_zeros)
        region_b.append({
            "triangle": triangle_number,
            "structural_zero_indices": [list(index)
                                        for index in structural_zeros],
            "minimum_positive_index": list(minimum_index),
            "minimum_positive_coefficient": _interval_record(
                minimum_interval),
        })

    # At y=1 put q=1-x.  Region B requires q >= q0=a/(1-b), and
    # q^2(alpha-V) is the following degree-four polynomial.  Map
    # [q0,1] to U in [0,1] and certify its Bernstein coefficients.
    q0 = a / (1 - b)
    q = q0 + (1 - q0) * U
    h = a + b * q
    boundary_polynomial = (Interval(alpha) * q ** 2
                           - h ** 2 * (3 * q - h * (1 + q)))
    boundary_coefficients = bernstein_coefficients(
        boundary_polynomial, 4)
    # At U=1, alpha=(a+b)^2(3-2(a+b)) makes the coefficient exactly zero.
    boundary_zero = (4, 0, 0)
    boundary_zero_interval = boundary_coefficients[boundary_zero]
    if not (boundary_zero_interval.lower <= 0
            <= boundary_zero_interval.upper):
        raise AssertionError(
            "boundary structural-zero enclosure does not contain zero")
    boundary_min_index, boundary_min = _minimum_positive(
        boundary_coefficients, (boundary_zero,))

    return {
        "alpha": alpha_text,
        "nominal_confidence": str(1 - alpha),
        "factor_bits": FACTOR_BITS,
        "factor_brackets": [
            _interval_record(Interval(lower, upper))
            for lower, upper in brackets
        ],
        "weight_brackets": {
            name: _interval_record(value)
            for name, value in zip("abcd", (a, b, c, d))
        },
        "region_c_am_gm_margins": {
            name: _interval_record(value)
            for name, value in region_c_margins.items()
        },
        "region_a": {
        "degree": 3,
        "structural_zero_index": list(region_a_zero),
        "structural_zero_identity": "alpha=(a+b+c)^3",
        "structural_zero_enclosure": _interval_record(
            region_a_zero_interval),
            "minimum_positive_index": list(region_a_min_index),
            "minimum_positive_coefficient": _interval_record(region_a_min),
        },
        "region_b_derivative": {
            "degree": 5,
            "triangles": region_b,
        },
        "region_b_y_equals_1_boundary": {
            "degree": 4,
            "structural_zero_index": list(boundary_zero),
            "structural_zero_identity":
                "alpha=(a+b)^2*(3-2*(a+b))",
            "structural_zero_enclosure": _interval_record(
                boundary_zero_interval),
            "minimum_positive_index": list(boundary_min_index),
            "minimum_positive_coefficient": _interval_record(boundary_min),
        },
        "conclusion": (
            "For n=3 at this alpha, the binomial Stringer bound "
            "pointwise dominates the one-sided Gaffke bound."
        ),
    }


def build_certificate():
    formulas, formula_sha256 = _load_formulas()
    levels = [certify_level(alpha, formulas) for alpha in ALPHAS]
    return {
        "schema_version": 1,
        "claim": (
            "At n=3 and nominal confidence 90%, 95%, or 99%, the "
            "binomial-factor Stringer bound pointwise dominates the "
            "one-sided Gaffke bound and is therefore distribution-free "
            "conservative under independent sampling."
        ),
        "arithmetic": (
            "Exact rational interval arithmetic; Clopper--Pearson factor "
            "brackets have integer-checked dyadic endpoint signs."
        ),
        "formula_file": str(FORMULA_PATH.relative_to(HERE.parents[2])),
        "formula_sha256": formula_sha256,
        "levels": levels,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        help="write the machine-readable certificate")
    args = parser.parse_args(argv)
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    for level in certificate["levels"]:
        rb = level["region_b_derivative"]["triangles"]
        print(
            f"alpha={level['alpha']}: certified; "
            f"region-A min={level['region_a']['minimum_positive_coefficient']['lower']['decimal']}; "
            f"region-B mins={rb[0]['minimum_positive_coefficient']['lower']['decimal']},"
            f"{rb[1]['minimum_positive_coefficient']['lower']['decimal']}; "
            f"boundary min={level['region_b_y_equals_1_boundary']['minimum_positive_coefficient']['lower']['decimal']}"
        )
    print("n3_gaffke_certificate: all exact sign checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
