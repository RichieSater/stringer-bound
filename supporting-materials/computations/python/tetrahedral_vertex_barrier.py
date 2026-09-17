"""Exact Bernstein-multiplier certificate for the tetrahedral barrier.

Let ``D`` be uniform on the four-coordinate simplex and consider a regular
two-versus-two section with centered knots

    (-1, -u, t*v, t),       0 < u,v < 1, t > 0.

Write ``p`` for the upper-cap volume, ``gamma`` for the section centroid,
and ``M = gamma[0] + gamma[1]``.  Orderedness of the two upper centroid
coordinates implies ``R(u,v,t) <= 0``.  This module verifies exactly that

    I_M(2,2) - p = t**2 * P / Delta,

where ``Delta`` is positive, and checks a rational Bernstein multiplier
``S`` for which

    H = P + R*S

has only nonnegative tensor Bernstein coefficients on the unit cube.  The
multiplier does too.  Hence ``P = H - R*S >= 0`` whenever ``R <= 0``.

The floating-point linear program that proposed ``S`` is not part of the
proof.  Every identity and sign used here is recomputed with SymPy integers
and Python ``Fraction`` arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import comb
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "tetrahedral-vertex-barrier-certificate.json"
)

VARIABLES = ("u", "v", "t")
R_DEGREE = (3, 2, 4)
P_DEGREE = (9, 9, 11)
MULTIPLIER_DEGREE = (6, 7, 7)
RESIDUAL_DEGREE = (9, 9, 11)
BOUNDARY_ORDER_DEGREE = (2, 0, 3)
BOUNDARY_MULTIPLIER_DEGREE = (2, 0, 2)
BOUNDARY_RESIDUAL_DEGREE = (4, 0, 5)

# Nonzero tensor-Bernstein coefficients of S, in variable order (u,v,t).
# This is the complete proof witness; omitted coefficients are zero.
MULTIPLIER_COEFFICIENTS = (
    ((0, 7, 7), Fraction(16)),
    ((3, 1, 7), Fraction(561, 10640)),
    ((3, 2, 7), Fraction(1749, 5320)),
    ((4, 1, 7), Fraction(11, 35)),
    ((4, 2, 7), Fraction(7796, 1995)),
    ((4, 3, 7), Fraction(4078, 525)),
    ((5, 0, 4), Fraction(2173, 2800)),
    ((5, 0, 5), Fraction(1873, 720)),
    ((5, 0, 6), Fraction(289, 56)),
    ((5, 0, 7), Fraction(25, 4)),
    ((5, 1, 3), Fraction(3460739, 537600)),
    ((5, 1, 4), Fraction(167807, 19600)),
    ((5, 1, 5), Fraction(85367, 8820)),
    ((5, 1, 6), Fraction(641, 49)),
    ((5, 1, 7), Fraction(122, 7)),
    ((5, 2, 1), Fraction(132193, 125440)),
    ((5, 2, 2), Fraction(738323, 161280)),
    ((5, 2, 3), Fraction(248303, 39200)),
    ((5, 2, 4), Fraction(316997, 29400)),
    ((5, 2, 5), Fraction(168013, 8820)),
    ((5, 2, 6), Fraction(1420, 49)),
    ((5, 2, 7), Fraction(299, 7)),
    ((5, 3, 1), Fraction(212979, 78400)),
    ((5, 3, 2), Fraction(90327, 19600)),
    ((5, 3, 3), Fraction(132423, 14000)),
    ((5, 3, 4), Fraction(214063, 12250)),
    ((5, 3, 5), Fraction(113236, 3675)),
    ((5, 3, 6), Fraction(13002, 245)),
    ((5, 3, 7), Fraction(618, 7)),
    ((5, 4, 0), Fraction(9, 10)),
    ((5, 4, 1), Fraction(72, 35)),
    ((5, 4, 2), Fraction(9, 2)),
    ((5, 4, 3), Fraction(333, 35)),
    ((5, 4, 4), Fraction(684, 35)),
    ((5, 4, 5), Fraction(1368, 35)),
    ((5, 4, 6), Fraction(2664, 35)),
    ((5, 4, 7), Fraction(144)),
    ((5, 5, 0), Fraction(3, 7)),
    ((5, 5, 1), Fraction(15, 14)),
    ((5, 5, 2), Fraction(18, 7)),
    ((5, 5, 3), Fraction(6)),
    ((5, 5, 4), Fraction(96, 7)),
    ((5, 5, 5), Fraction(216, 7)),
    ((5, 5, 6), Fraction(480, 7)),
    ((5, 5, 7), Fraction(1056, 7)),
    ((6, 0, 7), Fraction(6)),
    ((6, 1, 4), Fraction(23523, 2450)),
    ((6, 1, 5), Fraction(21464, 735)),
    ((6, 1, 6), Fraction(42)),
    ((6, 1, 7), Fraction(396, 7)),
    ((6, 2, 2), Fraction(882311, 70560)),
    ((6, 2, 3), Fraction(104311, 3675)),
    ((6, 2, 4), Fraction(169606, 3675)),
    ((6, 2, 5), Fraction(134912, 2205)),
    ((6, 2, 6), Fraction(4392, 49)),
    ((6, 2, 7), Fraction(888, 7)),
    ((6, 3, 1), Fraction(6327, 560)),
    ((6, 3, 2), Fraction(77806, 3675)),
    ((6, 3, 3), Fraction(214691, 6125)),
    ((6, 3, 4), Fraction(353344, 6125)),
    ((6, 3, 5), Fraction(348836, 3675)),
    ((6, 3, 6), Fraction(36864, 245)),
    ((6, 3, 7), Fraction(8208, 35)),
    ((6, 4, 0), Fraction(27, 5)),
    ((6, 4, 1), Fraction(351, 35)),
    ((6, 4, 2), Fraction(648, 35)),
    ((6, 4, 3), Fraction(1188, 35)),
    ((6, 4, 4), Fraction(432, 7)),
    ((6, 4, 5), Fraction(3888, 35)),
    ((6, 4, 6), Fraction(6912, 35)),
    ((6, 4, 7), Fraction(1728, 5)),
    ((6, 5, 0), Fraction(18, 7)),
    ((6, 5, 1), Fraction(36, 7)),
    ((6, 5, 2), Fraction(72, 7)),
    ((6, 5, 3), Fraction(144, 7)),
    ((6, 5, 4), Fraction(288, 7)),
    ((6, 5, 5), Fraction(576, 7)),
    ((6, 5, 6), Fraction(1152, 7)),
    ((6, 5, 7), Fraction(2304, 7)),
)


def _fraction(value: sp.Rational) -> Fraction:
    return Fraction(int(value.p), int(value.q))


def derive_polynomials():
    """Derive the cap, centroid-order factor, and comparison numerator."""
    u, v, t = sp.symbols("u v t", positive=True)
    A, B, C, D = sp.symbols("A B C D", positive=True)

    # Divided-difference upper-tail formula for knots (-A,-B,C,D).
    cap = sp.factor(
        C**3 / ((C + A) * (C + B) * (C - D))
        + D**3 / ((D + A) * (D + B) * (D - C))
    )
    density = sp.factor(
        3 * (
            C**2 / ((C + A) * (C + B) * (C - D))
            + D**2 / ((D + A) * (D + B) * (D - C))
        )
    )
    gradient = (
        -sp.diff(cap, A),
        -sp.diff(cap, B),
        sp.diff(cap, C),
        sp.diff(cap, D),
    )
    substitution = {A: 1, B: u, C: t * v, D: t}
    p = sp.factor(cap.subs(substitution))
    gamma = tuple(sp.factor(item.subs(substitution) / density.subs(substitution))
                  for item in gradient)
    if sp.factor(sum(gamma) - 1) != 0:
        raise AssertionError("section-centroid coordinates do not sum to one")

    order_difference = sp.factor(gamma[3] - gamma[2])
    order_numerator, order_denominator = map(
        sp.factor, sp.fraction(order_difference))
    R = sp.factor(order_numerator / (v - 1))
    expected_R = (
        t**4 * v**2 * (u**2 + u + 1)
        + 2 * t**3 * u * v * (u + 1) * (v + 1)
        + t**2 * u**2 * (v**2 + 4 * v + 1)
        - u**3
    )
    if sp.factor(R - expected_R) != 0:
        raise AssertionError("upper-centroid order factor changed")

    M = sp.factor(gamma[0] + gamma[1])
    comparison = sp.factor(3 * M**2 - 2 * M**3 - p)
    numerator, denominator = map(sp.factor, sp.fraction(comparison))
    if sp.rem(sp.Poly(numerator, t), sp.Poly(t**2, t)) != 0:
        raise AssertionError("comparison numerator lost its t^2 factor")
    P = sp.Poly(sp.cancel(numerator / t**2), u, v, t)
    expected_denominator = (
        27 * (t + 1)**3 * (t + u)**3
        * (t * v + 1)**3 * (t * v + u)**3
        * (t * u * v + t * v + u * v + u)**3
    )
    if sp.factor(denominator - expected_denominator) != 0:
        raise AssertionError("comparison denominator changed")
    if P.degree_list() != P_DEGREE:
        raise AssertionError("unexpected comparison tensor degree")
    return (u, v, t), p, gamma, sp.Poly(R, u, v, t), P, denominator


def power_to_bernstein(poly: sp.Poly, degree: tuple[int, int, int]):
    """Return exact tensor-Bernstein coefficients on ``[0,1]^3``."""
    result: dict[tuple[int, int, int], Fraction] = {}
    for powers, coefficient in poly.terms():
        coefficient = _fraction(coefficient)
        if any(powers[index] > degree[index] for index in range(3)):
            raise ValueError("requested Bernstein degree is too small")
        for i in range(powers[0], degree[0] + 1):
            factor_i = Fraction(comb(i, powers[0]), comb(degree[0], powers[0]))
            for j in range(powers[1], degree[1] + 1):
                factor_ij = factor_i * Fraction(
                    comb(j, powers[1]), comb(degree[1], powers[1]))
                for k in range(powers[2], degree[2] + 1):
                    index = (i, j, k)
                    result[index] = result.get(index, Fraction(0)) + (
                        coefficient * factor_ij
                        * Fraction(comb(k, powers[2]),
                                   comb(degree[2], powers[2]))
                    )
    return result


def multiply_bernstein(
    left: dict[tuple[int, int, int], Fraction],
    left_degree: tuple[int, int, int],
    right: dict[tuple[int, int, int], Fraction],
    right_degree: tuple[int, int, int],
):
    """Multiply two exact tensor-Bernstein expansions."""
    total_degree = tuple(left_degree[q] + right_degree[q] for q in range(3))
    result: dict[tuple[int, int, int], Fraction] = {}
    for left_index, left_value in left.items():
        if left_value == 0:
            continue
        for right_index, right_value in right.items():
            if right_value == 0:
                continue
            index = tuple(left_index[q] + right_index[q] for q in range(3))
            factor = Fraction(1)
            for q in range(3):
                factor *= Fraction(
                    comb(left_degree[q], left_index[q])
                    * comb(right_degree[q], right_index[q]),
                    comb(total_degree[q], index[q]),
                )
            result[index] = result.get(index, Fraction(0)) + (
                left_value * right_value * factor)
    return result


def _all_indices(degree: tuple[int, int, int]):
    for i in range(degree[0] + 1):
        for j in range(degree[1] + 1):
            for k in range(degree[2] + 1):
                yield i, j, k


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _digest(coefficients, degree):
    digest = hashlib.sha256()
    for index in _all_indices(degree):
        value = coefficients.get(index, Fraction(0))
        digest.update(("%d,%d,%d:%d/%d\n" % (
            *index, value.numerator, value.denominator)).encode())
    return digest.hexdigest()


def verify_n3_application():
    """Check the exact algebra locating the monotone Stringer-weight range."""
    x, r, z = sp.symbols("x r z")

    def f(value):
        return 1 - 3 * value**2 + 2 * value**3

    last_weight_identity = sp.factor(
        (1 - x)**3 - f(2 * x) - x * (-3 + 15 * x - 17 * x**2)
    )
    if last_weight_identity != 0:
        raise AssertionError("last Stringer-weight identity changed")

    middle_identity = sp.factor(
        4 * (f((1 - r + z) / 2) - r**3)
        - 3 * (z - r) * (z**2 + r**2 - 1)
    )
    if sp.expand(middle_identity + 2 * (r**3 + z**3 - 1)) != 0:
        raise AssertionError("middle Stringer-weight identity changed")

    first_weight_identity = sp.factor(
        1 - z**3 - f(2 * z - 1)
        - (1 - z) * (17 * z**2 - 19 * z + 5)
    )
    if first_weight_identity != 0:
        raise AssertionError("first Stringer-weight identity changed")

    x_star = (sp.Integer(15) - sp.sqrt(21)) / 34
    alpha_star = sp.factor((1 - x_star)**3)
    if sp.simplify(-3 + 15 * x_star - 17 * x_star**2) != 0:
        raise AssertionError("monotone-weight endpoint is not an exact root")
    if not (sp.ask(sp.Q.positive(x_star))
            and sp.ask(sp.Q.positive(sp.Rational(1, 2) - x_star))
            and sp.ask(sp.Q.positive(sp.Rational(1, 2) - alpha_star))):
        raise AssertionError("monotone-weight endpoint left its stated range")
    return {
        "x_star": "(15-sqrt(21))/34",
        "alpha_star": "((19+sqrt(21))/34)^3",
        "last_weight_identity": (
            "(1-x)^3-f(2*x)=x*(-3+15*x-17*x^2)"
        ),
        "middle_weight_identity": (
            "4*(f((1-r+z)/2)-r^3)=3*(z-r)*(z^2+r^2-1) "
            "when r^3+z^3=1"
        ),
        "first_weight_identity": (
            "1-z^3-f(2*z-1)=(1-z)*(17*z^2-19*z+5)"
        ),
        "pointwise_comparison_sharpness": (
            "For alpha>alpha_star, c_2>c_3 and the cap along knots "
            "(0,0,h,1) has derivative g_alpha*(c_3-c_2)/2<0 at h=1; "
            "therefore its value exceeds alpha for h<1 sufficiently close "
            "to one."
        ),
    }


def verify_certificate():
    """Recompute every identity and exact coefficient sign in the witness."""
    variables, p, gamma, R, P, denominator = derive_polynomials()
    n3_application = verify_n3_application()
    u, v, t = variables
    r_bernstein = power_to_bernstein(R, R_DEGREE)
    p_bernstein = power_to_bernstein(P, RESIDUAL_DEGREE)
    multiplier = dict(MULTIPLIER_COEFFICIENTS)
    if len(multiplier) != len(MULTIPLIER_COEFFICIENTS):
        raise AssertionError("duplicate multiplier index")
    if any(value <= 0 for value in multiplier.values()):
        raise AssertionError("multiplier coefficient is not positive")

    product = multiply_bernstein(
        r_bernstein, R_DEGREE, multiplier, MULTIPLIER_DEGREE)
    residual = {
        index: p_bernstein.get(index, Fraction(0))
        + product.get(index, Fraction(0))
        for index in _all_indices(RESIDUAL_DEGREE)
    }
    negative = {index: value for index, value in residual.items() if value < 0}
    if negative:
        first = min(negative.items(), key=lambda item: item[1])
        raise AssertionError(f"negative residual Bernstein coefficient: {first}")

    residual_positive = [value for value in residual.values() if value > 0]
    multiplier_positive = list(multiplier.values())
    if not residual_positive:
        raise AssertionError("the residual witness is identically zero")
    # When v=1 the two upper knots coincide, so their centroid coordinates
    # are equal and the factor (v-1)R supplies no condition.  The remaining
    # orderedness inequality gamma_1 <= gamma_2 closes that boundary with a
    # particularly small second exact multiplier.
    boundary_M = sp.factor((gamma[0] + gamma[1]).subs(v, 1))
    boundary_comparison = sp.factor(
        3 * boundary_M**2 - 2 * boundary_M**3 - p.subs(v, 1))
    boundary_numerator, boundary_denominator = map(
        sp.factor, sp.fraction(boundary_comparison))
    boundary_P = sp.Poly(
        sp.cancel(boundary_numerator / (t**2 * (1 - u)**2)), u, v, t)
    boundary_order, boundary_order_denominator = sp.fraction(
        sp.cancel((gamma[2] - gamma[1]).subs(v, 1)))
    boundary_order = sp.Poly(boundary_order, u, v, t)
    if boundary_P.degree_list() != (4, 0, 5):
        raise AssertionError("unexpected repeated-upper comparison degree")
    if boundary_order.degree_list() != BOUNDARY_ORDER_DEGREE:
        raise AssertionError("unexpected repeated-upper order degree")
    if boundary_denominator != (
        27 * (t + 1)**3 * (t + u)**3 * (t * u + t + 2 * u)**3
    ):
        raise AssertionError("repeated-upper comparison denominator changed")
    expected_boundary_order_denominator = (
        3 * (t + 1) * (t + u) * (t * u + t + 2 * u)
    )
    if sp.factor(
        boundary_order_denominator - expected_boundary_order_denominator
    ) != 0:
        raise AssertionError("repeated-upper order denominator changed")

    boundary_p_bernstein = power_to_bernstein(
        boundary_P, BOUNDARY_RESIDUAL_DEGREE)
    boundary_order_bernstein = power_to_bernstein(
        boundary_order, BOUNDARY_ORDER_DEGREE)
    # 2 B_{0,2}(u) B_{0,0}(v) B_{2,2}(t)
    boundary_multiplier = {(0, 0, 2): Fraction(2)}
    boundary_product = multiply_bernstein(
        boundary_order_bernstein,
        BOUNDARY_ORDER_DEGREE,
        boundary_multiplier,
        BOUNDARY_MULTIPLIER_DEGREE,
    )
    boundary_residual = {
        index: boundary_p_bernstein.get(index, Fraction(0))
        - boundary_product.get(index, Fraction(0))
        for index in _all_indices(BOUNDARY_RESIDUAL_DEGREE)
    }
    boundary_negative = {
        index: value for index, value in boundary_residual.items() if value < 0
    }
    if boundary_negative:
        first = min(boundary_negative.items(), key=lambda item: item[1])
        raise AssertionError(
            f"negative repeated-upper Bernstein coefficient: {first}")
    boundary_positive = [
        value for value in boundary_residual.values() if value > 0]

    return {
        "denominator": str(denominator),
        "multiplier": multiplier,
        "residual": residual,
        "multiplier_positive": multiplier_positive,
        "residual_positive": residual_positive,
        "boundary_multiplier": boundary_multiplier,
        "boundary_residual": boundary_residual,
        "boundary_positive": boundary_positive,
        "n3_application": n3_application,
    }


def build_certificate() -> dict[str, object]:
    verified = verify_certificate()
    multiplier = verified["multiplier"]
    residual = verified["residual"]
    multiplier_positive = verified["multiplier_positive"]
    residual_positive = verified["residual_positive"]
    boundary_multiplier = verified["boundary_multiplier"]
    boundary_residual = verified["boundary_residual"]
    boundary_positive = verified["boundary_positive"]
    n3_application = verified["n3_application"]
    total_multiplier = (MULTIPLIER_DEGREE[0] + 1) * (
        MULTIPLIER_DEGREE[1] + 1) * (MULTIPLIER_DEGREE[2] + 1)
    total_residual = (RESIDUAL_DEGREE[0] + 1) * (
        RESIDUAL_DEGREE[1] + 1) * (RESIDUAL_DEGREE[2] + 1)
    return {
        "schema_version": 1,
        "status": (
            "Exact rational certificate for the ordinary four-coordinate "
            "two-versus-two section-centroid value barrier."
        ),
        "theorem": {
            "normalized_knots": ["-1", "-u", "t*v", "t"],
            "domain": "0<u<1, 0<v<1, t>0",
            "hypothesis": (
                "The section centroid is ordered; its upper endpoint "
                "inequality implies R(u,v,t)<=0 and t<1."
            ),
            "conclusion": (
                "For p=Pr{x dot D>0} and M=E[D0+D1|x dot D=0], "
                "p<=I_M(2,2)=3*M^2-2*M^3."
            ),
            "global_consequence": (
                "Together with the elementary one-versus-three strata, "
                "this proves vertex maximization for every nondecreasing "
                "four-coordinate barycentric weight vector."
            ),
        },
        "identity": {
            "comparison": "I_M(2,2)-p=t^2*P/Delta with Delta>0",
            "order_factor": (
                "R=t^4*v^2*(u^2+u+1)+2*t^3*u*v*(u+1)*(v+1)"
                "+t^2*u^2*(v^2+4*v+1)-u^3"
            ),
            "multiplier": "H=P+R*S, equivalently P=H-R*S",
            "logic": "R<=0 and H,S>=0 imply P>=0",
            "repeated_upper_boundary": (
                "At v=1, gamma_2-gamma_1 has numerator G; "
                "P_boundary=H_boundary+G*2*(1-u)^2*t^2, with both "
                "G and H_boundary nonnegative under centroid ordering."
            ),
        },
        "bernstein_basis": {
            "variables": list(VARIABLES),
            "order_factor_degree": list(R_DEGREE),
            "comparison_degree": list(P_DEGREE),
            "multiplier_degree": list(MULTIPLIER_DEGREE),
            "residual_degree": list(RESIDUAL_DEGREE),
        },
        "multiplier": {
            "total_coefficients": total_multiplier,
            "positive_coefficients": len(multiplier_positive),
            "zero_coefficients": total_multiplier - len(multiplier_positive),
            "negative_coefficients": 0,
            "smallest_positive": _fraction_record(min(multiplier_positive)),
            "largest": _fraction_record(max(multiplier_positive)),
            "sha256_all_coefficients": _digest(multiplier, MULTIPLIER_DEGREE),
            "nonzero_coefficients": [
                {"index": list(index), **_fraction_record(value)}
                for index, value in MULTIPLIER_COEFFICIENTS
            ],
        },
        "residual": {
            "total_coefficients": total_residual,
            "positive_coefficients": len(residual_positive),
            "zero_coefficients": total_residual - len(residual_positive),
            "negative_coefficients": 0,
            "smallest_positive": _fraction_record(min(residual_positive)),
            "largest": _fraction_record(max(residual_positive)),
            "sha256_all_coefficients": _digest(residual, RESIDUAL_DEGREE),
        },
        "repeated_upper_boundary": {
            "variables": ["u", "v (degree zero)", "t"],
            "order_numerator_degree": list(BOUNDARY_ORDER_DEGREE),
            "multiplier_degree": list(BOUNDARY_MULTIPLIER_DEGREE),
            "multiplier_nonzero_coefficients": [
                {"index": list(index), **_fraction_record(value)}
                for index, value in boundary_multiplier.items()
            ],
            "residual_degree": list(BOUNDARY_RESIDUAL_DEGREE),
            "residual_positive_coefficients": len(boundary_positive),
            "residual_zero_coefficients": (
                (BOUNDARY_RESIDUAL_DEGREE[0] + 1)
                * (BOUNDARY_RESIDUAL_DEGREE[1] + 1)
                * (BOUNDARY_RESIDUAL_DEGREE[2] + 1)
                - len(boundary_positive)
            ),
            "residual_negative_coefficients": 0,
            "smallest_positive": _fraction_record(min(boundary_positive)),
            "largest": _fraction_record(max(boundary_positive)),
            "sha256_all_residual_coefficients": _digest(
                boundary_residual, BOUNDARY_RESIDUAL_DEGREE),
        },
        "n3_stringer_application": n3_application,
        "arithmetic": {
            "symbolic_identities": "SymPy exact integer polynomial arithmetic",
            "coefficient_signs": "Python Fraction arithmetic",
            "floating_point_role": "none in verification or certificate output",
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        "tetrahedral barrier certified: multiplier %d positive/%d total; "
        "residual %d positive/%d total"
        % (
            certificate["multiplier"]["positive_coefficients"],
            certificate["multiplier"]["total_coefficients"],
            certificate["residual"]["positive_coefficients"],
            certificate["residual"]["total_coefficients"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
