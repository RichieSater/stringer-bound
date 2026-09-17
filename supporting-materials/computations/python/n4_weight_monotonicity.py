"""Exact algebraic certificate for the sharp monotone-weight range at n=4.

Let ``p_j(alpha)`` be the one-sided Clopper--Pearson upper limits for four
trials and define the ascending-knot Stringer weights

    (1-p_3, p_3-p_2, p_2-p_1, p_1-p_0, p_0).

The terminal edge changes order at ``alpha_4=r_4**4``, where ``r_4`` is the
unique nontrivial root of

    (2*r-1)**3*(7-6*r)=r**4.

This module certifies the only additional algebraic fact needed to prove that
all five weights are nondecreasing for ``0<alpha<=alpha_4``.  The middle
spacing ``2*p_1-p_2-p_0`` can vanish only at a root of an explicit degree-eight
resultant factor.  Exact Sturm counts isolate its two real roots in ``(0,1)``;
neither can occur in the required low-alpha domain.  Binomial-tail symmetry
then handles the other two nonterminal spacings.  All arithmetic is exact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "n4-weight-monotonicity-certificate.json"
)


def _fraction_record(value) -> dict[str, str]:
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _polynomial_digest(poly: sp.Poly) -> str:
    digest = hashlib.sha256()
    digest.update(("variable=" + str(poly.gens[0]) + "\n").encode())
    for (power,), coefficient in poly.terms():
        value = Fraction(coefficient)
        digest.update(
            f"{power}:{value.numerator}/{value.denominator}\n".encode()
        )
    return digest.hexdigest()


def derive_certificate_data():
    x, y, q, r = sp.symbols("x y q r")

    # F_j(p)=P(Bin(4,p)<=j).  Here x=p_0 and y=p_1, so
    # F_0(x)=F_1(y)=alpha.  At equality of the next two Stringer spacings,
    # p_2=q=2*y-x and F_2(q)=alpha.
    F0 = (1 - x) ** 4
    F1 = (1 - y) ** 3 * (1 + 3 * y)
    F2 = (1 - q) ** 2 * (1 + 2 * q + 3 * q**2)
    relation = sp.Poly(sp.expand(F0 - F1), x, y)
    spacing_equality = sp.Poly(
        sp.expand(F1 - F2.subs(q, 2 * y - x)), x, y
    )
    resultant = sp.factor(
        sp.resultant(relation.as_expr(), spacing_equality.as_expr(), y)
    )
    resultant_factor = sp.Poly(
        43209 * x**8
        - 432090 * x**7
        + 1912545 * x**6
        - 4833544 * x**5
        + 7534696 * x**4
        - 7294848 * x**3
        + 4152976 * x**2
        - 1188352 * x
        + 110592,
        x,
    )
    expected_resultant = 144 * x**2 * (x - 1) ** 6 * resultant_factor.as_expr()
    if sp.expand(resultant - expected_resultant) != 0:
        raise AssertionError("middle-spacing resultant changed")

    first_interval = (sp.Rational(411, 2500), sp.Rational(329, 2000))
    second_interval = (sp.Rational(139, 200), sp.Rational(6951, 10000))
    if sp.count_roots(resultant_factor, 0, 1) != 2:
        raise AssertionError("resultant no longer has exactly two unit roots")
    if sp.count_roots(resultant_factor, *first_interval) != 1:
        raise AssertionError("first resultant root lost its isolating interval")
    if sp.count_roots(resultant_factor, *second_interval) != 1:
        raise AssertionError("second resultant root lost its isolating interval")
    if not (
        resultant_factor.eval(first_interval[0])
        * resultant_factor.eval(first_interval[1])
        < 0
        and resultant_factor.eval(second_interval[0])
        * resultant_factor.eval(second_interval[1])
        < 0
    ):
        raise AssertionError("resultant isolating signs changed")

    # The middle comparison q=2*y-x is automatically favorable when q>=1.
    # Along F0(x)=F1(y), its q=1 boundary is exactly x=11/19.
    q_one_identity = sp.factor(
        F0 - F1.subs(y, (1 + x) / 2)
    )
    expected_q_one = (x - 1) ** 3 * (19 * x - 11) / 16
    if sp.expand(q_one_identity - expected_q_one) != 0:
        raise AssertionError("q=1 boundary identity changed")

    # The terminal critical point has r=alpha^(1/4) and x=1-r.
    critical = sp.Poly(
        sp.expand((2 * r - 1) ** 3 * (7 - 6 * r) - r**4), r
    )
    nontrivial = sp.Poly(49 * r**3 - 79 * r**2 + 41 * r - 7, r)
    if sp.expand(critical.as_expr() + (r - 1) * nontrivial.as_expr()) != 0:
        raise AssertionError("n=4 terminal factorization changed")
    coarse_r_interval = (
        sp.Rational(7529, 10000),
        sp.Rational(753, 1000),
    )
    if sp.count_roots(nontrivial, sp.Rational(1, 2), sp.Rational(7, 8)) != 1:
        raise AssertionError("terminal cubic lost uniqueness")
    if not (
        nontrivial.eval(coarse_r_interval[0]) < 0
        < nontrivial.eval(coarse_r_interval[1])
    ):
        raise AssertionError("terminal cubic bracket signs changed")

    # Rational comparisons used in the written domain exclusion.
    if not (
        1 - coarse_r_interval[1] >= sp.Rational(247, 1000)
        and first_interval[1] < sp.Rational(247, 1000)
        and sp.Rational(11, 19) < second_interval[0]
    ):
        raise AssertionError("low-alpha root exclusion margins changed")
    if not (
        sp.Rational(21, 25) ** 4 < sp.Rational(1, 2)
        and sp.Rational(4, 25) < first_interval[0]
    ):
        raise AssertionError("high-alpha root exclusion margins changed")

    r_lower, r_upper = coarse_r_interval
    for _ in range(80):
        midpoint = (r_lower + r_upper) / 2
        if nontrivial.eval(midpoint) < 0:
            r_lower = midpoint
        else:
            r_upper = midpoint

    return {
        "relation": relation,
        "spacing_equality": spacing_equality,
        "resultant_factor": resultant_factor,
        "resultant": resultant,
        "q_one_identity": q_one_identity,
        "critical": critical,
        "nontrivial": nontrivial,
        "first_interval": first_interval,
        "second_interval": second_interval,
        "coarse_r_interval": coarse_r_interval,
        "r_interval": (r_lower, r_upper),
    }


def build_certificate() -> dict[str, object]:
    data = derive_certificate_data()
    P = data["resultant_factor"]
    first = data["first_interval"]
    second = data["second_interval"]
    r_interval = data["r_interval"]
    alpha_interval = (r_interval[0] ** 4, r_interval[1] ** 4)
    display_interval = (
        sp.Rational(3214292869970077, 10**16),
        sp.Rational(3214292869970078, 10**16),
    )
    if not (
        display_interval[0]
        < alpha_interval[0]
        < alpha_interval[1]
        < display_interval[1]
    ):
        raise AssertionError("displayed alpha_4 bracket is not rigorous")
    return {
        "schema_version": 1,
        "mathematical_status": (
            "ordinary monotonicity reduction with exact resultant and Sturm certificate"
        ),
        "scope": {
            "hypothesis": "Four-trial Clopper--Pearson factors at 0<alpha<=alpha_4.",
            "conclusion": (
                "The five ascending-knot Stringer threshold weights are "
                "nondecreasing. The endpoint alpha_4 is sharp for this property "
                "and for pointwise Stringer--Gaffke domination."
            ),
        },
        "critical_endpoint": {
            "definition": (
                "alpha_4=r_4^4, where r_4 is the unique root in (1/2,7/8) "
                "of 49*r^3-79*r^2+41*r-7"
            ),
            "r_lower": _fraction_record(r_interval[0]),
            "r_upper": _fraction_record(r_interval[1]),
            "alpha_lower": _fraction_record(alpha_interval[0]),
            "alpha_upper": _fraction_record(alpha_interval[1]),
            "alpha_display_interval": (
                "[0.3214292869970077, 0.3214292869970078]"
            ),
            "alpha_display_interval_verified_exactly": True,
            "lower_sign": "negative",
            "upper_sign": "positive",
        },
        "middle_spacing_resultant": {
            "meaning": (
                "If F_0(x)=F_1(y)=F_2(2y-x), then "
                "144*x^2*(x-1)^6*P(x)=0."
            ),
            "degree": P.degree(),
            "coefficients_descending": [
                str(Fraction(value)) for value in P.all_coeffs()
            ],
            "sha256": _polynomial_digest(P),
            "roots_in_open_unit_interval": 2,
            "isolating_intervals": [
                {
                    "lower": _fraction_record(first[0]),
                    "upper": _fraction_record(first[1]),
                    "lower_sign": str(sp.sign(P.eval(first[0]))),
                    "upper_sign": str(sp.sign(P.eval(first[1]))),
                    "root_count": 1,
                },
                {
                    "lower": _fraction_record(second[0]),
                    "upper": _fraction_record(second[1]),
                    "lower_sign": str(sp.sign(P.eval(second[0]))),
                    "upper_sign": str(sp.sign(P.eval(second[1]))),
                    "root_count": 1,
                },
            ],
        },
        "domain_exclusion": {
            "q_equals_one_identity": (
                "F_0(x)-F_1((1+x)/2)=(x-1)^3*(19*x-11)/16"
            ),
            "low_alpha": (
                "alpha<=alpha_4 gives x>247/1000. If 2y-x<1, then "
                "x<11/19. The resultant has no root in that interval."
            ),
            "high_alpha": (
                "alpha>=1/2 gives x<4/25, while the first resultant root "
                "is above 411/2500. Thus the middle spacing has no zero there."
            ),
            "symmetry": (
                "p_j(1-alpha)=1-p_{3-j}(alpha), so the remaining "
                "nonterminal spacings are reflected copies."
            ),
        },
        "arithmetic": {
            "resultant": "SymPy exact integer polynomial elimination",
            "root_counts": "SymPy exact Sturm counts over rational intervals",
            "floating_point_role": (
                "none; the displayed decimal interval denotes exact rationals"
            ),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        "n=4 Stringer weights certified monotone through "
        "alpha_4 in [0.3214292869970077, 0.3214292869970078]"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
