"""Exact algebraic certificate for the sharp monotone-weight range at n=5.

Let ``p_j(alpha)`` be the one-sided Clopper--Pearson upper limits for five
trials and define the ascending-knot Stringer weights

    (1-p_4, p_4-p_3, p_3-p_2, p_2-p_1, p_1-p_0, p_0).

The terminal edge changes order at ``alpha_5=r_5**5``, where ``r_5`` is the
unique root in ``(1/2, 9/10)`` of

    129*r**4 - 271*r**3 + 209*r**2 - 71*r + 9 = 0.

This module certifies that none of the other four spacings changes sign
first.  Two exact eliminations reduce the nonterminal comparisons to an
eleventh-degree factor and to a fifth-degree plus an eighth-degree factor.
Exact Sturm counts, rational domain exclusions, and rational sign witnesses
then prove that all six weights are nondecreasing precisely through
``alpha_5``.  No floating-point sign enters the certificate.
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
    HERE.parent / "certificates" / "n5-weight-monotonicity-certificate.json"
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


def _tail(j: int, value, *, trials: int = 5):
    """Return P(Bin(trials,value)<=j) as an exact expression."""

    return sp.expand(
        sum(
            sp.binomial(trials, k)
            * value**k
            * (1 - value) ** (trials - k)
            for k in range(j + 1)
        )
    )


def derive_certificate_data():
    x, y, r = sp.symbols("x y r")

    # d_3=2*p_1-p_0-p_2.  At equality put p_0=x, p_1=y, and
    # p_2=2*y-x, then eliminate y from the two common-tail equations.
    d3_relation = sp.Poly(_tail(0, x) - _tail(1, y), x, y)
    d3_equality = sp.Poly(
        _tail(1, y) - _tail(2, 2 * y - x), x, y
    )
    d3_resultant = sp.factor(
        sp.resultant(d3_relation.as_expr(), d3_equality.as_expr(), y)
    )
    d3_factor = sp.Poly(
        682456 * x**11
        - 8871928 * x**10
        + 53222168 * x**9
        - 193600856 * x**8
        + 472259515 * x**7
        - 807107237 * x**6
        + 980133686 * x**5
        - 838362386 * x**4
        + 487990767 * x**3
        - 179494185 * x**2
        + 35660800 * x
        - 2508800,
        x,
    )
    expected_d3_resultant = (
        400000 * x**2 * (x - 1) ** 12 * d3_factor.as_expr()
    )
    if sp.expand(d3_resultant - expected_d3_resultant) != 0:
        raise AssertionError("d3 resultant changed")

    d3_root_interval = (
        sp.Rational(1356, 10000),
        sp.Rational(1357, 10000),
    )
    if sp.count_roots(d3_factor, 0, 1) != 1:
        raise AssertionError("d3 factor no longer has one unit root")
    if sp.count_roots(d3_factor, *d3_root_interval) != 1:
        raise AssertionError("d3 root lost its isolating interval")
    if not (
        d3_factor.eval(d3_root_interval[0]) < 0
        < d3_factor.eval(d3_root_interval[1])
    ):
        raise AssertionError("d3 isolating signs changed")

    # d_2=2*p_2-p_1-p_3.  At equality put p_1=x, p_2=y, and
    # p_3=2*y-x.  One factor is exactly the alpha=1/2 tail equation;
    # the other has one unit root, on the inadmissible q>=1 branch.
    d2_relation = sp.Poly(_tail(1, x) - _tail(2, y), x, y)
    d2_equality = sp.Poly(
        _tail(2, y) - _tail(3, 2 * y - x), x, y
    )
    d2_resultant = sp.factor(
        sp.resultant(d2_relation.as_expr(), d2_equality.as_expr(), y)
    )
    half_factor = sp.Poly(
        8 * x**5 - 30 * x**4 + 40 * x**3 - 20 * x**2 + 1,
        x,
    )
    extraneous_factor = sp.Poly(
        1597696 * x**8
        - 14379264 * x**7
        + 57099936 * x**6
        - 128866064 * x**5
        + 177551071 * x**4
        - 148524550 * x**3
        + 69152775 * x**2
        - 13296000 * x
        - 304000,
        x,
    )
    expected_d2_resultant = (
        1600000
        * x**6
        * (x - 1) ** 6
        * half_factor.as_expr()
        * extraneous_factor.as_expr()
    )
    if sp.expand(d2_resultant - expected_d2_resultant) != 0:
        raise AssertionError("d2 resultant changed")
    if sp.expand(half_factor.as_expr() - (2 * _tail(1, x) - 1)) != 0:
        raise AssertionError("alpha=1/2 factor identity changed")

    extraneous_root_interval = (
        sp.Rational(8209, 10000),
        sp.Rational(821, 1000),
    )
    if sp.count_roots(extraneous_factor, 0, 1) != 1:
        raise AssertionError("extraneous factor no longer has one unit root")
    if sp.count_roots(extraneous_factor, *extraneous_root_interval) != 1:
        raise AssertionError("extraneous root lost its isolating interval")
    if not (
        extraneous_factor.eval(extraneous_root_interval[0]) < 0
        < extraneous_factor.eval(extraneous_root_interval[1])
    ):
        raise AssertionError("extraneous-root isolating signs changed")

    q_one_identity = sp.factor(
        _tail(1, x) - _tail(2, (1 + x) / 2)
    )
    expected_q_one = (x - 1) ** 3 * (67 * x**2 - 39 * x - 8) / 16
    if sp.expand(q_one_identity - expected_q_one) != 0:
        raise AssertionError("d2 q=1 boundary identity changed")
    if sp.expand((67 * x**2 - 39 * x - 8).subs(x, sp.Rational(3, 4))) <= 0:
        raise AssertionError("positive q=1 root is no longer below 3/4")
    if not extraneous_root_interval[0] > sp.Rational(3, 4):
        raise AssertionError("extraneous root exclusion margin changed")

    # Terminal critical point.  The displayed quartic is the nontrivial
    # factor of the all-dimensional terminal-edge equation at n=5.
    critical = sp.Poly(
        sp.expand((2 * r - 1) ** 4 * (9 - 8 * r) - r**5), r
    )
    nontrivial = sp.Poly(
        129 * r**4 - 271 * r**3 + 209 * r**2 - 71 * r + 9,
        r,
    )
    if sp.expand(
        critical.as_expr() + (r - 1) * nontrivial.as_expr()
    ) != 0:
        raise AssertionError("n=5 terminal factorization changed")
    coarse_r_interval = (
        sp.Rational(7932, 10000),
        sp.Rational(7933, 10000),
    )
    if sp.count_roots(
        nontrivial, sp.Rational(1, 2), sp.Rational(9, 10)
    ) != 1:
        raise AssertionError("terminal quartic lost uniqueness")
    if not (
        nontrivial.eval(coarse_r_interval[0]) < 0
        < nontrivial.eval(coarse_r_interval[1])
    ):
        raise AssertionError("terminal quartic bracket signs changed")

    # Exact domain margins used to exclude the d3 resultant root on both
    # the low-alpha interval and its reflected high-alpha interval.
    if not (
        1 - coarse_r_interval[1] == sp.Rational(2067, 10000)
        and sp.Rational(2067, 10000) > d3_root_interval[1]
    ):
        raise AssertionError("low-alpha d3 exclusion margin changed")
    if not coarse_r_interval[1] ** 5 < sp.Rational(1, 3):
        raise AssertionError("alpha_5<1/3 margin changed")
    if not sp.Rational(9, 10) ** 5 < sp.Rational(2, 3):
        raise AssertionError("reflected p0<1/10 margin changed")
    if not sp.Rational(1, 10) < d3_root_interval[0]:
        raise AssertionError("high-alpha d3 exclusion margin changed")

    # Rational sign witnesses.  At alpha=(3/4)^5, p0=1/4.
    alpha_witness = sp.Rational(3, 4) ** 5
    low_checks = {
        "F1_at_46_over_100_minus_alpha": sp.factor(
            _tail(1, sp.Rational(46, 100)) - alpha_witness
        ),
        "F2_at_65_over_100_minus_alpha": sp.factor(
            _tail(2, sp.Rational(65, 100)) - alpha_witness
        ),
        "F1_at_463_over_1000_minus_alpha": sp.factor(
            _tail(1, sp.Rational(463, 1000)) - alpha_witness
        ),
        "F2_at_648_over_1000_minus_alpha": sp.factor(
            _tail(2, sp.Rational(648, 1000)) - alpha_witness
        ),
        "F3_at_813_over_1000_minus_alpha": sp.factor(
            _tail(3, sp.Rational(813, 1000)) - alpha_witness
        ),
    }
    expected_low_signs = (1, -1, -1, 1, -1)
    if tuple(sp.sign(value) for value in low_checks.values()) != expected_low_signs:
        raise AssertionError("low-alpha rational sign witness changed")
    if not (
        2 * sp.Rational(46, 100)
        - sp.Rational(1, 4)
        - sp.Rational(65, 100)
        == sp.Rational(1, 50)
    ):
        raise AssertionError("d3 positive witness margin changed")
    if not (
        2 * sp.Rational(648, 1000)
        - sp.Rational(463, 1000)
        - sp.Rational(813, 1000)
        == sp.Rational(1, 50)
    ):
        raise AssertionError("d2 positive witness margin changed")

    # At beta=(19/20)^5, p0=1/20 and d3 is strictly negative.
    beta_witness = sp.Rational(19, 20) ** 5
    high_checks = {
        "F1_at_182_over_1000_minus_beta": sp.factor(
            _tail(1, sp.Rational(182, 1000)) - beta_witness
        ),
        "F2_at_344_over_1000_minus_beta": sp.factor(
            _tail(2, sp.Rational(344, 1000)) - beta_witness
        ),
    }
    if tuple(sp.sign(value) for value in high_checks.values()) != (-1, 1):
        raise AssertionError("high-alpha rational sign witness changed")
    if not (
        2 * sp.Rational(182, 1000)
        - sp.Rational(1, 20)
        - sp.Rational(344, 1000)
        == -sp.Rational(3, 100)
    ):
        raise AssertionError("d3 negative witness margin changed")
    if not (
        alpha_witness < coarse_r_interval[0] ** 5
        and beta_witness > 1 - coarse_r_interval[0] ** 5
    ):
        raise AssertionError("sign witnesses left their required domains")

    r_lower, r_upper = coarse_r_interval
    for _ in range(90):
        midpoint = (r_lower + r_upper) / 2
        if nontrivial.eval(midpoint) < 0:
            r_lower = midpoint
        else:
            r_upper = midpoint

    return {
        "d3_relation": d3_relation,
        "d3_equality": d3_equality,
        "d3_resultant": d3_resultant,
        "d3_factor": d3_factor,
        "d3_root_interval": d3_root_interval,
        "d2_relation": d2_relation,
        "d2_equality": d2_equality,
        "d2_resultant": d2_resultant,
        "half_factor": half_factor,
        "extraneous_factor": extraneous_factor,
        "extraneous_root_interval": extraneous_root_interval,
        "q_one_identity": q_one_identity,
        "critical": critical,
        "nontrivial": nontrivial,
        "coarse_r_interval": coarse_r_interval,
        "r_interval": (r_lower, r_upper),
        "alpha_witness": alpha_witness,
        "beta_witness": beta_witness,
        "low_checks": low_checks,
        "high_checks": high_checks,
    }


def build_certificate() -> dict[str, object]:
    data = derive_certificate_data()
    r_interval = data["r_interval"]
    alpha_interval = (r_interval[0] ** 5, r_interval[1] ** 5)
    display_interval = (
        sp.Rational(3141689898050253, 10**16),
        sp.Rational(3141689898050254, 10**16),
    )
    if not (
        display_interval[0]
        < alpha_interval[0]
        < alpha_interval[1]
        < display_interval[1]
    ):
        raise AssertionError("displayed alpha_5 bracket is not rigorous")

    d3_factor = data["d3_factor"]
    d3_root = data["d3_root_interval"]
    extra = data["extraneous_factor"]
    extra_root = data["extraneous_root_interval"]
    return {
        "schema_version": 1,
        "mathematical_status": (
            "ordinary monotonicity reduction with exact resultants, Sturm "
            "counts, domain exclusions, and rational sign witnesses"
        ),
        "scope": {
            "hypothesis": "Five-trial Clopper--Pearson factors at 0<alpha<=alpha_5.",
            "conclusion": (
                "The six ascending-knot Stringer threshold weights are "
                "nondecreasing. The endpoint alpha_5 is sharp for this property "
                "and for the terminal-edge pointwise comparison mechanism."
            ),
            "dependency": (
                'Coverage throughout the weight-order interval would follow from the unresolved six-coordinate monotone cap theorem; the exact weight calculation alone does not prove it. The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.'
            ),
        },
        "critical_endpoint": {
            "definition": (
                "alpha_5=r_5^5, where r_5 is the unique root in (1/2,9/10) "
                "of 129*r^4-271*r^3+209*r^2-71*r+9"
            ),
            "r_lower": _fraction_record(r_interval[0]),
            "r_upper": _fraction_record(r_interval[1]),
            "alpha_lower": _fraction_record(alpha_interval[0]),
            "alpha_upper": _fraction_record(alpha_interval[1]),
            "alpha_display_interval": (
                "[0.3141689898050253, 0.3141689898050254]"
            ),
            "alpha_display_interval_verified_exactly": True,
            "confidence_display_interval": (
                "[0.6858310101949746, 0.6858310101949747]"
            ),
            "lower_sign": "negative",
            "upper_sign": "positive",
        },
        "d3_resultant": {
            "meaning": (
                "If F_0(x)=F_1(y)=F_2(2y-x), then "
                "400000*x^2*(x-1)^12*P_11(x)=0."
            ),
            "degree": d3_factor.degree(),
            "coefficients_descending": [
                str(Fraction(value)) for value in d3_factor.all_coeffs()
            ],
            "sha256": _polynomial_digest(d3_factor),
            "roots_in_open_unit_interval": 1,
            "isolating_interval": {
                "lower": _fraction_record(d3_root[0]),
                "upper": _fraction_record(d3_root[1]),
                "lower_sign": str(sp.sign(d3_factor.eval(d3_root[0]))),
                "upper_sign": str(sp.sign(d3_factor.eval(d3_root[1]))),
                "root_count": 1,
            },
            "domain_exclusion": (
                "For alpha<=alpha_5, p_0>2067/10000; for the reflected "
                "tail beta>=1-alpha_5>2/3, p_0<1/10. The sole root lies "
                "between 1356/10000 and 1357/10000."
            ),
        },
        "d2_resultant": {
            "meaning": (
                "If F_1(x)=F_2(y)=F_3(2y-x), then "
                "1600000*x^6*(x-1)^6*A_5(x)*B_8(x)=0."
            ),
            "half_factor": {
                "degree": data["half_factor"].degree(),
                "identity": "A_5(x)=2*F_1(x)-1",
                "sha256": _polynomial_digest(data["half_factor"]),
            },
            "extraneous_factor": {
                "degree": extra.degree(),
                "coefficients_descending": [
                    str(Fraction(value)) for value in extra.all_coeffs()
                ],
                "sha256": _polynomial_digest(extra),
                "roots_in_open_unit_interval": 1,
                "isolating_interval": {
                    "lower": _fraction_record(extra_root[0]),
                    "upper": _fraction_record(extra_root[1]),
                    "lower_sign": str(sp.sign(extra.eval(extra_root[0]))),
                    "upper_sign": str(sp.sign(extra.eval(extra_root[1]))),
                    "root_count": 1,
                },
            },
            "q_equals_one_identity": (
                "F_1(x)-F_2((1+x)/2)="
                "(x-1)^3*(67*x^2-39*x-8)/16"
            ),
            "domain_exclusion": (
                "At a genuine d_2 zero, q=2*p_2-p_1=p_3<1, which forces "
                "p_1 below the positive root of 67*x^2-39*x-8 and hence "
                "below 3/4. The sole B_8 root is above 8209/10000, while "
                "A_5=0 is exactly alpha=1/2."
            ),
        },
        "rational_sign_witnesses": {
            "low_alpha": {
                "alpha": _fraction_record(data["alpha_witness"]),
                "p0": _fraction_record(sp.Rational(1, 4)),
                "d3_lower_bound": _fraction_record(sp.Rational(1, 50)),
                "d2_lower_bound": _fraction_record(sp.Rational(1, 50)),
                "tail_differences": {
                    key: _fraction_record(value)
                    for key, value in data["low_checks"].items()
                },
            },
            "reflected_high_alpha": {
                "beta": _fraction_record(data["beta_witness"]),
                "p0": _fraction_record(sp.Rational(1, 20)),
                "d3_upper_bound": _fraction_record(-sp.Rational(3, 100)),
                "tail_differences": {
                    key: _fraction_record(value)
                    for key, value in data["high_checks"].items()
                },
            },
        },
        "symmetry": {
            "identity": "p_j(1-alpha)=1-p_(4-j)(alpha)",
            "spacing_relations": (
                "d_0(alpha)=-d_4(1-alpha), "
                "d_1(alpha)=-d_3(1-alpha), "
                "d_2(alpha)=-d_2(1-alpha)"
            ),
        },
        "arithmetic": {
            "resultants": "SymPy exact integer polynomial elimination",
            "root_counts": "SymPy exact Sturm counts over rational intervals",
            "sign_witnesses": "exact rational binomial-tail evaluations",
            "floating_point_role": (
                "none; displayed decimal intervals denote exact rationals"
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
        "n=5 Stringer weights certified monotone through "
        "alpha_5 in [0.3141689898050253, 0.3141689898050254]"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
