"""Exact checks for the open Dirichlet--Poissonization comparison.

The all-sample-size Poisson program compares

    T_y = sum_i y_i D_i,       Z_y = sum_i y_i E_i,

where ``D`` is uniform on the ``n``-simplex and the ``E_i`` are independent
unit exponentials.  At the normalized threshold used in the research note,
the desired inequality is

    P(Z_y > n) >= P(T_y > 1).

This module does *not* claim that inequality in full.  It records two exact
pieces of progress:

* the comparison is true for every equal-block vector (a common positive
  coordinate repeated ``k`` times, followed by zeros), by the published
  Anderson--Samuels binomial--Poisson inequality; and
* a fully explicit ``1/n``-concave law shows why generic one-dimensional
  s-concave localization is too broad to prove the comparison at alpha=1/2.

The latter check is formal.  Its probability above one is rational, and its
gamma-smoothed probability is a rational multiple of ``exp(-15/2)``.  The
exponential is enclosed by exact rational alternating-series bounds from
``stringer.exact_exp_neg_bounds``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from stringer import EXACT_EXP_PAIRS, exact_exp_neg_bounds


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "dirichlet-poissonization-certificate.json"
)


def verify_saffine_symbolic_identities() -> None:
    """Derive the two closed forms used in the localization obstruction."""
    from sympy import Rational, diff, exp, factorial, integrate, simplify, symbols

    n = 10
    y = symbols("y", positive=True)
    upper = Rational(4, 3)
    raw_density = (Rational(3, 4) + Rational(3, 16) * y) ** (n - 1)
    normalizer = integrate(raw_density, (y, 0, upper))
    assert normalizer == Rational(989527, 1966080)
    density = raw_density / normalizer

    tail = simplify(integrate(density, (y, 1, upper)))
    assert tail == Rational(47532839741, 94326751232)

    # H is C^(n-1) at one.  Integrating density*H^(n)/n! by parts n
    # times therefore leaves this endpoint expression; every lower-endpoint
    # term vanishes because y^n exp(-n/y) is flat at zero.
    h_upper_piece = y**n * exp(-Rational(n, 1) / y) - (y - 1) ** n
    boundary = sum(
        (-1) ** order
        * diff(density, y, order)
        * diff(h_upper_piece, y, n - 1 - order)
        for order in range(n)
    )
    smoothed_minus_tail = simplify(
        boundary.subs(y, upper) / factorial(n)
    )
    expected_smoothed = (
        Rational(6134560249, 6926689) * exp(-Rational(15, 2))
    )
    assert simplify(tail + smoothed_minus_tail - expected_smoothed) == 0


def _decimal(value: Fraction, digits: int = 18) -> str:
    value = Fraction(value)
    with localcontext() as context:
        context.prec = digits + 30
        result = Decimal(value.numerator) / Decimal(value.denominator)
        return format(result, f".{digits}f")


def _fraction_record(value: Fraction, digits: int = 18) -> dict[str, str]:
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": _decimal(value, digits),
    }


def _large_fraction_record(value: Fraction, digits: int = 18) -> dict[str, object]:
    """Record a huge exact rational compactly but reproducibly.

    Alternating-series bounds followed by repeated squaring can have tens of
    thousands of decimal digits.  Embedding each integer in JSON would make
    the research certificate needlessly large.  The source regenerates the
    exact value; this record preserves its bit lengths, a canonical binary
    SHA-256 digest, and a readable decimal summary.
    """
    value = Fraction(value)
    digest = hashlib.sha256()
    for integer in (value.numerator, value.denominator):
        sign = b"-" if integer < 0 else b"+"
        magnitude = abs(integer)
        encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
        digest.update(sign)
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return {
        "numerator_bits": abs(value.numerator).bit_length(),
        "denominator_bits": value.denominator.bit_length(),
        "sha256_canonical_binary": digest.hexdigest(),
        "decimal": _decimal(value, digits),
    }


def _binomial_lower_tail(n: int, j: int, p: Fraction) -> Fraction:
    """Exact ``P(Bin(n,p) <= j)`` for rational ``p``."""
    if not 0 <= j <= n:
        raise ValueError("j must lie in {0,...,n}")
    if not 0 <= p <= 1:
        raise ValueError("p must lie in [0,1]")
    return sum(
        Fraction(math.comb(n, ell))
        * p**ell
        * (1 - p) ** (n - ell)
        for ell in range(j + 1)
    )


def equal_block_check(n: int, k: int, lam: Fraction) -> dict[str, object]:
    """Certify one rational equal-block binomial--Poisson comparison.

    The coefficient vector has ``k`` entries equal to ``n/lam`` and
    ``n+1-k`` zero entries.  The condition ``lam>=k`` is exactly the
    coefficient-sum condition.  The Dirichlet tail is a binomial lower tail;
    the exponential tail is a Poisson lower tail.

    This function is only a regression check.  The proof for all real
    ``lam>=k`` is the Anderson--Samuels theorem plus continuity at ``lam=k``.
    """
    lam = Fraction(lam)
    if n < 1 or not 1 <= k <= n:
        raise ValueError("require n>=1 and 1<=k<=n")
    if not Fraction(k) <= lam <= n:
        raise ValueError("require k<=lam<=n")

    p = lam / n
    dirichlet_tail = _binomial_lower_tail(n, k - 1, p)
    polynomial = sum(lam**ell / math.factorial(ell) for ell in range(k))
    exp_lower, exp_upper = exact_exp_neg_bounds(lam, EXACT_EXP_PAIRS)
    poisson_lower = polynomial * exp_lower
    poisson_upper = polynomial * exp_upper
    if poisson_lower < dirichlet_tail:
        raise AssertionError("equal-block comparison was not certified")
    return {
        "n": n,
        "positive_coordinates": k,
        "lambda": _fraction_record(lam),
        "common_positive_coordinate": _fraction_record(Fraction(n, 1) / lam),
        "dirichlet_tail": _fraction_record(dirichlet_tail),
        "poisson_tail_lower_bound": _large_fraction_record(poisson_lower),
        "poisson_tail_upper_bound": _large_fraction_record(poisson_upper),
        "certified_lower_margin": _large_fraction_record(
            poisson_lower - dirichlet_tail
        ),
    }


def saffine_localization_obstruction() -> dict[str, object]:
    r"""Return an exact obstruction to the over-broad localization route.

    Put ``n=10`` and give ``Y`` the density proportional to

        (3/4 + 3y/16)^9,       0 <= y <= 4/3.

    This is a ``1/10``-concave probability law because its density to the
    power ``1/9`` is affine.  Direct integration gives

        P(Y>1) = 47532839741 / 94326751232 > 1/2,

    while, for ``S ~ Gamma(11,1)`` independent of ``Y``,

        P(SY>10) = (6134560249 / 6926689) exp(-15/2) < 1/2.

    Thus checking every s-affine localization extremizer cannot establish
    the desired implication at alpha=1/2.  The actual Dirichlet-average
    class is narrower and is not refuted by this example.
    """
    n = 10
    support_upper = Fraction(4, 3)
    density_normalizer = Fraction(989527, 1966080)
    tail = Fraction(47532839741, 94326751232)
    smoothed_coefficient = Fraction(6134560249, 6926689)
    exp_lower, exp_upper = exact_exp_neg_bounds(Fraction(15, 2))
    smoothed_lower = smoothed_coefficient * exp_lower
    smoothed_upper = smoothed_coefficient * exp_upper
    half = Fraction(1, 2)

    if not tail > half:
        raise AssertionError("the exact s-affine tail is not above one half")
    if not smoothed_upper < half:
        raise AssertionError("the smoothed tail is not certified below one half")
    if not smoothed_upper < tail:
        raise AssertionError("the localization obstruction was not certified")

    return {
        "n": n,
        "law": {
            "support": ["0", str(support_upper)],
            "unnormalized_density": "(3/4 + 3*y/16)^9",
            "density_integral": _fraction_record(density_normalizer),
            "s_concavity_parameter": "1/10",
        },
        "tail_probability_P_Y_gt_1": _fraction_record(tail),
        "smoothed_probability": {
            "identity": "P(S*Y>10)=(6134560249/6926689)*exp(-15/2)",
            "coefficient": _fraction_record(smoothed_coefficient),
            "exp_argument": "15/2",
            "exp_lower_bound": _large_fraction_record(exp_lower),
            "exp_upper_bound": _large_fraction_record(exp_upper),
            "probability_lower_bound": _large_fraction_record(smoothed_lower),
            "probability_upper_bound": _large_fraction_record(smoothed_upper),
        },
        "exact_comparisons": {
            "P_Y_gt_1": "> 1/2",
            "P_S_Y_gt_10": "< 1/2",
            "conclusion": (
                "The generic s-affine localization class is too broad; "
                "this is not a counterexample for a Dirichlet average."
            ),
        },
    }


def build_certificate() -> dict[str, object]:
    verify_saffine_symbolic_identities()
    checks = []
    # Include both the boundary lambda=k and strict lambda>k cases.  These
    # finite checks guard the beta/binomial and gamma/Poisson translations;
    # they are not a substitute for the cited all-parameter theorem.
    for n in (2, 3, 5, 10, 20):
        for k in range(1, n + 1):
            checks.append(equal_block_check(n, k, Fraction(k)))
            if k < n:
                checks.append(equal_block_check(n, k, Fraction(k + n, 2)))

    return {
        "schema_version": 1,
        "status": (
            "Research certificate: exact reductions and an obstruction, "
            "not an all-sample-size coverage certificate."
        ),
        "divided_difference_target": (
            "For y_i>=0 with sum_i y_i<=n, prove "
            "[y_0,...,y_n]H_n>=0, where "
            "H_n(y)=y^n exp(-n/y)-(y-1)_+^n."
        ),
        "equal_block_profiles": {
            "analytic_basis": (
                "Anderson--Samuels (1967), with continuity at lambda=k"
            ),
            "regression_checks": checks,
        },
        "localization_obstruction": saffine_localization_obstruction(),
        "exponential_series_pairs": EXACT_EXP_PAIRS,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args(argv)
    certificate = build_certificate()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    obstruction = certificate["localization_obstruction"]
    print("Dirichlet--Poissonization reduction checks: PASS")
    print(
        "equal-block rational regression cases:",
        len(certificate["equal_block_profiles"]["regression_checks"]),
    )
    print(
        "s-affine obstruction: P(Y>1)=%s; P(SY>10)<%s"
        % (
            obstruction["tail_probability_P_Y_gt_1"]["decimal"],
            obstruction["smoothed_probability"]["probability_upper_bound"][
                "decimal"
            ],
        )
    )
    print("general Dirichlet-average inequality: OPEN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
