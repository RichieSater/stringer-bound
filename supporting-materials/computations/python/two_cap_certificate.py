"""Exact certificate for a two-cap Stringer--Gaffke collar through n=200.

See theory/TWO-CAP-COMPARISON.md for the analytic cap lemma and the proof
that each collar expression has no interior maximum.  This program checks
only the nontrivial endpoint q=q_0.  Clopper--Pearson factors are enclosed
dyadically with exact integer CDF signs; all subsequent inequalities are
exact Fraction comparisons.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from stringer import exact_binomial_factor_brackets


CONFIG = (
    ("0.10", Fraction(5, 12)),
    ("0.05", Fraction(1, 3)),
    ("0.01", Fraction(1, 4)),
)
FACTOR_BITS = 64
N_MAX = 200
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE.parent / "certificates" / "two-cap-certificate.json"

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _decimal(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits + 30
        result = Decimal(value.numerator) / Decimal(value.denominator)
        return format(result, f".{digits}f")


def _fraction_record(value: Fraction) -> dict[str, str]:
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": _decimal(value),
    }


def _update_integer(digest: "hashlib._Hash", value: int) -> None:
    sign = b"-" if value < 0 else b"+"
    magnitude = abs(value)
    encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    digest.update(sign)
    digest.update(len(encoded).to_bytes(8, "big"))
    digest.update(encoded)


def _update_fraction(digest: "hashlib._Hash", value: Fraction) -> None:
    _update_integer(digest, value.numerator)
    _update_integer(digest, value.denominator)


def geometric_sum(n: int, x: Fraction) -> Fraction:
    return sum((x ** k for k in range(n)), Fraction(0))


def derivative_coefficients(n: int, r: int, A: Fraction,
                            b: Fraction) -> tuple[Fraction, ...]:
    """Coefficients of (A+bx)S_n'(x)-rb S_n(x), trailing zero removed."""
    if n < 2 or not 1 <= r < n or A <= 0 or b <= 0:
        raise ValueError("require n>=2, 1<=r<n, and positive A,b")
    values = [A * (j + 1) + b * (j - r) for j in range(n - 1)]
    values.append(b * (n - 1 - r))
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def two_cap_tail_exact(knots: tuple[Fraction, ...], s: Fraction) -> Fraction:
    """Exact two-active-knot spline tail for distinct rational knots."""
    if len(knots) < 3 or tuple(sorted(knots)) != knots or knots[-1] != 1:
        raise ValueError("knots must be ordered and end at one")
    if len(set(knots)) != len(knots):
        raise ValueError("this helper requires distinct knots")
    n = len(knots) - 1
    if not knots[-3] <= s <= knots[-2] or s >= 1:
        raise ValueError("threshold is not in the two-cap region")
    top = (1 - s) ** n
    for value in knots[:-1]:
        top /= 1 - value
    second = (knots[-2] - s) ** n
    denominator = (1 - knots[-2])
    for value in knots[:-2]:
        denominator *= knots[-2] - value
    return top - second / denominator


def two_cap_uniform_bound(weights: tuple[Fraction, ...],
                          q: Fraction) -> Fraction:
    """Return the dimension-free right side of the two-cap lemma."""
    if len(weights) < 3 or any(value <= 0 for value in weights):
        raise ValueError("at least three positive weights are required")
    if sum(weights) != 1 or not 0 < q <= 1:
        raise ValueError("weights must sum to one and q must lie in (0,1]")
    n = len(weights) - 1
    K = weights[n] + weights[n - 1] * (1 - q)
    g = geometric_sum(n, 1 - q)
    prefix = Fraction(0)
    terms = []
    for r in range(1, n):
        prefix += weights[r - 1]
        terms.append(g * (prefix / (prefix + K)) ** r)
    return max(terms)


def certify_alpha(alpha: str, q_min: Fraction, n_max: int = N_MAX,
                  factor_bits: int = FACTOR_BITS) -> dict[str, object]:
    target = Fraction(alpha)
    x = 1 - q_min
    digest = hashlib.sha256()
    checked = 0
    worst: tuple[Fraction, int, int] | None = None
    maximum_width = Fraction(0)

    for n in range(2, n_max + 1):
        brackets = exact_binomial_factor_brackets(n, alpha, factor_bits)
        for lower, upper in brackets:
            if not lower <= upper:
                raise AssertionError("factor bracket endpoints are reversed")
            maximum_width = max(maximum_width, upper - lower)
        for j in range(n):
            if not brackets[j][1] < brackets[j + 1][0]:
                raise AssertionError("factor monotonicity was not certified")

        c_n_lower = brackets[0][0]
        c_nm1_lower = brackets[1][0] - brackets[0][1]
        if c_n_lower <= 0 or c_nm1_lower <= 0:
            raise AssertionError("positive terminal weights were not certified")
        g = geometric_sum(n, x)

        for r in range(1, n):
            C_upper = 1 - brackets[n - r][0]
            denominator_lower = C_upper + c_n_lower + c_nm1_lower * x
            upper = g * (C_upper / denominator_lower) ** r
            if not upper < target:
                raise AssertionError(
                    f"two-cap collar failed: alpha={alpha}, n={n}, r={r}"
                )
            checked += 1
            _update_integer(digest, n)
            _update_integer(digest, r)
            _update_fraction(digest, upper)
            if worst is None or upper > worst[0]:
                worst = (upper, n, r)

    if worst is None:
        raise AssertionError("no collar inequalities were checked")
    value, n, r = worst
    return {
        "alpha": alpha,
        "nominal_confidence": _decimal(1 - target, 2),
        "q_min": _fraction_record(q_min),
        "sample_sizes_checked": f"2 through {n_max}",
        "endpoint_inequalities_checked": checked,
        "worst_endpoint_case": {
            "n": n,
            "r": r,
            "certified_upper_bound": _fraction_record(value),
            "strict_margin_below_alpha": _fraction_record(target - value),
        },
        "maximum_factor_bracket_width": _fraction_record(maximum_width),
        "sha256_all_endpoint_upper_bounds": digest.hexdigest(),
    }


def build_certificate(n_max: int = N_MAX,
                      factor_bits: int = FACTOR_BITS) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": (
            "Exact finite-range certificate for the stated two-cap collars; "
            "not a certificate for the residual two-cap region or ordinary "
            "Stringer coverage."
        ),
        "theorem": {
            "sample_condition": (
                "The binomial Stringer threshold s lies between the second-"
                "largest and largest observed taints, and "
                "q=(1-largest_taint)/(1-s) is at least q_min."
            ),
            "conclusion": (
                "The binomial Stringer value is at least the valid Gaffke "
                "upper limit, so the pre-specified maximum has zero uplift."
            ),
            "analytic_bound": (
                "cap <= g_n(q)*max_{1<=r<n} "
                "(C_r/(C_r+c_n+c_(n-1)*(1-q)))^r"
            ),
            "collar_maximization": (
                "Each r-term has no interior maximum in x=1-q by a one-sign-"
                "change derivative polynomial; q=1 follows from the analytic "
                "one-cap theorem and this certificate checks q=q_min."
            ),
        },
        "arithmetic": {
            "factor_bracket_bits": factor_bits,
            "factor_endpoint_signs": "exact integer binomial-CDF comparisons",
            "collar_endpoint_comparisons": "exact Fraction arithmetic",
            "floating_point_role": "decimal rendering only",
        },
        "results": [
            certify_alpha(alpha, q_min, n_max, factor_bits)
            for alpha, q_min in CONFIG
        ],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-max", type=int, default=N_MAX)
    parser.add_argument("--factor-bits", type=int, default=FACTOR_BITS)
    args = parser.parse_args(argv)
    certificate = build_certificate(args.n_max, args.factor_bits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    for result in certificate["results"]:
        worst = result["worst_endpoint_case"]
        print(
            "alpha=%s q_min=%s: %d endpoints; worst n=%d r=%d upper=%s margin=%s"
            % (
                result["alpha"], result["q_min"]["decimal"],
                result["endpoint_inequalities_checked"], worst["n"], worst["r"],
                worst["certified_upper_bound"]["decimal"],
                worst["strict_margin_below_alpha"]["decimal"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
