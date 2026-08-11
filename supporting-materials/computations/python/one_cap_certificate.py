"""Exact certificate for the one-upper-knot Stringer--Gaffke comparison.

Let ``x_0 <= ... <= x_n = 1`` be knots, let ``c_i > 0`` sum to one,
and put ``s = sum(c_i*x_i)``.  When ``s >= x_(n-1)``, the upper tail of
the uniform Dirichlet average has only one active spline knot.  A concavity
argument then gives

    P(sum(x_i*D_i) > s)
        <= max_{1 <= r <= n} (C_r/(C_r+c_n))**r,

where ``C_r = c_0 + ... + c_(r-1)``.

For the ascending-knot binomial Stringer weights,

    C_r = 1 - p_n(n-r),       c_n = p_n(0).

This program proves that every term in the displayed maximum is at most
``alpha`` for ``1 <= n <= 200`` and ``alpha`` in ``{0.01,0.05,0.10}``.
The ``r=n`` term is exactly alpha by the defining equation for ``p_n(0)``.
For ``r<n``, all Clopper--Pearson roots are enclosed on a dyadic grid and
their binomial-CDF endpoint signs are checked with integer arithmetic.
The resulting rational upper bound is then compared exactly with alpha.

Floating point is used only to render readable decimal summaries.  The
certificate is an independent finite-range proof and regression for the
stronger analytic all-``n`` theorem.  It is not a proof of the unrestricted
Stringer conjecture.
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


ALPHAS = ("0.01", "0.05", "0.10")
FACTOR_BITS = 64
N_MAX = 200
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE.parent / "certificates" / "one-cap-certificate.json"


# The worst-case rational powers at n=200 can exceed Python's defensive
# integer-to-decimal conversion limit.  Those integers are intentional
# certificate output, not untrusted input.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


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


def _update_integer(digest: "hashlib._Hash", value: int) -> None:
    """Add a signed integer to a canonical length-delimited digest."""

    sign = b"-" if value < 0 else b"+"
    magnitude = abs(value)
    encoded = magnitude.to_bytes(
        max(1, (magnitude.bit_length() + 7) // 8), "big"
    )
    digest.update(sign)
    digest.update(len(encoded).to_bytes(8, "big"))
    digest.update(encoded)


def _update_fraction(digest: "hashlib._Hash", value: Fraction) -> None:
    value = Fraction(value)
    _update_integer(digest, value.numerator)
    _update_integer(digest, value.denominator)


def certify_alpha(alpha: str, n_max: int = N_MAX,
                  factor_bits: int = FACTOR_BITS) -> dict[str, object]:
    """Certify every nonterminal cap vertex for one alpha value."""

    target = Fraction(alpha)
    digest = hashlib.sha256()
    checked = 0
    worst: tuple[Fraction, int, int, Fraction] | None = None
    maximum_width = Fraction(0)

    for n in range(1, n_max + 1):
        brackets = exact_binomial_factor_brackets(n, alpha, factor_bits)
        if len(brackets) != n + 1:
            raise AssertionError("factor bracket count is inconsistent")

        for lower, upper in brackets:
            if not lower <= upper:
                raise AssertionError("factor bracket endpoints are reversed")
            maximum_width = max(maximum_width, upper - lower)

        # The CP roots are strictly increasing.  Enclosing adjacent roots
        # in disjoint intervals supplies an additional exact regression
        # check on the ascending-knot Stringer weights used in the proof.
        for j in range(n):
            if not brackets[j][1] < brackets[j + 1][0]:
                raise AssertionError(
                    f"factor monotonicity was not certified: alpha={alpha}, "
                    f"n={n}, j={j}"
                )

        p0_lower = brackets[0][0]
        if p0_lower <= 0:
            raise AssertionError("the p_0 lower endpoint must be positive")

        # r=n is analytic equality:
        #   C_n=1-p_0 and (C_n/(C_n+p_0))^n=(1-p_0)^n=alpha.
        # For r<n, C_r=1-p_n(n-r).  Taking the upper endpoint of C_r
        # and the lower endpoint of p_0 gives a rigorous upper bound because
        # C/(C+p_0) increases in C and decreases in p_0.
        for r in range(1, n):
            c_upper = 1 - brackets[n - r][0]
            upper = (c_upper / (c_upper + p0_lower)) ** r
            if not upper < target:
                raise AssertionError(
                    "one-cap inequality was not certified: "
                    f"alpha={alpha}, n={n}, r={r}, upper={upper}"
                )

            checked += 1
            _update_integer(digest, n)
            _update_integer(digest, r)
            _update_fraction(digest, upper)
            margin = target - upper
            if worst is None or upper > worst[0]:
                worst = (upper, n, r, margin)

    # n=1 has no nonterminal vertex, but n>=2 in the committed range does.
    if n_max >= 2 and worst is None:
        raise AssertionError("no nonterminal cap vertices were checked")

    worst_record = None
    if worst is not None:
        upper, n, r, margin = worst
        worst_record = {
            "n": n,
            "r": r,
            "certified_upper_bound": _fraction_record(upper),
            "strict_margin_below_alpha": _fraction_record(margin),
            "upper_bound_divided_by_alpha": _fraction_record(upper / target),
        }

    return {
        "alpha": alpha,
        "nominal_confidence": _decimal(1 - target, 2),
        "sample_sizes_checked": f"1 through {n_max}",
        "strict_nonterminal_inequalities_checked": checked,
        "analytic_terminal_equalities": n_max,
        "worst_nonterminal_case": worst_record,
        "sha256_all_nonterminal_upper_bounds": digest.hexdigest(),
        "maximum_factor_bracket_width": _fraction_record(maximum_width),
    }


def build_certificate(n_max: int = N_MAX,
                      factor_bits: int = FACTOR_BITS) -> dict[str, object]:
    if n_max < 1:
        raise ValueError("n_max must be positive")
    if factor_bits < 1:
        raise ValueError("factor_bits must be positive")

    return {
        "schema_version": 1,
        "status": (
            "Exact finite-range certificate for the one-upper-knot region; "
            "not an unrestricted coverage certificate."
        ),
        "theorem": {
            "sample_sizes": f"1 <= n <= {n_max}",
            "alphas": list(ALPHAS),
            "sample_condition": (
                "The binomial Stringer value is at least the largest "
                "observed taint."
            ),
            "conclusion": (
                "The binomial Stringer value is at least the valid Gaffke "
                "upper limit; hence the pre-specified Stringer-Gaffke "
                "safeguard has zero uplift on that sample."
            ),
            "analytic_cap_bound": (
                "P(sum_i x_i D_i > sum_i c_i x_i) <= "
                "max_{1<=r<=n} (C_r/(C_r+c_n))^r when the threshold is "
                "at least x_(n-1), where C_r=sum_{i<r} c_i."
            ),
            "stringer_specialization": (
                "C_r=1-p_n(n-r), c_n=p_n(0), and the r=n term equals "
                "alpha analytically."
            ),
        },
        "arithmetic": {
            "factor_bracket_bits": factor_bits,
            "factor_endpoint_signs": (
                "exact integer binomial-CDF comparisons"
            ),
            "cap_comparisons": "exact Fraction arithmetic",
            "floating_point_role": "decimal rendering only",
        },
        "results": [
            certify_alpha(alpha, n_max, factor_bits) for alpha in ALPHAS
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
    with args.out.open("w") as handle:
        json.dump(certificate, handle, indent=2, sort_keys=True)
        handle.write("\n")

    for result in certificate["results"]:
        worst = result["worst_nonterminal_case"]
        if worst is None:
            print("alpha=%s: no nonterminal vertices" % result["alpha"])
        else:
            print(
                "alpha=%s: %d strict vertices; worst n=%d r=%d upper=%s "
                "margin=%s"
                % (
                    result["alpha"],
                    result["strict_nonterminal_inequalities_checked"],
                    worst["n"],
                    worst["r"],
                    worst["certified_upper_bound"]["decimal"],
                    worst["strict_margin_below_alpha"]["decimal"],
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
