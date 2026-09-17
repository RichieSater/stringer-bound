"""Exact certificate for the terminal Clopper--Pearson edge bifurcation.

For sample size n, let p_0 and p_1 be the one-sided Clopper--Pearson
upper limits after zero and one observed errors.  Along the ordered-knot edge

    (0, ..., 0, 1-epsilon, 1),

the cap derivative at epsilon=0 has the sign of p_1-2*p_0.  The nontrivial
critical level alpha_n is characterized by a single integer polynomial after
writing r=alpha**(1/n).  This module verifies the symbolic identities and
constructs exact rational brackets for the finite-n roots and the limiting
root.  Floating point is not used for any sign decision.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE.parent / "certificates" / "terminal-edge-bifurcation-certificate.json"


def critical_polynomial(n: int, r: Fraction) -> Fraction:
    """Return F_n(r); its nontrivial zero is r_n."""
    return (2 * r - 1) ** (n - 1) * (
        2 * n - 1 - 2 * (n - 1) * r
    ) - r**n


def bracket_finite_root(n: int, digits: int = 18) -> tuple[Fraction, Fraction]:
    """Bisect the unique nontrivial root using exact rational signs."""
    lo = Fraction(1, 2)
    hi = Fraction(2 * n - 1, 2 * n)  # unique maximum of the likelihood ratio
    if not (critical_polynomial(n, lo) < 0 < critical_polynomial(n, hi)):
        raise AssertionError("finite-root endpoint signs changed")
    target_width = Fraction(1, 10**digits)
    while hi - lo > target_width:
        mid = (lo + hi) / 2
        if critical_polynomial(n, mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def log_bounds(x: Fraction, terms: int = 160) -> tuple[Fraction, Fraction]:
    """Rigorous rational bounds for log(x), for positive x.

    Uses log(x)=2*atanh((x-1)/(x+1)).  In the range used here the argument is
    negative, so the omitted tail has a known sign and a geometric bound.
    """
    if x <= 0:
        raise ValueError("log input must be positive")
    z = (x - 1) / (x + 1)
    total = Fraction(0)
    power = z
    for k in range(terms):
        total += power / (2 * k + 1)
        power *= z * z
    partial = 2 * total
    remainder = (
        2 * abs(z) ** (2 * terms + 1)
        / ((2 * terms + 1) * (1 - z * z))
    )
    if z < 0:
        return partial - remainder, partial
    return partial, partial + remainder


def limiting_function_bounds(alpha: Fraction) -> tuple[Fraction, Fraction]:
    """Bound alpha*(1-2 log(alpha))-1 exactly."""
    log_lo, log_hi = log_bounds(alpha)
    # Multiplication by -2*alpha reverses the log bounds.
    lower = alpha * (1 - 2 * log_hi) - 1
    upper = alpha * (1 - 2 * log_lo) - 1
    return lower, upper


def bracket_limiting_root(digits: int = 18) -> tuple[Fraction, Fraction]:
    """Bisect the nontrivial limiting root with exact log enclosures."""
    lo = Fraction(1, 4)
    hi = Fraction(1, 3)
    if not (limiting_function_bounds(lo)[1] < 0):
        raise AssertionError("lower limiting endpoint is not negative")
    if not (limiting_function_bounds(hi)[0] > 0):
        raise AssertionError("upper limiting endpoint is not positive")
    target_width = Fraction(1, 10**digits)
    while hi - lo > target_width:
        mid = (lo + hi) / 2
        lower, upper = limiting_function_bounds(mid)
        if upper < 0:
            lo = mid
        elif lower > 0:
            hi = mid
        else:
            # The enclosure is much narrower than the requested bracket at
            # 160 terms; this branch guards against silent loss of rigor.
            raise AssertionError("log enclosure did not decide the root sign")
    return lo, hi


def decimal_interval(lo: Fraction, hi: Fraction, digits: int = 16) -> str:
    """Human-readable enclosure; the rational endpoints remain authoritative."""
    scale = 10**digits
    lower = lo.numerator * scale // lo.denominator
    upper = -((-hi.numerator * scale) // hi.denominator)
    return f"[{lower / scale:.{digits}f}, {upper / scale:.{digits}f}]"


def rational_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def verify_symbolic_identities() -> dict[str, str]:
    n, r = sp.symbols("n r", positive=True)
    A = 2 * n - 1 - 2 * (n - 1) * r
    ratio = (2 * r - 1) ** (n - 1) * A / r**n
    log_derivative = sp.factor(sp.diff(sp.log(ratio), r))
    expected = sp.factor(
        n * (2 * n * r - 2 * n + 1)
        / (r * (2 * r - 1) * (2 * n * r - 2 * n - 2 * r + 1))
    )
    if sp.simplify(log_derivative - expected) != 0:
        raise AssertionError("likelihood-ratio derivative identity changed")

    r3 = (sp.Integer(19) + sp.sqrt(21)) / 34
    finite3 = (2 * r3 - 1) ** 2 * (5 - 4 * r3) - r3**3
    if sp.simplify(finite3) != 0:
        raise AssertionError("n=3 algebraic endpoint changed")
    if critical_polynomial(2, Fraction(3, 5)) != 0:
        raise AssertionError("n=2 rational endpoint changed")

    # Conditional on S, W is uniform.  Differentiating the moving threshold
    # [p1-epsilon*(p1-p0)]/[1-epsilon*W] at zero and averaging W gives
    # p1/2-p0; the upper-tail derivative has its negative, hence the sign
    # displayed below.
    epsilon, p0, p1, W = sp.symbols("epsilon p0 p1 W")
    threshold = (p1 - epsilon * (p1 - p0)) / (1 - epsilon * W)
    averaged_motion = sp.expand(sp.diff(threshold, epsilon).subs(epsilon, 0).subs(W, sp.Rational(1, 2)))
    if sp.simplify(-averaged_motion - (p1 / 2 - p0)) != 0:
        raise AssertionError("terminal-edge derivative identity changed")

    return {
        "edge_derivative": "f_Beta(2,n-1)(p_1)*(p_1/2-p_0)",
        "finite_critical_equation": "(2r-1)^(n-1)*(2n-1-2(n-1)r)=r^n",
        "likelihood_ratio_log_derivative": str(expected),
        "n3_endpoint": "r_3=(19+sqrt(21))/34",
        "limit_equation": "alpha*(1-2*log(alpha))=1",
    }


def build_certificate() -> dict[str, object]:
    identities = verify_symbolic_identities()
    finite = []
    previous_alpha_hi: Fraction | None = None
    for n in range(2, 31):
        r_lo, r_hi = bracket_finite_root(n)
        alpha_lo, alpha_hi = r_lo**n, r_hi**n
        if not (
            critical_polynomial(n, r_lo) < 0
            and critical_polynomial(n, r_hi) > 0
        ):
            raise AssertionError("finite critical bracket lost its signs")
        # This finite table also certifies the observed strict decrease
        # through n=30; the theorem only needs existence and convergence.
        if previous_alpha_hi is not None and not alpha_hi < previous_alpha_hi:
            raise AssertionError("finite critical table stopped decreasing")
        previous_alpha_hi = alpha_hi
        finite.append(
            {
                "n": n,
                "r_lower": rational_record(r_lo),
                "r_upper": rational_record(r_hi),
                "alpha_lower": rational_record(alpha_lo),
                "alpha_upper": rational_record(alpha_hi),
                "alpha_display_interval": decimal_interval(alpha_lo, alpha_hi),
                "lower_polynomial_sign": "negative",
                "upper_polynomial_sign": "positive",
            }
        )

    limit_lo, limit_hi = bracket_limiting_root()
    lower_sign = limiting_function_bounds(limit_lo)
    upper_sign = limiting_function_bounds(limit_hi)
    if not (lower_sign[1] < 0 < upper_sign[0]):
        raise AssertionError("limiting critical bracket lost its signs")

    return {
        "schema_version": 1,
        "mathematical_status": "ordinary theorem with exact rational root brackets",
        "scope": {
            "sample_sizes": "all integers n>=2",
            "comparison": "pointwise binomial Stringer--Gaffke domination",
            "nonclaim": "failure of pointwise domination does not imply Stringer undercoverage",
        },
        "symbolic_identities": identities,
        "finite_root_brackets_n2_through_n30": finite,
        "limiting_root": {
            "lower": rational_record(limit_lo),
            "upper": rational_record(limit_hi),
            "display_interval": decimal_interval(limit_lo, limit_hi),
            "lower_function_sign": "negative",
            "upper_function_sign": "positive",
        },
        "consequence": (
            "For alpha>alpha_n the terminal ordered-knot edge rises above "
            "the calibrated step cap. For every alpha above the limiting "
            "root, this occurs for all sufficiently large n."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        "terminal edge bifurcation: exact symbolic identities; "
        "29 finite root brackets; one exact limiting-root bracket"
    )


if __name__ == "__main__":
    main()
