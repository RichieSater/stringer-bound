"""Checks supporting the analytic all-n one-cap comparison.

The proof is written in ``theory/ONE-CAP-COMPARISON.md``.  It reduces the
Clopper--Pearson vertex inequality to a binomial upper-tail comparison and
then uses Anderson--Samuels fixed-mean monotonicity, Pinelis's binomial-mean
bound, and an elementary boundary estimate.  This program checks the
algebraic identities and exact rational constant comparisons used in that
argument.  It also runs a floating-point regression grid, which is a bug
detector rather than part of the proof.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction

import numpy as np
import sympy as sp
from scipy.special import betainc


DEFAULT_N_MAX = 500
REGRESSION_ALPHAS = (
    ("0.25", 0.25),
    ("0.10", 0.10),
    ("0.05", 0.05),
    ("0.01", 0.01),
)


def exact_constant_checks() -> dict[str, object]:
    """Verify the rational comparisons used to avoid decimal arguments."""

    # The geometric tail estimate after 1/4! gives e<3<4, hence log(4)>1.
    e_partial = sum(Fraction(1, math.factorial(k)) for k in range(5))
    assert e_partial == Fraction(65, 24)
    e_tail_upper = Fraction(1, math.factorial(5)) / (
        1 - Fraction(1, 6)
    )
    assert e_tail_upper == Fraction(1, 100)
    e_upper = e_partial + e_tail_upper
    assert e_upper < 3 < 4

    # At x=1, the common rational lower bound for the boundary multiplier
    # is already greater than one for r=2,3,4.  For r>=5 a still simpler
    # lower bound is minimized at r=5.
    finite_margins = {}
    expected = {
        2: Fraction(17, 3025),
        3: Fraction(885001, 16693124),
        4: Fraction(1841131309, 25297477335),
    }
    for r in range(2, 5):
        t = Fraction(1, r)
        d = t**2 / (1 + t)
        delta_upper = (
            (d / (1 - d / 2)) / (t + t**2 / 2)
        )
        a_lower = 1 - delta_upper
        bracket_lower = Fraction(3 * r, r + 1)
        margin = a_lower**r * bracket_lower - 1
        assert margin == expected[r]
        assert margin > 0
        finite_margins[str(r)] = str(margin)

    large_r_margin = 3 * Fraction(5, 6) ** 6 - 1
    assert large_r_margin == Fraction(73, 15552)
    assert large_r_margin > 0

    return {
        "log_four": "e < 65/24 + 1/100 < 3 < 4, hence log(4)>1",
        "finite_boundary_margins": finite_margins,
        "r_ge_5_margin": str(large_r_margin),
    }


def symbolic_identity_checks() -> dict[str, str]:
    """Check the differentiations and boundary-tail identities."""

    x, z, r, b, A = sp.symbols("x z r b A", positive=True)
    increasing_mean = z * (1 - sp.exp(-x / z))
    derivative = sp.simplify(sp.diff(increasing_mean, z))
    assert sp.simplify(
        derivative - (1 - (1 + x / z) * sp.exp(-x / z))
    ) == 0

    # The r=1 estimate uses alpha=s^2 and alpha<=1/4, hence s<=1/2.
    s = sp.symbols("s", positive=True)
    r_one_gap = sp.factor(
        2 / (1 + s) - 1 / (1 - s**2)
    )
    assert sp.simplify(
        r_one_gap - (1 - 2 * s) / (1 - s**2)
    ) == 0

    # psi(z)=z/(exp(z)-1) is decreasing but has derivative greater than -1.
    # These two facts prove that A(x) increases and b(x)A(x) decreases.
    u = sp.symbols("u", positive=True)
    bernoulli_ratio = u / (sp.exp(u) - 1)
    ratio_derivative = sp.factor(sp.diff(bernoulli_ratio, u))
    expected_derivative = (
        sp.exp(u) * (1 - u) - 1
    ) / (sp.exp(u) - 1) ** 2
    assert sp.simplify(ratio_derivative - expected_derivative) == 0
    ratio_derivative_plus_one = sp.factor(ratio_derivative + 1)
    assert sp.simplify(
        ratio_derivative_plus_one
        - sp.exp(u) * (sp.exp(u) - 1 - u)
        / (sp.exp(u) - 1) ** 2
    ) == 0

    # P(Bin(r+1,bA)>=r) has only its r and r+1 mass terms.
    boundary_tail = (
        (r + 1) * (b * A) ** r * (1 - b * A)
        + (b * A) ** (r + 1)
    )
    assert sp.simplify(
        boundary_tail - b**r * A**r * (r + 1 - r * b * A)
    ) == 0

    # At x=1, write t=1/r and d=t^2/(1+t).  The exact relation
    # 1-A=(exp(d)-1)/(exp(t)-1) and elementary exponential bounds give
    # the rational lower bounds used in the paper.
    t = sp.symbols("t", positive=True)
    d = t**2 / (1 + t)
    delta_upper = sp.factor(
        (d / (1 - d / 2)) / (t + t**2 / 2)
    )
    a_lower = sp.factor(1 - delta_upper)
    a_lower_in_r = sp.factor(a_lower.subs(t, 1 / r))
    expected_a_lower = (
        (2 * r - 1) * (2 * r**2 + 2 * r + 1)
        / ((2 * r + 1) * (2 * r**2 + 2 * r - 1))
    )
    assert sp.simplify(a_lower_in_r - expected_a_lower) == 0

    bracket_lower = sp.factor(
        1 + (1 - t / 2)
        + (1 - t) * (1 - t / 2) / (1 + t)
    )
    assert sp.simplify(bracket_lower - 3 / (1 + t)) == 0

    return {
        "mean_derivative": str(derivative),
        "r_one_gap": str(r_one_gap),
        "bernoulli_ratio_derivative": str(ratio_derivative),
        "bernoulli_ratio_derivative_plus_one": str(
            ratio_derivative_plus_one
        ),
        "boundary_tail": "(b A)^r (r+1-r b A)",
        "finite_r_a_lower": str(a_lower_in_r),
        "boundary_bracket_lower": str(bracket_lower),
    }


def numerical_regression(n_max: int = DEFAULT_N_MAX) -> list[dict[str, object]]:
    """Stress the proved tail inequality on a finite grid.

    No floating-point result from this routine is used in the proof.
    ``betainc(r,n-r+1,q)`` is the binomial probability ``P(X>=r)``.
    """

    if n_max < 2:
        raise ValueError("n_max must be at least 2")

    out: list[dict[str, object]] = []
    for alpha_text, alpha in REGRESSION_ALPHAS:
        x = -math.log(alpha)
        minimum = (math.inf, 0, 0)
        minimum_boundary_log = (math.inf, 0)

        for n in range(2, n_max + 1):
            rs = np.arange(1, n, dtype=float)
            q = -math.expm1(-x / n) / np.expm1(x / rs)
            tails = betainc(rs, n - rs + 1, q)
            ratios = tails / alpha
            index = int(np.argmin(ratios))
            candidate = (float(ratios[index]), n, index + 1)
            if candidate[0] < minimum[0]:
                minimum = candidate

        rs = np.arange(2, max(3, n_max + 1), dtype=float)
        us = x / (rs + 1)
        vs = x / rs
        log_a = np.log(-np.expm1(-us)) - np.log(-np.expm1(-vs))
        q_boundary = np.exp(-vs + log_a)
        boundary_logs = (
            rs * log_a + np.log(rs + 1 - rs * q_boundary)
        )
        index = int(np.argmin(boundary_logs))
        minimum_boundary_log = (float(boundary_logs[index]), index + 2)

        if minimum[0] < 1 - 2e-12:
            raise AssertionError(
                f"tail regression failed at alpha={alpha_text}: {minimum}"
            )
        if minimum_boundary_log[0] <= 0:
            raise AssertionError(
                "boundary multiplier regression failed at "
                f"alpha={alpha_text}: {minimum_boundary_log}"
            )

        out.append({
            "alpha": alpha_text,
            "n_max": n_max,
            "minimum_tail_divided_by_alpha": minimum[0],
            "minimum_tail_case": {"n": minimum[1], "r": minimum[2]},
            "minimum_boundary_log_multiplier": minimum_boundary_log[0],
            "minimum_boundary_r": minimum_boundary_log[1],
        })

    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-max", type=int, default=DEFAULT_N_MAX)
    args = parser.parse_args(argv)

    constants = exact_constant_checks()
    identities = symbolic_identity_checks()
    regression = numerical_regression(args.n_max)

    print("exact rational constant checks: passed")
    for key, value in constants.items():
        print(f"  {key}: {value}")
    print("symbolic identity checks: passed")
    for key, value in identities.items():
        print(f"  {key}: {value}")
    print("finite numerical regression (not part of proof): passed")
    for row in regression:
        case = row["minimum_tail_case"]
        print(
            "  alpha=%s n<=%d min tail/alpha=%.15f at (n,r)=(%d,%d); "
            "min boundary log multiplier=%.15f at r=%d"
            % (
                row["alpha"],
                row["n_max"],
                row["minimum_tail_divided_by_alpha"],
                case["n"],
                case["r"],
                row["minimum_boundary_log_multiplier"],
                row["minimum_boundary_r"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
