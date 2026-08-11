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
    ("exp(-2)", math.exp(-2)),
    ("0.10", 0.10),
    ("0.05", 0.05),
    ("0.01", 0.01),
)


def exact_constant_checks() -> dict[str, str]:
    """Verify the rational comparisons used to avoid decimal arguments."""

    # The exponential series gives 8/3 < e < 87/32 < 11/4.  The first
    # bound is its partial sum through 1/3!, and the second follows by
    # bounding the remaining tail geometrically after that same partial sum.
    e_series_upper = Fraction(8, 3) + Fraction(5, 96)
    assert e_series_upper == Fraction(87, 32)
    assert e_series_upper < Fraction(11, 4)

    # The r=2 boundary case uses e^(3/4) < 9/4.
    assert Fraction(11, 4) ** 3 < Fraction(9, 4) ** 4

    return {
        "e_lower": "8/3 < e (positive exponential-series remainder)",
        "e_upper": "e < 87/32 < 11/4",
        "r_2_log": "e^(3/4) < 9/4",
    }


def symbolic_identity_checks() -> dict[str, str]:
    """Check the differentiations and the two-term boundary-tail identity."""

    x, z, r, b, A = sp.symbols("x z r b A", positive=True)
    increasing_mean = z * (1 - sp.exp(-x / z))
    derivative = sp.simplify(sp.diff(increasing_mean, z))
    assert sp.simplify(
        derivative - (1 - (1 + x / z) * sp.exp(-x / z))
    ) == 0

    # The r=1 estimate uses alpha=s^2 and alpha<=1/4, hence s<=1/2.
    s = sp.symbols("s", positive=True)
    r_one_gap = sp.factor(
        2 / (1 + s) - 1 / (1 - s ** 2)
    )
    assert sp.simplify(
        r_one_gap - (1 - 2 * s) / (1 - s ** 2)
    ) == 0

    u = sp.symbols("u", positive=True)
    bernoulli_ratio = u / (sp.exp(u) - 1)
    ratio_derivative = sp.factor(sp.diff(bernoulli_ratio, u))
    expected_derivative = (
        sp.exp(u) * (1 - u) - 1
    ) / (sp.exp(u) - 1) ** 2
    assert sp.simplify(ratio_derivative - expected_derivative) == 0

    # P(Bin(r+1,bA)>=r) has only its r and r+1 mass terms.
    boundary_tail = ((r + 1) * (b * A) ** r * (1 - b * A)
                     + (b * A) ** (r + 1))
    assert sp.simplify(
        boundary_tail - b ** r * A ** r * (r + 1 - r * b * A)
    ) == 0

    # At x=2, the r>=3 estimate compares a standard logarithm lower
    # bound with the exponential-series upper bound for u/(exp(u)-1).
    r_ge_3_gap = sp.factor(
        2 * (r - 1) / (2 * r - 1) - (r + 1) / (r + 2)
    )
    assert sp.simplify(
        r_ge_3_gap - (r - 3) / ((2 * r - 1) * (r + 2))
    ) == 0

    y = sp.symbols("y", positive=True)
    log_lower_gap_derivative = sp.factor(sp.diff(
        sp.log(y) - 2 * (y - 1) / (y + 1), y
    ))
    assert sp.simplify(
        log_lower_gap_derivative
        - (y - 1) ** 2 / (y * (y + 1) ** 2)
    ) == 0

    return {
        "mean_derivative": str(derivative),
        "r_one_gap": str(r_one_gap),
        "bernoulli_ratio_derivative": str(ratio_derivative),
        "boundary_tail": "(b A)^r (r+1-r b A)",
        "r_ge_3_gap": str(r_ge_3_gap),
        "log_lower_gap_derivative": str(log_lower_gap_derivative),
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
        boundary_logs = (
            np.log1p(rs * (-np.expm1(-x / rs)))
            - us / np.expm1(us)
        )
        index = int(np.argmin(boundary_logs))
        minimum_boundary_log = (float(boundary_logs[index]), index + 2)

        if minimum[0] < 1 - 2e-12:
            raise AssertionError(
                f"tail regression failed at alpha={alpha_text}: {minimum}"
            )
        if minimum_boundary_log[0] <= 0:
            raise AssertionError(
                "boundary lower bound regression failed at "
                f"alpha={alpha_text}: {minimum_boundary_log}"
            )

        out.append({
            "alpha": alpha_text,
            "n_max": n_max,
            "minimum_tail_divided_by_alpha": minimum[0],
            "minimum_tail_case": {"n": minimum[1], "r": minimum[2]},
            "minimum_boundary_log_lower_bound": minimum_boundary_log[0],
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
            "min H=%.15f at r=%d"
            % (
                row["alpha"],
                row["n_max"],
                row["minimum_tail_divided_by_alpha"],
                case["n"],
                case["r"],
                row["minimum_boundary_log_lower_bound"],
                row["minimum_boundary_r"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
