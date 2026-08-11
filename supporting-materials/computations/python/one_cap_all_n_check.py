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
    ("exp(1-e)", math.exp(1 - math.e)),
    ("0.10", 0.10),
    ("0.05", 0.05),
    ("0.01", 0.01),
)


def exact_constant_checks() -> dict[str, str]:
    """Verify the rational comparisons used to avoid decimal arguments."""

    # The positive exponential series through 1/4! gives 65/24 < e.  A
    # geometric bound on the remaining terms also gives e < 3.  These two
    # rational bounds are all that the boundary proof needs.
    e_series_lower = sum(
        Fraction(1, math.factorial(k)) for k in range(5)
    )
    assert e_series_lower == Fraction(65, 24)
    e_series_upper = e_series_lower + Fraction(1, 100)
    assert e_series_upper < 3
    assert Fraction(8, 3) ** 3 > 16

    # At r=2 the written proof substitutes x=41/24 into a rational lower
    # bound.  This is positive exactly, so no decimal comparison is hidden.
    x = Fraction(41, 24)
    z = x / 2
    exp_complement_lower = (
        z - z**2 / 2 + z**3 / 6 - z**4 / 24
    )
    y = 1 + 2 * exp_complement_lower
    atanh_argument = (y - 1) / (y + 1)
    log_lower = 2 * (
        atanh_argument + atanh_argument**3 / 3
    )
    u = x / 3
    bernoulli_ratio_upper = 1 / (1 + u / 2 + u**2 / 6)
    r_two_margin = log_lower - bernoulli_ratio_upper
    assert r_two_margin > Fraction(1, 100)

    # For r>=3 the endpoint t=1/3 is handled after writing x=e-1.  The
    # lower bound is increasing in e and is already positive at 65/24.
    e0 = Fraction(65, 24)
    endpoint_argument = (
        54 * (e0 - 1) / (e0**2 + 70 * e0 + 37)
    )
    endpoint_margin = (
        2 * (endpoint_argument + endpoint_argument**3 / 3)
        - 8 / (e0 + 7)
    )
    assert endpoint_margin == Fraction(
        3187019776783200, 569923468569077849
    )
    assert endpoint_margin > Fraction(1, 200)

    return {
        "e_lower": "65/24 < e (positive exponential-series remainder)",
        "e_upper": "e < 65/24 + 1/100 < 3 (geometric tail bound)",
        "r_2_margin": str(r_two_margin),
        "r_ge_3_endpoint_margin": str(endpoint_margin),
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

    y = sp.symbols("y", positive=True)
    log_lower_gap_derivative = sp.factor(sp.diff(
        sp.log(y) - 2 * (y - 1) / (y + 1), y
    ))
    assert sp.simplify(
        log_lower_gap_derivative
        - (y - 1) ** 2 / (y * (y + 1) ** 2)
    ) == 0

    # A [2/2] Pade bound supplies
    #   1-exp(-z) >= z/(1+z/2+z^2/12).
    # The displayed derivative proves that the rational function on the
    # right of exp(-z) is an upper bound.
    z = sp.symbols("z", positive=True)
    pade_log_gap = (
        sp.log((z**2 + 6 * z + 12) / (z**2 - 6 * z + 12)) - z
    )
    pade_derivative = sp.factor(sp.diff(pade_log_gap, z))
    assert sp.simplify(
        pade_derivative
        + z**4 / ((z**2 - 6 * z + 12) * (z**2 + 6 * z + 12))
    ) == 0

    # For r>=3 put t=1/r and x=e-1.  After the Pade and Bernoulli-ratio
    # bounds, the derivative of the remaining comparison has the sign of
    # this quartic P(t).  Its derivative is negative on 0<t<=1/3 when
    # 0<x<2, so the comparison has no interior minimum.
    t = sp.symbols("t", positive=True)
    denominator = 1 + x * t / 2 + x**2 * t**2 / 12
    comparison = (
        sp.log(1 + x / denominator)
        - 1 / (1 + x * t / (2 * (1 + t)))
    )
    comparison_derivative = sp.factor(sp.diff(comparison, t))
    quartic = (
        t**4 * x**4
        - 12 * t**3 * x**4
        - 36 * t**3 * x**3
        - 48 * t**3 * x**2
        - 72 * t**2 * x**3
        - 180 * t**2 * x**2
        - 144 * t**2 * x
        - 120 * t * x**2
        - 144 * t * x
        + 144
    )
    expected_comparison_derivative = (
        2 * x * quartic
        / (
            (t * x + 2 * t + 2) ** 2
            * (t**2 * x**2 + 6 * t * x + 12)
            * (t**2 * x**2 + 6 * t * x + 12 * x + 12)
        )
    )
    assert sp.simplify(
        comparison_derivative - expected_comparison_derivative
    ) == 0

    quartic_derivative = sp.Poly(sp.diff(quartic, t), t)
    assert quartic_derivative.degree() == 3
    decreasing_decomposition = (
        -4 * t**2 * x**4 * (9 - t)
        - 108 * t**2 * x**3
        - 144 * t**2 * x**2
        - 144 * t * x**3
        - 360 * t * x**2
        - 288 * t * x
        - 120 * x**2
        - 144 * x
    )
    assert sp.simplify(
        sp.diff(quartic, t) - decreasing_decomposition
    ) == 0

    # The endpoint transformation used in the exact rational check above.
    e = sp.symbols("e", positive=True)
    endpoint_z = sp.factor(
        ((1 + (e - 1) / (1 + (e - 1) / 6 + (e - 1)**2 / 108)) - 1)
        / ((1 + (e - 1) / (1 + (e - 1) / 6 + (e - 1)**2 / 108)) + 1)
    )
    assert sp.simplify(
        endpoint_z - 54 * (e - 1) / (e**2 + 70 * e + 37)
    ) == 0
    endpoint_z_derivative = sp.factor(sp.diff(endpoint_z, e))
    assert sp.simplify(
        endpoint_z_derivative
        - 54 * (-e**2 + 2 * e + 107) / (e**2 + 70 * e + 37)**2
    ) == 0

    return {
        "mean_derivative": str(derivative),
        "r_one_gap": str(r_one_gap),
        "bernoulli_ratio_derivative": str(ratio_derivative),
        "boundary_tail": "(b A)^r (r+1-r b A)",
        "log_lower_gap_derivative": str(log_lower_gap_derivative),
        "pade_derivative": str(pade_derivative),
        "boundary_comparison_derivative": str(comparison_derivative),
        "boundary_quartic": str(quartic),
        "boundary_quartic_derivative": str(decreasing_decomposition),
        "endpoint_atanh_argument": str(endpoint_z),
        "endpoint_atanh_derivative": str(endpoint_z_derivative),
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
