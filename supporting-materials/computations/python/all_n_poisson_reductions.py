"""Checks for the exact reductions in ALL-N-POISSON-PROGRAM.md.

This script does *not* prove either open all-n quantile inequality.  It checks
the algebra that isolates them:

* the ordered-weight Poisson Stringer identity;
* the exact two-exponential quantile-convexity threshold ``4*exp(-3)``;
* the gamma-kernel antiderivative identity; and
* a decisive numerical counterexample to the tempting SymPol shortcut.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from mpmath import mp
from sympy import (Rational, diff, exp, factorial, log, simplify, symbols,
                   series)

from stringer import factor_prefix


def check_ordered_weight_identity() -> None:
    """Symbolically compare the two forms of Poisson Stringer."""
    for n in range(1, 9):
        p = symbols(f"p0:{n + 1}")
        t = symbols(f"t1:{n + 1}")
        direct = p[0] + sum(
            (p[j] - p[j - 1]) * t[j - 1]
            for j in range(1, n + 1)
        )
        ordered = ((1 - t[0]) * p[0]
                   + sum((t[j - 1] - t[j]) * p[j]
                         for j in range(1, n))
                   + t[n - 1] * p[n])
        assert simplify(direct - ordered) == 0


def check_two_exponential_obstruction() -> None:
    """Verify q'(1) and q''(1) by an implicit tail expansion."""
    h, x = symbols("h x", positive=True)
    t = 1 + h
    # Upper tail of E_1+tE_2, continued analytically through t=1.
    tail = (exp(-x) - t * exp(-x / t)) / (1 - t)
    expansion = series(tail, h, 0, 3).removeO()
    expected = (exp(-x) * (1 + x)
                + h * x**2 * exp(-x) / 2
                + h**2 * exp(-x) * (x**3 / 6 - x**2 / 2))
    assert simplify(expansion - expected) == 0

    q1 = x / 2
    q2 = x * (x - 3) / 12
    q_path = x + q1 * h + q2 * h**2 / 2
    implicit = series(
        expansion.subs(x, q_path), h, 0, 3).removeO()
    baseline = exp(-x) * (1 + x)
    assert simplify(implicit - baseline) == 0


def check_two_exponential_global_convexity() -> None:
    """Verify the exact identities behind the global two-weight theorem.

    The written proof in ``TWO-EXPONENTIAL-QUANTILE.md`` supplies the sign
    arguments.  This routine independently checks the nontrivial symbolic
    factorization and the coefficient formula used there.
    """

    z = symbols("z", positive=True)
    ez = exp(z)
    b = (1 - exp(-z)) / z
    a = 1 - exp(-z) - z**2 / (ez - 1)
    h = simplify(z * diff(a, z) / a)
    expected_h = (
        z * (z * ez - ez + 1)**2
        / ((ez - 1) * ((ez - 1)**2 - z**2 * ez))
    )
    assert simplify(h - expected_h) == 0

    margin = 3 - h + log(1 + h * b) - log(4)
    k = (2 * z**2 * ez + z * ez**2 - z
         - 4 * ez**2 + 8 * ez - 4)
    expected_derivative = (
        (ez - 1 - z)**2 * (z * ez - ez + 1)**2 * k * ez
        / (
            (ez - 1)**2
            * (ez**2 - 1 - 2 * z * ez)
            * (z**2 * ez - ez**2 + 2 * ez - 1)**2
        )
    )
    assert simplify(diff(margin, z) - expected_derivative) == 0

    # K(z) has zero coefficients through degree five.  For m>=2 its
    # z^m coefficient, multiplied by m!, is the following integer.  It is
    # positive at m=6,7 and manifestly positive for every m>=8.
    def scaled_k_coefficient(m: int) -> int:
        if m == 0 or m == 1:
            return 0
        return 2 * m * (m - 1) + (m - 8) * 2 ** (m - 1) + 8

    k_series = series(k, z, 0, 25).removeO().expand()
    assert k_series.coeff(z, 0) == 0
    assert k_series.coeff(z, 1) == 0
    for m in range(2, 25):
        assert simplify(
            factorial(m) * k_series.coeff(z, m)
            - scaled_k_coefficient(m)
        ) == 0
    assert [scaled_k_coefficient(m) for m in range(6)] == [0] * 6
    assert scaled_k_coefficient(6) == 4
    assert scaled_k_coefficient(7) == 28
    assert all(scaled_k_coefficient(m) > 0 for m in range(8, 65))

    # These series also check the removable z=0 endpoint used in the proof:
    # h(0)=3 and the comparison margin starts positively at order two.
    assert simplify(series(h, z, 0, 3).removeO()
                    - (3 - z / 2 - z**2 / 60)) == 0
    assert simplify(series(margin, z, 0, 3).removeO()
                    - 3 * z**2 / 40) == 0


def check_kernel_identity() -> None:
    """Verify the nth-derivative identity for several symbolic n."""
    y = symbols("y", positive=True)
    for n in range(1, 11):
        left = diff(y**n * exp(-Rational(n, 1) / y), y, n) / factorial(n)
        right = (exp(-Rational(n, 1) / y)
                 * sum((Rational(n, 1) / y)**k / factorial(k)
                       for k in range(n + 1)))
        assert simplify(left - right) == 0
        # For y>1, the nth derivative of (y-1)^n/n! is one; for
        # 0<y<1, the positive-part term is identically zero.
        assert diff((y - 1)**n, y, n) / factorial(n) == 1


def sympol_shortcut_counterexample() -> tuple[mp.mpf, mp.mpf, mp.mpf, int]:
    """Return the n=15, alpha=.05 one-full-taint comparison.

    The calculation is deliberately labeled numerical.  Its margin is over
    2e-3, so it is more than enough to rule out this proposed proof route.
    """
    n = 15
    errors = 1
    alpha = Fraction(1, 20)
    with mp.workdps(80):
        poisson_p1 = factor_prefix(
            n, str(float(alpha)), errors, "poisson", 80)[errors][0]
        candidates = []
        nonzero_complements = n - errors
        for k in range(1, nonzero_complements + 1):
            elementary_mean = (mp.mpf(comb(nonzero_complements, k))
                               / comb(n, k))
            candidates.append(
                (mp.mpf(alpha.numerator) / alpha.denominator
                 * elementary_mean) ** (mp.mpf(1) / k))
        index = max(range(len(candidates)), key=candidates.__getitem__)
        sympol_upper = 1 - candidates[index]
        gap = poisson_p1 - sympol_upper
        assert gap < mp.mpf("-0.002")
        return +poisson_p1, +sympol_upper, +gap, index + 1


def main() -> int:
    check_ordered_weight_identity()
    check_two_exponential_obstruction()
    check_two_exponential_global_convexity()
    check_kernel_identity()
    poisson, sympol, gap, active_k = sympol_shortcut_counterexample()
    print("all-n reduction algebra: PASS")
    print("exact two-exponential convexity threshold: 4*exp(-3) =",
          mp.nstr(4 * mp.e**-3, 16))
    print("SymPol shortcut counterexample (numerical):")
    print("  n=15 alpha=0.05 one full taint")
    print("  active elementary-symmetric order:", active_k)
    print("  Poisson Stringer:", mp.nstr(poisson, 16))
    print("  SymPol upper:", mp.nstr(sympol, 16))
    print("  Stringer - SymPol:", mp.nstr(gap, 16))
    print("higher-dimensional inequality A and general inequality B: "
          "NOT CLAIMED AS PROVED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
