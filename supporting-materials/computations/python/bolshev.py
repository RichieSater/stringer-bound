"""Bimpeh's coverage lower bound for the Stringer bound -- and its gap.

Bimpeh (PhD thesis, DCU 2008, ch. 5) writes the Stringer bound as the
integral of one minus the stepwise lower Clopper--Pearson confidence band
built at the order statistics, and derives (his (5.16))

    P(SB >= theta)  >=  Pbar_n := P( U_{i:n} <= p_n(i) for all i ),

with ``U_{i:n}`` uniform order statistics and ``p_n(i)`` the upper
Clopper--Pearson limits (``p_n(n) = 1``).  ``Pbar_n`` is computed exactly
by Bolshev's recursion (Shorack--Wellner, Section 9.1): with boundaries
``a_1 <= ... <= a_m``,

    Q_0 = 1,   Q_m = 1 - sum_{i=1}^{m} C(m, i) (1 - a_{m-i+1})^i Q_{m-i},

where ``Q_m`` uses the first ``m`` boundaries of the fixed sequence.

Status of the containment (established in this repository; see
``audit/BIMPEH-GAP.md``): (5.16) fails for continuous F (hand
counterexamples at n = 1, 2, corroborated by
``bimpeh_continuous_check.py``) and for atomic F (exact machine
counterexample at n = 5: coverage 31/32 = 0.96875 < Pbar_5 = 0.98746,
reproduced by ``--crosscheck``).  The derivation substitutes the band's
left limit ``q_n(i-1)`` for its value ``q_n(i)`` at the order
statistics, dropping in particular the top constraint
``F(t_{n:n}) >= alpha^{1/n}``; the corrected containment probability is
at most ``P(U_{n:n} >= alpha^{1/n}) = 1 - alpha`` for every continuous
F.  Thus the corrected event does not provide the claimed numerical
cushion above nominal coverage.

``Pbar_n`` itself is computed correctly here (Bimpeh's Table 5.1
reproduces to all printed digits); it is simply not a coverage bound.
This module exists to document that finding.

Usage:
    python3 bolshev.py --alpha 0.05 --n-max 25 --crosscheck
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction

from mpmath import mp, mpf, binomial

from coverage import coverage_exact
from stringer import factor_values


def pbar(n: int, alpha: str, dps: int = 50):
    """Bolshev recursion for P(U_{i:n} <= p_n(i), all i)."""
    with mp.workdps(dps):
        a = factor_values(n, alpha, "binomial", dps)[1:]  # a_i = p_n(i), i>=1
        q = [mpf(1)]
        for m in range(1, n + 1):
            s = mpf(0)
            for i in range(1, m + 1):
                s += binomial(m, i) * (1 - a[m - i]) ** i * q[m - i]
            q.append(1 - s)
        return q[n]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="0.05")
    ap.add_argument("--n-max", type=int, default=25)
    ap.add_argument("--crosscheck", action="store_true",
                    help="compare Pbar_n with exact atomic coverage")
    args = ap.parse_args(argv)

    nominal = 1 - mpf(args.alpha)
    frontier = 0
    print("n    Pbar_n           Pbar_n >= 1-alpha?")
    for n in range(1, args.n_max + 1):
        p = pbar(n, args.alpha)
        ok = p >= nominal
        if ok:
            frontier = n
        print("%-4d %.12f   %s" % (n, float(p), "yes" if ok else "no"))
    print("\nPbar_n >= 1-alpha for n <= %d at alpha = %s, matching Bimpeh's "
          "Table 5.1 (n <= 11 at 0.05). NOTE: Pbar_n is not a lower bound "
          "on Stringer coverage; see the reassessment of (5.16) in "
          "audit/BIMPEH-GAP.md."
          % (frontier, args.alpha))

    if args.crosscheck:
        print("\ncomparison with exact coverage for selected atomic "
              "distributions:")
        for n in (5, 8, 10):
            p = pbar(n, args.alpha)
            for v, q in (((Fraction(3, 4), Fraction(1, 4)),
                          (Fraction(1, 5), Fraction(2, 5))),
                         ((Fraction(1, 1), Fraction(1, 10)),
                          (Fraction(1, 2), Fraction(1, 4)))):
                cov, _t, _m = coverage_exact(list(v), list(q), n, args.alpha)
                holds = mpf(cov.numerator) / cov.denominator >= p
                tag = ("proposed lower-bound inequality holds" if holds
                       else "proposed lower-bound inequality fails")
                print("  n=%-3d v=%s q=%s cov=%.10f vs Pbar=%.10f "
                      "[%s]"
                      % (n, tuple(map(str, v)), tuple(map(str, q)),
                         float(cov), float(p), tag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
