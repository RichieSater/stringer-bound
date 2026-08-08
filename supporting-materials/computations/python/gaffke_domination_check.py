"""Does the n=2 Stringer bound pointwise dominate Gaffke's bound?

Gaffke's upper confidence bound for the mean of a [0,1]-valued variable
(Gaffke 2005; finite-sample validity for independent variables proved by
Vlassis & Thomas, arXiv:2607.08415) is, for a sample of n taints, the
(1-alpha)-quantile of  sum_i w_i z_i  with z = (t_(1), ..., t_(n), 1)
ascending and w ~ Dirichlet(1, ..., 1).

If SB >= Gaffke pointwise on every sample, Stringer conservatism at that
(n, alpha) would follow from Gaffke validity.  This script computes the
exact quantile at n = 2 (the survival function of a linear form in
uniform-simplex weights is a standard piecewise-quadratic formula) and
scans the sample square.

Finding (recorded in audit/NOVELTY.md and the manuscript): at
alpha = 0.05 the domination HOLDS with equality exactly on the two
extremal families (equal taints; one taint = 1); at alpha = 0.7 it FAILS
(gap -0.032) while the n=2 theorem still holds.  So the wedge-inequality
proof is strictly stronger than the Gaffke route.

Usage:
    python3 gaffke_domination_check.py
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.optimize import brentq


def surv(cs, s):
    """P(sum c_i w_i > s), w uniform on the 2-simplex, distinct c_i."""
    tot = 0.0
    for i, ci in enumerate(cs):
        if ci > s:
            denom = 1.0
            for j, cj in enumerate(cs):
                if j != i:
                    denom *= ci - cj
            tot += (ci - s) ** 2 / denom
    return tot


def gaffke2(lo, hi, alpha):
    cs = [lo + 0e-9, hi + 1e-9, 1.0 + 2e-9]   # break ties infinitesimally
    return brentq(lambda s: surv(cs, s) - alpha, min(cs) - 1e-8, max(cs))


def scan(alpha, n_lo=101, n_hi=60):
    b, c = alpha ** 0.5, (1 - alpha) ** 0.5
    p0, p1 = 1 - b, c
    worst, arg = np.inf, None
    for lo in np.linspace(0, 1, n_lo):
        for hi in np.linspace(lo, 1, n_hi):
            sb = p0 + (p1 - p0) * hi + (1 - p1) * lo
            d = sb - gaffke2(lo, hi, alpha)
            if d < worst:
                worst, arg = d, (round(lo, 4), round(hi, 4))
    return worst, arg


def main():
    ok = True
    for alpha, expect_dom in ((0.05, True), (0.10, True), (0.7, False)):
        worst, arg = scan(alpha)
        dom = worst > -1e-6
        print("alpha=%.2f  min[SB - Gaffke] = %+.4e at %s  -> %s"
              % (alpha, worst, arg,
                 "SB dominates" if dom else "domination fails"))
        ok &= (dom == expect_dom)
    print("gaffke_domination_check:",
          "matches recorded findings" if ok else "UNEXPECTED CHANGE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
