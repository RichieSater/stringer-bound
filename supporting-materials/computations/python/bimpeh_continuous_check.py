"""Monte Carlo confirmation of the hand counterexamples to Bimpeh (5.16).

F is continuous: mass 0.99 uniform on (0.99, 1), mass 0.01 uniform on
(0, 0.01); mu = 0.9851.  Hand calculations (audit/BIMPEH-GAP.md):

    n = 1:  CP = 0.99    <  Pbar_1 = 1
    n = 2:  CP = 0.9801  <  Pbar_2 = 2 sqrt(0.95) - 0.95 = 0.999359

so Bimpeh's inequality CP >= Pbar_n fails for continuous F.  This script
only corroborates arithmetic that is already rigorous by hand; the Monte
Carlo is not part of any proof.

Usage:
    python3 bimpeh_continuous_check.py
"""

from __future__ import annotations

import sys

import numpy as np

from bolshev import pbar
from stringer import factor_values

N = 4_000_000
MU = 0.99 * 0.995 + 0.01 * 0.005
HAND = {1: 0.99, 2: 0.9801}


def draw(rng, n):
    hi = rng.random((N, n)) < 0.99
    u = rng.random((N, n))
    return np.where(hi, 0.99 + 0.01 * u, 0.01 * u)


def main():
    rng = np.random.default_rng(1)
    ok = True
    for n in (1, 2):
        p = np.array([float(x) for x in factor_values(n, "0.05")])
        t = np.sort(draw(rng, n), axis=1)[:, ::-1]
        sb = p[0] + ((p[1:n + 1] - p[0:n]) * t).sum(axis=1)
        cp = (sb >= MU).mean()
        moe = 3 * np.sqrt(cp * (1 - cp) / N)
        pb = float(pbar(n, "0.05"))
        agree = abs(cp - HAND[n]) <= moe
        refutes = HAND[n] < pb
        ok = ok and agree and refutes
        print("n=%d  MC CP = %.5f +/- %.5f  hand CP = %.5f  Pbar = %.6f  "
              "[MC agrees: %s; CP < Pbar: %s]"
              % (n, cp, moe, HAND[n], pb, agree, refutes))
    print("bimpeh_continuous_check:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
