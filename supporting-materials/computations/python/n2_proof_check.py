"""Machine checks for the n = 2 theorem (theory/N2-PROOF.md).

Four layers:
1. the reduction identity SB = 1 - (A min U + B max U) against the
   direct Stringer definition (random samples);
2. every symbolic identity used in the proof (sympy): beta^2 - B^2 =
   2 gamma (1-gamma); h(B) = B; h(beta) = 0; k(B) = k(beta) = beta;
   h'' = -gamma^2/(1-g)^3; beta - B = A; and A = B at alpha = 16/25;
3. the potential bound k(g) >= beta on [B, beta] numerically across
   alpha (concavity + equal endpoints makes this rigorous; the numeric
   scan guards against algebra slips);
4. the theorem itself, stress-tested against tens of thousands of
   atomic distributions per alpha, including the boundary families.

This script corroborates the written proof; the proof stands on its own.

Usage:
    python3 n2_proof_check.py
"""

from __future__ import annotations

import itertools
import sys

import numpy as np
import sympy as sp

from stringer import factor_values, stringer_bound


def check_reduction(trials: int = 20) -> bool:
    rng = np.random.default_rng(7)
    for alpha in ("0.05", "0.3", "0.8"):
        p = [float(x) for x in factor_values(2, alpha)]
        A, B = p[1] - p[0], 1 - p[1]
        for _ in range(trials):
            t = rng.random(2)
            sb = float(stringer_bound(list(t), 2, alpha))
            u = 1 - t
            sb2 = 1 - (A * u.min() + B * u.max())
            if abs(sb - sb2) > 1e-12:
                return False
    return True


def check_identities() -> bool:
    a = sp.symbols("alpha", positive=True)
    g = sp.symbols("g")
    beta, gamma = sp.sqrt(a), sp.sqrt(1 - a)
    A, B = beta + gamma - 1, 1 - gamma
    h = (beta ** 2 - g ** 2) / (2 * (1 - g))
    k = g + (A / B) * h
    checks = [
        sp.simplify(beta ** 2 - B ** 2 - 2 * gamma * (1 - gamma)) == 0,
        sp.simplify(h.subs(g, B) - B) == 0,
        sp.simplify(h.subs(g, beta)) == 0,
        sp.simplify(k.subs(g, B) - beta) == 0,
        sp.simplify(k.subs(g, beta) - beta) == 0,
        sp.simplify(sp.diff(h, g, 2) + gamma ** 2 / (1 - g) ** 3) == 0,
        sp.simplify(beta - B - A) == 0,
        sp.simplify((A - B).subs(a, sp.Rational(16, 25))) == 0,
        sp.simplify(sp.expand((beta + gamma) ** 2 - 1
                              - 2 * sp.sqrt(a) * sp.sqrt(1 - a))) == 0,
    ]
    return all(checks)


def check_potential(alphas=(0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.95)) -> bool:
    for a in alphas:
        b, c = a ** 0.5, (1 - a) ** 0.5
        A, B = b + c - 1, 1 - c
        gs = np.linspace(B, b, 4001)
        k = gs + (A / B) * (b ** 2 - gs ** 2) / (2 * (1 - gs))
        if k.min() < b - 1e-12:
            return False
    return True


def closed_wedge_tail_prob(atoms, probs, A, B, w):
    """Probability of the closed wedge tail, enlarged by roundoff tolerance."""
    p = 0.0
    for (i, u), (j, v) in itertools.product(enumerate(atoms), repeat=2):
        if A * min(u, v) + B * max(u, v) >= w - 1e-15:
            p += probs[i] * probs[j]
    return p


def stress(alphas=(0.05, 0.1, 0.5, 0.7, 0.9), trials: int = 40000) -> float:
    rng = np.random.default_rng(3)
    worst = -1.0
    for a in alphas:
        b, c = a ** 0.5, (1 - a) ** 0.5
        A, B = b + c - 1, 1 - c
        for trial in range(trials):
            k = rng.integers(2, 5)
            atoms = np.sort(rng.random(k))
            if trial % 4 == 0:
                atoms[0] = 0.0
            if trial % 7 == 0:
                atoms[-1] = 1.0
            pr = rng.dirichlet(np.ones(k) * rng.choice([0.3, 1.0, 3.0]))
            w = float(atoms @ pr)
            if w <= 1e-9:
                continue
            worst = max(
                worst, closed_wedge_tail_prob(atoms, pr, A, B, w) - a)
        for q in np.linspace(0.001, 0.999, 400):
            for u in np.linspace(0.01, 1.0, 50):
                worst = max(
                    worst,
                    closed_wedge_tail_prob(
                        [0.0, u], [1 - q, q], A, B, q * u) - a)
    return worst


def main():
    ok = True
    r = check_reduction()
    print("reduction identity SB = 1 - (A minU + B maxU):",
          "OK" if r else "FAIL")
    ok &= r
    r = check_identities()
    print("symbolic identities (h(B)=B, k endpoints, h'' concavity, ...):",
          "OK" if r else "FAIL")
    ok &= r
    r = check_potential()
    print("potential k(g) >= beta on [B, beta] across alpha:",
          "OK" if r else "FAIL")
    ok &= r
    excess = stress()
    print("theorem stress test: max[P(Am+BM >= w) - alpha] over atomic F "
          "= %+.2e %s"
          % (excess, "OK" if excess <= 1e-12 else "FAIL"))
    ok &= excess <= 1e-12
    print("n2_proof_check:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
