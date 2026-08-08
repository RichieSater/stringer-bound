"""Screening search over taint supports with ``m`` nonzero values.

Generalises ``search_two_value.py``: the support is
``{v_1 > v_2 > ... > v_m > 0}`` plus an atom at 0.  For fixed values and
``n`` the Stringer bound of a sample depends only on the count vector
``(k_1, ..., k_m)``; all ``C(n+m, m)`` vectors are enumerated once, their
bounds precomputed, and coverage for a probability vector ``q`` is a single
vectorised multinomial sum.  Probability vectors are optimised by
Nelder--Mead from Dirichlet-style restarts on a simplex grid.

Everything printed is float64 screening; candidates below nominal go to a
JSON file for ``certify.py`` (which re-certifies in exact arithmetic via
``coverage_exact`` for any ``m``).

Usage:
    python3 search_multi_value.py --alpha 0.05 --m 3 --n 10 20 30 \
        --out candidates3.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from typing import List, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln

from coverage import float_factors


class MultiValueGrid:
    """Precomputed Stringer bounds for every count vector at fixed values."""

    def __init__(self, values: Sequence[float], n: int, p: np.ndarray):
        self.values = np.asarray(values, dtype=float)
        m = len(values)
        self.n, self.m = n, m
        counts = [c for c in itertools.product(range(n + 1), repeat=m)
                  if sum(c) <= n]
        self.k = np.array(counts, dtype=np.int64)          # (N, m)
        self.k0 = n - self.k.sum(axis=1)                   # (N,)
        # telescoped bound: SB = p_0 + sum_i v_i (p_{c_i} - p_{c_{i-1}})
        cum = np.cumsum(self.k, axis=1)                    # (N, m)
        prev = np.hstack([np.zeros((len(counts), 1), np.int64), cum[:, :-1]])
        self.sb = p[0] + (self.values * (p[cum] - p[prev])).sum(axis=1)
        self.logcoef = (gammaln(n + 1)
                        - gammaln(self.k + 1.0).sum(axis=1)
                        - gammaln(self.k0 + 1.0))

    def coverage(self, q: np.ndarray) -> float:
        q0 = 1.0 - q.sum()
        if q0 < -1e-12 or (q < 0).any():
            return 1.0
        q0 = max(q0, 0.0)
        theta = float(self.values @ q)
        logw = self.logcoef.copy()
        for i in range(self.m):
            logw += _xlogy(self.k[:, i], q[i])
        logw += _xlogy(self.k0, q0)
        w = np.exp(logw)
        return float(w[self.sb >= theta].sum())


def _xlogy(k, q):
    k = np.asarray(k, dtype=float)
    if q <= 0:
        return np.where(k > 0, -np.inf, 0.0)
    return k * np.log(q)


def simplex_grid(m: int, steps: int) -> List[np.ndarray]:
    """Interior points of the probability simplex over m nonzero atoms."""
    axis = np.linspace(0.0, 1.0, steps + 1)[1:]
    out = []
    for combo in itertools.product(axis, repeat=m):
        if sum(combo) <= 1.0:
            out.append(np.array(combo))
    return out


def value_tuples(m: int, steps: int) -> List[Tuple[float, ...]]:
    grid = sorted(set(
        [round(x, 6) for x in np.linspace(0.05, 1.0, steps)]
        + [1.0, 0.99, 0.95, 0.9, 0.75, 0.5, 0.25, 0.1, 0.05]
    ), reverse=True)
    return [t for t in itertools.combinations(grid, m)]


def refine(grid: MultiValueGrid, q0: np.ndarray) -> Tuple[np.ndarray, float]:
    def obj(x):
        if (x <= 0).any() or x.sum() >= 1.0:
            return 1.0
        return grid.coverage(x)

    res = minimize(obj, q0, method="Nelder-Mead",
                   options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 4000})
    x = res.x
    if (x <= 0).any() or x.sum() >= 1.0:
        return q0, grid.coverage(q0)
    return x, float(res.fun)


def sweep(ns, alpha: str, method: str, m: int, value_steps: int,
          prob_steps: int, keep: int):
    nominal = 1.0 - float(alpha)
    tuples = value_tuples(m, value_steps)
    results = []
    for n in ns:
        p = float_factors(n, alpha, method)
        qs = simplex_grid(m, prob_steps)
        best = (np.inf, None, None)
        below = []
        for values in tuples:
            grid = MultiValueGrid(values, n, p)
            cand = sorted((grid.coverage(q), tuple(q)) for q in qs)[:keep]
            for c, q in cand:
                rq, rc = refine(grid, np.array(q))
                if rc < best[0]:
                    best = (rc, values, tuple(float(x) for x in rq))
                if rc < nominal:
                    below.append({"n": n, "alpha": alpha, "method": method,
                                  "values": list(values),
                                  "probs": [float(x) for x in rq],
                                  "coverage": rc})
        results.append({"n": n, "alpha": alpha, "method": method, "m": m,
                        "min_coverage": best[0], "argmin_values": best[1],
                        "argmin_probs": best[2], "nominal": nominal,
                        "below": below})
        print("n=%-4d alpha=%s m=%d min coverage = %.10f (nominal %.4f, "
              "slack %+.3e) at v=%s q=%s"
              % (n, alpha, m, best[0], nominal, best[0] - nominal,
                 best[1], tuple(round(x, 5) for x in best[2])), flush=True)
    return results


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="0.05")
    ap.add_argument("--method", default="binomial",
                    choices=["binomial", "poisson"])
    ap.add_argument("--m", type=int, default=3)
    ap.add_argument("--n", nargs="+", type=int, required=True)
    ap.add_argument("--value-steps", type=int, default=8)
    ap.add_argument("--prob-steps", type=int, default=10)
    ap.add_argument("--keep", type=int, default=2)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    res = sweep(args.n, args.alpha, args.method, args.m, args.value_steps,
                args.prob_steps, args.keep)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)
    worst = min(r["min_coverage"] for r in res)
    print("\nglobal minimum screening coverage: %.10f (nominal %.4f)"
          % (worst, 1 - float(args.alpha)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
