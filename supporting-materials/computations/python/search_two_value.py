"""Search for taint distributions on which the Stringer bound under-covers.

Scope of this search: taint distributions supported on ``{v1 > v2 > 0}`` plus
an atom at 0.  Distributions with a single nonzero taint value are excluded on
purpose -- ``two_point_lemma.py`` shows they can never under-cover -- so
``{v1, v2, 0}`` is the smallest support that can possibly carry a
counterexample.

The search is a coarse grid over ``(v1, v2, q1, q2)`` followed by a
Nelder--Mead refinement of the worst grid cells.  Everything it prints is a
float64 screening number; any candidate below the nominal level is written to
``candidates.json`` and must be re-certified by ``certify.py`` in exact
rational arithmetic before it is claimed.

Usage:
    python3 search_two_value.py --alpha 0.05 --n 5 60 --out candidates.json
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import List, Tuple

import numpy as np
from scipy.optimize import minimize

from coverage import TwoValueGrid, float_factors


def value_pairs(steps: int) -> List[Tuple[float, float]]:
    """Grid of (v1, v2) with 0 < v2 < v1 <= 1, refined near the extremes."""
    grid = sorted(set(
        [round(x, 6) for x in np.linspace(0.02, 1.0, steps)]
        + [1.0, 0.99, 0.95, 0.9, 0.75, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01]
    ))
    out = []
    for i, v1 in enumerate(grid):
        for v2 in grid[:i]:
            out.append((v1, v2))
    return out


def prob_grid(steps: int) -> List[Tuple[float, float]]:
    axis = np.linspace(0.0, 1.0, steps + 1)[1:]
    out = []
    for q1 in axis:
        for q2 in axis:
            if q1 + q2 <= 1.0:
                out.append((float(q1), float(q2)))
    return out


def refine(grid: TwoValueGrid, q1: float, q2: float) -> Tuple[float, float, float]:
    """Local Nelder--Mead descent on coverage in (q1, q2), values fixed."""

    def obj(x):
        a, b = x
        if a <= 0 or b <= 0 or a + b >= 1.0:
            return 1.0
        return grid.coverage(a, b)

    res = minimize(obj, np.array([q1, q2]), method="Nelder-Mead",
                   options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 2000})
    a, b = res.x
    if a <= 0 or b <= 0 or a + b >= 1.0:
        return q1, q2, grid.coverage(q1, q2)
    return float(a), float(b), float(res.fun)


def sweep(ns, alpha: str, method: str, value_steps: int, prob_steps: int,
          keep: int, verbose: bool = True):
    nominal = 1.0 - float(alpha)
    pairs = value_pairs(value_steps)
    qs = prob_grid(prob_steps)
    results = []
    for n in ns:
        p = float_factors(n, alpha, method)
        best = (np.inf, None)
        below = []
        for v1, v2 in pairs:
            grid = TwoValueGrid(v1, v2, n, p)
            cand = []
            for q1, q2 in qs:
                c = grid.coverage(q1, q2)
                cand.append((c, q1, q2))
            cand.sort()
            for c, q1, q2 in cand[:keep]:
                rq1, rq2, rc = refine(grid, q1, q2)
                if rc < best[0]:
                    best = (rc, (v1, v2, rq1, rq2))
                if rc < nominal:
                    below.append({"n": n, "alpha": alpha, "method": method,
                                  "v1": v1, "v2": v2, "q1": rq1, "q2": rq2,
                                  "coverage": rc})
        rec = {"n": n, "alpha": alpha, "method": method,
               "min_coverage": best[0], "argmin": best[1],
               "nominal": nominal, "n_below": len(below), "below": below}
        results.append(rec)
        if verbose:
            v = best[1]
            print("n=%-4d alpha=%s %-8s min coverage = %.10f  "
                  "(nominal %.4f, slack %+.3e)  at v=(%.4f, %.4f) q=(%.6f, %.6f)"
                  % (n, alpha, method, best[0], nominal, best[0] - nominal,
                     v[0], v[1], v[2], v[3]), flush=True)
    return results


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="0.05")
    ap.add_argument("--method", default="binomial", choices=["binomial", "poisson"])
    ap.add_argument("--n", nargs="+", type=int, required=True,
                    help="explicit sample sizes, or 'lo hi step' with --range")
    ap.add_argument("--range", action="store_true")
    ap.add_argument("--value-steps", type=int, default=12)
    ap.add_argument("--prob-steps", type=int, default=24)
    ap.add_argument("--keep", type=int, default=3)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    ns = (list(range(args.n[0], args.n[1] + 1, args.n[2] if len(args.n) > 2 else 1))
          if args.range else args.n)
    res = sweep(ns, args.alpha, args.method, args.value_steps,
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
