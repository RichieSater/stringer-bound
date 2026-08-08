"""Recertify screening candidates in exact rational arithmetic.

``search_two_value.py`` screens in float64 and writes any parameter point
whose screened coverage falls below the nominal level to a JSON file.  This
program is the only place a claim is allowed to come from:

* the taint values and probabilities are converted to exact rationals
  (``Fraction(float)`` is exact -- every float64 is a dyadic rational);
* the coverage is recomputed as an exact ``Fraction`` by full enumeration
  of multinomial count vectors (``coverage.coverage_exact``);
* every ``SB >= theta`` comparison is certified: the smallest observed
  ``|SB - theta|`` must exceed a rigorous bound on the numerical error of
  ``SB``, namely ``(sum of factor bracket widths) * max taint <= (n+1) w``
  with ``w`` the largest bisection bracket width (``~1e-40`` at default
  precision), so no comparison can have been decided inside the noise;
* the verdict compares the exact rational coverage to the exact rational
  nominal level.

A candidate is CONFIRMED only if exact coverage < 1 - alpha and the margin
certificate holds.  Anything else is reported and dropped.

Usage:
    python3 certify.py candidates.json
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction

from mpmath import mp, mpf

from coverage import coverage_exact
from stringer import max_factor_width


def certify(cand: dict, dps: int = 50) -> dict:
    n = cand["n"]
    alpha = cand["alpha"]
    method = cand.get("method", "binomial")
    if "values" in cand:
        values = [Fraction(v) for v in cand["values"]]
        probs = [Fraction(q) for q in cand["probs"]]
    else:
        values = [Fraction(cand["v1"]), Fraction(cand["v2"])]
        probs = [Fraction(cand["q1"]), Fraction(cand["q2"])]

    cov, theta, min_margin = coverage_exact(values, probs, n, alpha,
                                            method, dps)
    nominal = 1 - Fraction(alpha)
    with mp.workdps(dps):
        err = (n + 1) * max_factor_width(n, alpha, method, dps)
        margin_ok = min_margin is not None and min_margin > err

    confirmed = cov < nominal and margin_ok
    return {
        "input": cand,
        "exact_coverage": "%d/%d" % (cov.numerator, cov.denominator),
        "exact_coverage_float": float(cov),
        "nominal": float(nominal),
        "shortfall": float(nominal - cov),
        "theta": "%d/%d" % (theta.numerator, theta.denominator),
        "min_margin": float(min_margin) if min_margin is not None else None,
        "margin_error_bound": float(err),
        "margin_certified": bool(margin_ok),
        "confirmed_counterexample": bool(confirmed),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates", help="JSON file from search_two_value.py")
    ap.add_argument("--dps", type=int, default=50)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    with open(args.candidates) as fh:
        data = json.load(fh)
    cands = []
    for rec in data:
        cands.extend(rec.get("below", []))
    if not cands:
        print("no candidates below nominal in %s; nothing to certify"
              % args.candidates)
        return 0

    results, confirmed = [], 0
    for cand in cands:
        r = certify(cand, args.dps)
        results.append(r)
        tag = "CONFIRMED" if r["confirmed_counterexample"] else "rejected "
        print("[%s] n=%-4d %s alpha=%s cov=%s (%.10f, shortfall %+.3e) "
              "margin=%.2e>err=%.2e:%s"
              % (tag, cand["n"], cand.get("method", "binomial"),
                 cand["alpha"], r["exact_coverage"],
                 r["exact_coverage_float"], -r["shortfall"],
                 r["min_margin"] or 0.0, r["margin_error_bound"],
                 r["margin_certified"]))
        confirmed += r["confirmed_counterexample"]

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(results, fh, indent=1)
    print("\n%d/%d candidates confirmed as exact counterexamples"
          % (confirmed, len(results)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
