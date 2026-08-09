"""Certify binomial screening candidates in exact rational arithmetic.

``search_two_value.py`` screens in float64 and writes any parameter point
whose screened coverage falls below the nominal level to a JSON file.  This
program is the only place a claim is allowed to come from:

* the taint values and probabilities are converted to exact rationals
  (``Fraction(float)`` is exact -- every float64 is a dyadic rational);
* the coverage is recomputed as an exact ``Fraction`` by full enumeration
  of multinomial count vectors (``coverage.coverage_exact``);
* every binomial confidence factor is enclosed between dyadic rationals;
  the binomial-CDF signs at both endpoints are evaluated exactly with
  integer arithmetic;
* those rational factor intervals are propagated through every Stringer
  comparison, and a candidate is rejected if an interval overlaps the exact
  rational mean;
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

from coverage import coverage_exact
from stringer import EXACT_FACTOR_BITS


# Exact n=400 multinomial probabilities can have more than 4,300 decimal
# digits.  Python 3.11's defensive integer-to-string limit is inappropriate
# for a program whose requested output is precisely those audited integers.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _fraction_text(value: Fraction) -> str:
    return "%d/%d" % (value.numerator, value.denominator)


def certify(cand: dict, factor_bits: int = EXACT_FACTOR_BITS) -> dict:
    n = cand["n"]
    alpha = cand["alpha"]
    method = cand.get("method", "binomial")
    if "values" in cand:
        values = [Fraction(v) for v in cand["values"]]
        probs = [Fraction(q) for q in cand["probs"]]
    else:
        values = [Fraction(cand["v1"]), Fraction(cand["v2"])]
        probs = [Fraction(cand["q1"]), Fraction(cand["q2"])]

    if method != "binomial":
        raise ValueError("formal certification is implemented only for "
                         "binomial factors")

    cov, theta, min_gap = coverage_exact(
        values, probs, n, alpha, method, factor_bits)
    nominal = 1 - Fraction(str(alpha))
    margin_ok = min_gap is not None and min_gap >= 0
    confirmed = cov < nominal and margin_ok
    factor_width = Fraction(1, 1 << factor_bits)
    q0 = Fraction(1) - sum(probs)
    return {
        "certificate_version": 2,
        "input": cand,
        "exact_distribution": {
            "positive_taint_values": [_fraction_text(v) for v in values],
            "probabilities": [_fraction_text(q) for q in probs],
            "zero_taint_probability": _fraction_text(q0),
        },
        "exact_coverage": _fraction_text(cov),
        "exact_coverage_float": float(cov),
        "exact_nominal": _fraction_text(nominal),
        "nominal": float(nominal),
        "shortfall": float(nominal - cov),
        "theta": _fraction_text(theta),
        "factor_interval_method": (
            "dyadic brackets with exact integer binomial-CDF signs"),
        "factor_bracket_bits": factor_bits,
        "max_factor_bracket_width": _fraction_text(factor_width),
        "max_factor_bracket_width_float": float(factor_width),
        "min_certified_comparison_gap": (
            _fraction_text(min_gap) if min_gap is not None else None),
        "min_certified_comparison_gap_float": (
            float(min_gap) if min_gap is not None else None),
        "margin_certified": bool(margin_ok),
        "confirmed_counterexample": bool(confirmed),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates", help="JSON file from search_two_value.py")
    ap.add_argument("--factor-bits", type=int, default=EXACT_FACTOR_BITS)
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
        r = certify(cand, args.factor_bits)
        results.append(r)
        tag = "CONFIRMED" if r["confirmed_counterexample"] else "rejected "
        print("[%s] n=%-4d %s alpha=%s coverage=%.10f "
              "shortfall=%+.3e certified-gap=%.2e "
              "factor-width<=%.2e:%s"
              % (tag, cand["n"], cand.get("method", "binomial"),
                 cand["alpha"], r["exact_coverage_float"], r["shortfall"],
                 r["min_certified_comparison_gap_float"] or 0.0,
                 r["max_factor_bracket_width_float"],
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
