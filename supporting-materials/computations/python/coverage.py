"""Finite-sample coverage of the Stringer bound.

Model (the standard monetary-unit-sampling model in which the bound is stated
and in which every published conservatism claim is made):

* ``n`` dollar units are drawn i.i.d. from the population;
* the taint ``T`` of a drawn dollar unit has distribution ``F`` on ``[0, 1]``;
* the target parameter is ``theta = E[T]``;
* coverage is ``P(SB >= theta)``.

For a taint distribution with finite support the coverage is a finite sum and
can be evaluated **exactly**: the Stringer bound depends on the sample only
through the multiset of observed taints, i.e. through the multinomial count
vector, so

    coverage = sum over count vectors c of  1[SB(c) >= theta] * Multinomial(c; q).

Two evaluators are provided.

``coverage_exact``  rational arithmetic throughout the binomial certificate:
                    exact multinomial weights, dyadic factor intervals whose
                    binomial-CDF endpoint signs are evaluated with integers,
                    and rational interval propagation through ``SB``.  It
                    returns the exact coverage and the smallest certified gap
                    between a bound interval and ``theta``.

``coverage_fast``   float64/numpy, specialised to support ``{v1 > v2 > 0}``
                    plus an atom at 0.  Used for search; every candidate it
                    reports is re-checked with ``coverage_exact``.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, lcm
from typing import Sequence, Tuple

import numpy as np
from stringer import (EXACT_FACTOR_BITS, exact_binomial_factor_brackets,
                      factors)


def float_factors(n: int, alpha: str, method: str = "binomial") -> np.ndarray:
    return np.array([float(p) for p in
                     (v for v, _ in factors(n, alpha, method))], dtype=float)


# --------------------------------------------------------------------------
# exact evaluator
# --------------------------------------------------------------------------

def _count_vectors(n: int, m: int):
    """All (c_1, ..., c_m) with c_i >= 0 and sum c_i <= n."""
    if m == 0:
        yield ()
        return
    for head in range(n + 1):
        for tail in _count_vectors(n - head, m - 1):
            yield (head,) + tail


class UncertifiedComparison(RuntimeError):
    """A rational factor interval overlaps the comparison threshold."""


def _bound_interval_from_counts(
        values: Sequence[Fraction], counts: Sequence[int],
        brackets: Sequence[Tuple[Fraction, Fraction]]
) -> Tuple[Fraction, Fraction]:
    """Rational enclosure of the Stringer bound for one count vector."""
    # Only the cumulative-count indices can have nonzero coefficients.  A
    # sparse representation matters here: for three-point populations there
    # are at most three such indices, whereas a dense length-(n+1) vector in
    # every one of O(n^2) count-vector iterations is needlessly expensive.
    coefficients = {0: Fraction(1)}
    cumulative = 0
    for value, count in zip(values, counts):
        if count:
            previous = cumulative
            cumulative += count
            coefficients[cumulative] = (coefficients.get(
                cumulative, Fraction(0)) + value)
            coefficients[previous] = (coefficients.get(
                previous, Fraction(0)) - value)

    lower = Fraction(0)
    upper = Fraction(0)
    for index, coefficient in coefficients.items():
        if coefficient == 0:
            continue
        factor_lower, factor_upper = brackets[index]
        if coefficient >= 0:
            lower += coefficient * factor_lower
            upper += coefficient * factor_upper
        else:
            lower += coefficient * factor_upper
            upper += coefficient * factor_lower
    return lower, upper


def coverage_exact(values: Sequence[Fraction], probs: Sequence[Fraction],
                   n: int, alpha: str, method: str = "binomial",
                   factor_bits: int = EXACT_FACTOR_BITS
                   ) -> Tuple[Fraction, Fraction, Fraction]:
    """Exact coverage for a finitely supported taint distribution.

    ``values``  strictly decreasing positive rational taints.
    ``probs``   their probabilities; ``1 - sum(probs)`` is the atom at 0.

    Returns ``(coverage, theta, min_gap)``.  The first two values are exact
    ``Fraction`` objects.  ``min_gap`` is the smallest rational separation
    between ``theta`` and the certified Stringer-bound interval on the
    appropriate side of the comparison.

    Formal interval certification is currently implemented for binomial
    factors only.  Poisson computations remain high-precision numerical
    calculations in :mod:`stringer`.
    """
    if method != "binomial":
        raise ValueError("exact certification is implemented only for "
                         "binomial factors")
    values = [Fraction(v) for v in values]
    probs = [Fraction(p) for p in probs]
    if any(values[i] <= values[i + 1] for i in range(len(values) - 1)):
        raise ValueError("values must be strictly decreasing")
    if any(v <= 0 or v > 1 for v in values):
        raise ValueError("values must lie in (0, 1]")
    q0 = 1 - sum(probs)
    if q0 < 0 or any(p < 0 for p in probs):
        raise ValueError("probabilities must be nonnegative and sum to <= 1")

    theta = sum(v * p for v, p in zip(values, probs))
    brackets = exact_binomial_factor_brackets(n, alpha, factor_bits)

    m = len(values)

    # Put every multinomial probability over one common denominator.  Adding
    # Fraction objects inside the count-vector loop repeatedly computes gcds
    # of integers with thousands of digits.  The common-denominator form is
    # the same exact arithmetic but reduces only once, after enumeration.
    probability_denominator = 1
    for probability in [q0] + probs:
        probability_denominator = lcm(
            probability_denominator, probability.denominator)
    probability_numerators = [
        int(probability * probability_denominator)
        for probability in probs
    ]
    q0_numerator = int(q0 * probability_denominator)
    powers = [
        [numerator ** exponent for exponent in range(n + 1)]
        for numerator in probability_numerators
    ]
    q0_powers = [q0_numerator ** exponent for exponent in range(n + 1)]
    common_weight_denominator = probability_denominator ** n

    covered_numerator = 0
    min_gap = None
    for counts in _count_vectors(n, m):
        k0 = n - sum(counts)
        multinomial = 1
        rem = n
        for c in counts:
            multinomial *= comb(rem, c)
            rem -= c
        weight_numerator = multinomial * q0_powers[k0]
        for c, probability_powers in zip(counts, powers):
            weight_numerator *= probability_powers[c]
        if weight_numerator == 0:
            continue

        sb_lower, sb_upper = _bound_interval_from_counts(
            values, counts, brackets)
        if sb_lower >= theta:
            covered_numerator += weight_numerator
            gap = sb_lower - theta
        elif sb_upper < theta:
            gap = theta - sb_upper
        else:
            raise UncertifiedComparison(
                "factor interval overlaps theta for counts=%r; increase "
                "factor_bits above %d" % (counts, factor_bits))
        if min_gap is None or gap < min_gap:
            min_gap = gap
    covered = Fraction(covered_numerator, common_weight_denominator)
    return covered, theta, min_gap


# --------------------------------------------------------------------------
# fast evaluator: support {v1 > v2 > 0} u {0}
# --------------------------------------------------------------------------

class TwoValueGrid:
    """Precomputed Stringer bounds for every (k1, k2) at fixed (v1, v2, n)."""

    def __init__(self, v1: float, v2: float, n: int, p: np.ndarray):
        self.n = n
        self.v1, self.v2 = v1, v2
        k1 = np.arange(n + 1)[:, None]
        k2 = np.arange(n + 1)[None, :]
        tot = k1 + k2
        self.valid = tot <= n
        safe = np.minimum(tot, n)
        # SB = p_0 + v1*(p_{k1} - p_0) + v2*(p_{k1+k2} - p_{k1})
        self.sb = p[0] + v1 * (p[k1] - p[0]) + v2 * (p[safe] - p[k1])
        self.sb = np.where(self.valid, self.sb, np.inf)
        self.k1, self.k2, self.k0 = k1, k2, n - safe
        self._logcoef = (_lgamma_fact(n) - _lgamma_fact(k1) - _lgamma_fact(k2)
                         - _lgamma_fact(np.where(self.valid, n - tot, 0)))

    def coverage(self, q1: float, q2: float) -> float:
        q0 = 1.0 - q1 - q2
        if q0 < -1e-15 or q1 < 0 or q2 < 0:
            raise ValueError("invalid probability vector")
        q0 = max(q0, 0.0)
        theta = q1 * self.v1 + q2 * self.v2
        with np.errstate(divide="ignore", invalid="ignore"):
            logw = (self._logcoef
                    + _xlogy(self.k1, q1) + _xlogy(self.k2, q2)
                    + _xlogy(np.where(self.valid, self.k0, 0), q0))
        w = np.where(self.valid, np.exp(logw), 0.0)
        return float(w[self.sb >= theta].sum())

    def failure_mass(self, q1: float, q2: float) -> np.ndarray:
        """Weight matrix restricted to the non-covering count vectors."""
        q0 = max(1.0 - q1 - q2, 0.0)
        theta = q1 * self.v1 + q2 * self.v2
        with np.errstate(divide="ignore", invalid="ignore"):
            logw = (self._logcoef
                    + _xlogy(self.k1, q1) + _xlogy(self.k2, q2)
                    + _xlogy(np.where(self.valid, self.k0, 0), q0))
        w = np.where(self.valid, np.exp(logw), 0.0)
        return np.where(self.sb < theta, w, 0.0)


def _lgamma_fact(k):
    from scipy.special import gammaln
    return gammaln(np.asarray(k, dtype=float) + 1.0)


def _xlogy(k, q):
    k = np.asarray(k, dtype=float)
    if q <= 0:
        return np.where(k > 0, -np.inf, 0.0)
    return k * np.log(q)


__all__ = [
    "float_factors",
    "coverage_exact",
    "TwoValueGrid",
    "UncertifiedComparison",
]
