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

``coverage_exact``  rational arithmetic (``fractions.Fraction``) for the
                    multinomial weights, high-precision mpmath for ``SB``.
                    Returns the coverage as an exact rational together with
                    the smallest ``|SB - theta|`` seen, which certifies that
                    every comparison was decided far outside the numerical
                    error of the confidence factors.

``coverage_fast``   float64/numpy, specialised to support ``{v1 > v2 > 0}``
                    plus an atom at 0.  Used for search; every candidate it
                    reports is re-checked with ``coverage_exact``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb
from typing import List, Sequence, Tuple

import numpy as np
from mpmath import mp, mpf

from stringer import bound_from_counts, factors


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


def coverage_exact(values: Sequence[Fraction], probs: Sequence[Fraction],
                   n: int, alpha: str, method: str = "binomial",
                   dps: int = 50) -> Tuple[Fraction, object, object]:
    """Exact coverage for a finitely supported taint distribution.

    ``values``  strictly decreasing positive rational taints.
    ``probs``   their probabilities; ``1 - sum(probs)`` is the atom at 0.

    Returns ``(coverage, theta, min_margin)`` where ``coverage`` is an exact
    ``Fraction``, ``theta`` is the exact mean taint, and ``min_margin`` is
    ``min |SB - theta|`` over all enumerated count vectors.
    """
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
    p = [v for v, _ in factors(n, alpha, method, dps)]

    m = len(values)
    covered = Fraction(0)
    min_margin = None
    with mp.workdps(dps):
        theta_mp = mpf(theta.numerator) / mpf(theta.denominator)
        for counts in _count_vectors(n, m):
            k0 = n - sum(counts)
            weight = Fraction(1)
            rem = n
            for c in counts:
                weight *= comb(rem, c)
                rem -= c
            for c, pr in zip(counts, probs):
                weight *= pr ** c
            weight *= q0 ** k0
            if weight == 0:
                continue
            sb = bound_from_counts(values, counts, p)
            margin = abs(sb - theta_mp)
            if min_margin is None or margin < min_margin:
                min_margin = margin
            if sb >= theta_mp:
                covered += weight
    return covered, theta, min_margin


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


__all__ = ["float_factors", "coverage_exact", "TwoValueGrid"]
