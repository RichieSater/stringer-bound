"""Stringer bound: confidence factors and the bound itself.

Reference definition (Stringer 1963; Leslie--Teitlebaum--Anderson 1979, ch. 4).

Monetary-unit sampling draws ``n`` dollar units at random (with replacement)
from a population whose taint at a dollar unit is ``T in [0, 1]``, where
``T = (book value - audited value) / book value`` of the line item that owns
the dollar unit.  The parameter of interest is the population overstatement
fraction ``theta = E[T]``.

With ``t_(1) >= t_(2) >= ... >= t_(n)`` the observed taints sorted in
decreasing order, the Stringer upper bound at nominal confidence ``1 - alpha``
is

    SB = p_0 + sum_{j >= 1} (p_j - p_{j-1}) * t_(j),

where ``p_j`` is the upper ``1 - alpha`` confidence limit for a binomial
proportion after observing ``j`` errors in ``n`` trials, i.e. the root of

    sum_{i=0}^{j} C(n, i) p^i (1 - p)^{n-i} = alpha.

Audit practice usually substitutes the Poisson limits (the AICPA "MUS
factors"), ``p_j = m_j / n`` with ``sum_{i<=j} e^{-m} m^i / i! = alpha``.
Both are provided; every claim in this repository states which one it uses.

All factors are computed in mpmath at high precision by bisection on a
provably monotone function, so each ``p_j`` is returned together with a
rigorous enclosing interval width.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable, List, Sequence, Tuple

from mpmath import mp, mpf, binomial, exp, factorial

DEFAULT_DPS = 50
#: bisection is stopped when the bracket is narrower than 10**-BISECT_DIGITS
BISECT_DIGITS = 40

METHODS = ("binomial", "poisson")


def _tail_binomial(p, n: int, j: int):
    """P(Bin(n, p) <= j), computed exactly in the current mpmath precision."""
    s = mpf(0)
    q = 1 - p
    for i in range(j + 1):
        s += binomial(n, i) * p ** i * q ** (n - i)
    return s


def _tail_poisson(m, j: int):
    """P(Poisson(m) <= j)."""
    s = mpf(0)
    for i in range(j + 1):
        s += m ** i / factorial(i)
    return s * exp(-m)


def _bisect(f, lo, hi, digits: int):
    """Bisect a strictly decreasing f with f(lo) > 0 > f(hi).

    Returns (midpoint, bracket_width).
    """
    lo, hi = mpf(lo), mpf(hi)
    target = mpf(10) ** (-digits)
    flo, fhi = f(lo), f(hi)
    if flo < 0 or fhi > 0:
        raise ValueError("bisection bracket does not straddle the root")
    while hi - lo > target:
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2, hi - lo


@lru_cache(maxsize=None)
def factors(n: int, alpha: str, method: str = "binomial",
            dps: int = DEFAULT_DPS) -> Tuple[Tuple, ...]:
    """Confidence factors ``p_0, ..., p_n`` and their bracket widths.

    ``alpha`` is passed as a string (e.g. ``"0.05"``) so the cache key and the
    decimal value are both exact.  Returns a tuple of ``(p_j, width)`` pairs.
    """
    if method not in METHODS:
        raise ValueError("method must be one of %r" % (METHODS,))
    with mp.workdps(dps):
        a = mpf(alpha)
        if not (0 < a < 1):
            raise ValueError("alpha must lie strictly between 0 and 1")
        out = []
        for j in range(n + 1):
            if method == "binomial":
                if j == n:
                    # sum_{i<=n} ... == 1 > alpha for every p: the bound is 1.
                    out.append((mpf(1), mpf(0)))
                    continue
                f = lambda p, j=j: _tail_binomial(p, n, j) - a
                root, width = _bisect(f, 0, 1, BISECT_DIGITS)
            else:
                f = lambda m, j=j: _tail_poisson(m, j) - a
                hi = mpf(1)
                while f(hi) > 0:
                    hi *= 2
                root, width = _bisect(f, 0, hi, BISECT_DIGITS)
                root, width = root / n, width / n
            out.append((root, width))
        return tuple(out)


def factor_values(n: int, alpha: str, method: str = "binomial",
                  dps: int = DEFAULT_DPS) -> List:
    return [p for p, _ in factors(n, alpha, method, dps)]


def max_factor_width(n: int, alpha: str, method: str = "binomial",
                     dps: int = DEFAULT_DPS):
    return max(w for _, w in factors(n, alpha, method, dps))


def stringer_bound(taints: Iterable[float], n: int, alpha: str,
                   method: str = "binomial", dps: int = DEFAULT_DPS):
    """Stringer bound from a raw sample of ``n`` taints (zeros may be omitted).

    Taints need not be sorted.  Zero taints contribute nothing and, because
    the sort is decreasing, do not displace any nonzero taint, so omitting
    them leaves the bound unchanged.
    """
    p = factor_values(n, alpha, method, dps)
    ts = sorted((mpf(t) for t in taints), reverse=True)
    if len(ts) > n:
        raise ValueError("more taints than sample units")
    with mp.workdps(dps):
        total = p[0]
        for j, t in enumerate(ts, start=1):
            total += (p[j] - p[j - 1]) * t
        return total


def bound_from_counts(values: Sequence, counts: Sequence[int],
                      p: Sequence) -> object:
    """Stringer bound when the sample contains ``counts[i]`` taints of size
    ``values[i]``, with ``values`` strictly decreasing and positive.

    Uses the telescoped form
        SB = p_0 + sum_i v_i * (p_{c_i} - p_{c_{i-1}}),
    with ``c_i`` the number of sampled taints greater than or equal to
    ``v_i``.  This is exactly the sorted-taint definition, regrouped.
    """
    total = p[0]
    c = 0
    for v, k in zip(values, counts):
        if k:
            prev, c = c, c + k
            total += v * (p[c] - p[prev])
    return total


__all__ = [
    "factors",
    "factor_values",
    "max_factor_width",
    "stringer_bound",
    "bound_from_counts",
    "METHODS",
]
