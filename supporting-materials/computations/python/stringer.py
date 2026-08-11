"""Stringer bound: confidence factors and the bound itself.

Reference definition (Stringer 1963; Leslie--Teitlebaum--Anderson 1979, ch. 4).

Monetary-unit sampling draws ``n`` dollar units at random (with replacement)
from a population whose taint at a dollar unit is ``T in [0, 1]``, where
``T = (book value - audited value) / book value`` of the account item
associated with the dollar unit.  The parameter of interest is the
population overstatement fraction ``theta = E[T]``.

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
For ``alpha < exp(-1)``, Anderson--Samuels (1967) implies that the Poisson
factors dominate the binomial factors for every ``n`` and ``j``; see
``supporting-materials/theory/POISSON-DOMINATION.md``.

The search routines compute numerical factors in ``mpmath``.  Formal
binomial certificates use a separate exact routine: it bisects on a dyadic
grid, evaluates the binomial-CDF sign with integer arithmetic, and returns
rational lower and upper bounds for every factor.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import ceil
from typing import Iterable, List, Sequence, Tuple

from mpmath import mp, mpf, binomial, exp, factorial

DEFAULT_DPS = 50
#: bisection is stopped when the bracket is narrower than 10**-BISECT_DIGITS
BISECT_DIGITS = 40
#: denominator of the default exact dyadic factor brackets is ``2**bits``.
EXACT_FACTOR_BITS = 80

METHODS = ("binomial", "poisson")


def _tail_binomial(p, n: int, j: int):
    """Numerical ``P(Bin(n, p) <= j)`` at the current mpmath precision."""
    if p == 0:
        return mpf(1)
    if p == 1:
        return mpf(1 if j == n else 0)
    q = 1 - p
    term = q ** n
    s = term
    for i in range(j):
        term *= mpf(n - i) * p / (mpf(i + 1) * q)
        s += term
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


def _binomial_cdf_sign_dyadic(k: int, bits: int, n: int, j: int,
                              alpha: Fraction) -> int:
    """Exact sign of ``P(Bin(n, k/2**bits) <= j) - alpha``.

    Every term is put over the common denominator ``2**(bits*n)``.
    The shorter of the lower tail and its complementary upper tail is
    accumulated.  Thus the returned sign depends only on integer arithmetic,
    not on a floating-point evaluation of the binomial CDF.
    """
    denominator = 1 << bits
    if not 0 <= k <= denominator:
        raise ValueError("dyadic numerator is outside [0, 2**bits]")
    if not 0 <= j < n:
        raise ValueError("j must satisfy 0 <= j < n")

    if k == 0:
        return 1  # 1 - alpha > 0
    if k == denominator:
        return -1  # 0 - alpha < 0

    q = denominator - k
    scale = 1 << (bits * n)
    aden, anum = alpha.denominator, alpha.numerator

    if j + 1 <= n - j:
        # Lower tail, starting from the i=0 term q**n.
        term = q ** n
        total = term
        for i in range(j):
            term, remainder = divmod(
                term * (n - i) * k, (i + 1) * q)
            if remainder:
                raise ArithmeticError("nonintegral binomial recurrence")
            total += term
        delta = total * aden - anum * scale
    else:
        # Upper complement, starting from the i=n term k**n.
        term = k ** n
        upper = term
        for i in range(n, j + 1, -1):
            term, remainder = divmod(
                term * i * q, (n - i + 1) * k)
            if remainder:
                raise ArithmeticError("nonintegral binomial recurrence")
            upper += term
        delta = (aden - anum) * scale - upper * aden

    return (delta > 0) - (delta < 0)


def _approx_binomial_root(n: int, j: int, alpha: Fraction, bits: int):
    """Locate a binomial factor numerically before exact sign checks.

    The approximation affects speed only.  ``exact_binomial_factor_brackets``
    verifies the final endpoints by exact integer arithmetic and falls back
    to an exact grid bisection if this locator is inaccurate.
    """
    from scipy.special import betaincinv

    dps = max(DEFAULT_DPS, ceil(bits * 0.30103) + 25)
    with mp.workdps(dps):
        a = mpf(alpha.numerator) / alpha.denominator
        guess = betaincinv(j + 1, n - j, 1.0 - float(alpha))
        p = mpf(float(guess))
        if not 0 < p < 1:
            p = mpf(j + 1) / (n + 1)
        for _ in range(4):
            tail = _tail_binomial(p, n, j)
            derivative = (-(n - j) * binomial(n, j) * p ** j
                          * (1 - p) ** (n - j - 1))
            candidate = p - (tail - a) / derivative
            if not 0 < candidate < 1:
                break
            p = candidate
        return +p


@lru_cache(maxsize=None)
def exact_binomial_factor_brackets(
        n: int, alpha: str, bits: int = EXACT_FACTOR_BITS
) -> Tuple[Tuple[Fraction, Fraction], ...]:
    """Return rigorous rational brackets for ``p_0, ..., p_n``.

    For ``j < n``, the endpoints lie on the dyadic grid with denominator
    ``2**bits`` and their binomial-CDF signs are evaluated exactly.  The
    numerical root locator is never trusted: if its adjacent grid points do
    not straddle the root, a full exact bisection of the grid is used.
    ``p_n`` is returned as the exact interval ``[1, 1]``.
    """
    if n < 1:
        raise ValueError("n must be positive")
    if bits < 1:
        raise ValueError("bits must be positive")
    # Nominal alpha is a decimal design parameter, not a sampled float.
    # Converting through ``str`` preserves that intended decimal value even
    # if a caller supplied a Python float rather than the preferred string.
    a = Fraction(str(alpha))
    if not 0 < a < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")

    denominator = 1 << bits
    out = []
    for j in range(n):
        approximate = _approx_binomial_root(n, j, a, bits)
        with mp.workdps(max(DEFAULT_DPS, ceil(bits * 0.30103) + 25)):
            k = int(mp.floor(approximate * denominator))
        k = min(max(k, 0), denominator - 1)
        sign_lo = _binomial_cdf_sign_dyadic(k, bits, n, j, a)
        sign_hi = _binomial_cdf_sign_dyadic(k + 1, bits, n, j, a)

        if sign_lo == 0:
            lo = hi = k
        elif sign_hi == 0:
            lo = hi = k + 1
        elif sign_lo > 0 and sign_hi < 0:
            lo, hi = k, k + 1
        else:
            # The approximate locator is only an optimization.  This branch
            # establishes the bracket from scratch using exact signs.
            lo, hi = 0, denominator
            while hi - lo > 1:
                mid = (lo + hi) // 2
                sign_mid = _binomial_cdf_sign_dyadic(
                    mid, bits, n, j, a)
                if sign_mid == 0:
                    lo = hi = mid
                    break
                if sign_mid > 0:
                    lo = mid
                else:
                    hi = mid

        out.append((Fraction(lo, denominator),
                    Fraction(hi, denominator)))
    out.append((Fraction(1), Fraction(1)))
    return tuple(out)


@lru_cache(maxsize=None)
def factor_prefix(n: int, alpha: str, maximum_j: int,
                  method: str = "binomial",
                  dps: int = DEFAULT_DPS) -> Tuple[Tuple, ...]:
    """Numerical confidence factors through ``maximum_j``.

    ``alpha`` is passed as a string (e.g. ``"0.05"``) so the cache key and the
    decimal value are both exact.  These values are used for screening and
    Poisson calculations, not for formal binomial certificates.  Use
    :func:`exact_binomial_factor_brackets` for the latter.

    Computing only a prefix matters in audit samples with few nonzero
    taints: the Stringer formula then needs only ``p_0,...,p_k`` rather than
    all ``n+1`` factors.
    """
    if n < 1:
        raise ValueError("n must be positive")
    if not 0 <= maximum_j <= n:
        raise ValueError("maximum_j must lie between zero and n")
    if method not in METHODS:
        raise ValueError("method must be one of %r" % (METHODS,))
    with mp.workdps(dps):
        a = mpf(alpha)
        if not (0 < a < 1):
            raise ValueError("alpha must lie strictly between 0 and 1")
        out = []
        for j in range(maximum_j + 1):
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


@lru_cache(maxsize=None)
def factors(n: int, alpha: str, method: str = "binomial",
            dps: int = DEFAULT_DPS) -> Tuple[Tuple, ...]:
    """All numerical confidence factors and final bisection widths."""
    return factor_prefix(n, alpha, n, method, dps)


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
    raw_taints = [mpf(t) for t in taints]
    if len(raw_taints) > n:
        raise ValueError("more taints than sample units")
    ts = sorted((t for t in raw_taints if t != 0), reverse=True)
    p = [value for value, _ in factor_prefix(
        n, alpha, len(ts), method, dps)]
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
    "factor_prefix",
    "factor_values",
    "max_factor_width",
    "stringer_bound",
    "bound_from_counts",
    "exact_binomial_factor_brackets",
    "EXACT_FACTOR_BITS",
    "METHODS",
]
