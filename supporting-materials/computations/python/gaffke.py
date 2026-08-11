"""Gaffke upper bound and an all-sample-size Stringer safeguard.

For observations ``x_1, ..., x_n`` in ``[0, 1]``, the one-sided Gaffke--
Learned-Miller--Thomas upper confidence limit at tail probability ``alpha``
is the ``(1-alpha)`` quantile of

    D_0 + sum_i x_i D_i,

where ``(D_0, ..., D_n)`` is uniform on the ``n``-simplex, equivalently
Dirichlet(1, ..., 1).  Vlassis and Thomas (2026) proved that this is a
distribution-free upper confidence limit for the mean of independent
``[0, 1]``-valued observations.

Consequently,

    max(Stringer bound, Gaffke upper bound)

has finite-sample coverage at least ``1-alpha`` for every ``n``.  This
``safeguarded_stringer_bound`` is not a proof that the ordinary Stringer
bound is conservative.  It is an operational fallback: it changes the
reported value only when the independently valid Gaffke limit is larger.

The Dirichlet-average density is the normalized B-spline with knots
``(x_1, ..., x_n, 1)``.  SciPy supplies a fast numerical root locator.  The
reported Gaffke endpoint is nevertheless certified without trusting that
locator: decimal inputs are converted to exact rationals, the Dirichlet tail
is evaluated as a confluent divided difference in ``Fraction`` arithmetic,
and the two endpoints of a dyadic quantile bracket have their signs checked
exactly.  The upper endpoint is therefore a conservative rational enclosure
of the mathematically defined Gaffke quantile.  This exact-sign layer is what
makes the command-line safeguard suitable for reproducible reporting.

Command-line example (zero taints may be omitted):

    uv run --frozen python \
      supporting-materials/computations/python/gaffke.py \
      --n 100 --alpha 0.05 --method poisson --taints 1,0.4,0.1
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import asdict, dataclass
from fractions import Fraction
from math import comb, isfinite
from typing import Iterable, Sequence

import numpy as np
from scipy.interpolate import BSpline

from stringer import METHODS, stringer_bound


DEFAULT_CERTIFICATE_BITS = 48


def _as_fraction(value: object) -> Fraction:
    """Convert a decimal-like input without inheriting its binary expansion."""
    if isinstance(value, Fraction):
        return value
    # ``str`` preserves an intended command-line decimal and gives a compact
    # decimal interpretation of Python and NumPy floating-point values.
    return Fraction(str(value))


def _complete_sample_fraction(
        taints: Iterable[object], n: int) -> tuple[Fraction, ...]:
    """Validate a zero-omitted sample and return all taints exactly."""
    if n < 1:
        raise ValueError("n must be positive")
    values = [_as_fraction(value) for value in taints]
    if len(values) > n:
        raise ValueError("more taints than sample units")
    if any(value < 0 or value > 1 for value in values):
        raise ValueError("every taint must lie in [0, 1]")
    values.extend([Fraction(0)] * (n - len(values)))
    return tuple(sorted(values))


def _fraction_to_float_up(value: Fraction) -> float:
    """Convert a finite rational to binary64 without rounding downward."""
    converted = float(value)
    if Fraction.from_float(converted) < value:
        converted = float(np.nextafter(converted, np.inf))
    return converted


def _fraction_to_float_down(value: Fraction) -> float:
    """Convert a finite rational to binary64 without rounding upward."""
    converted = float(value)
    if Fraction.from_float(converted) > value:
        converted = float(np.nextafter(converted, -np.inf))
    return converted


def dirichlet_average_cdf(value: float, knots: Sequence[float]) -> float:
    """CDF of a uniform-Dirichlet average of the supplied knots.

    If ``m`` knots are supplied, the weights have the
    ``Dirichlet(1, ..., 1)`` law of length ``m``.  The implementation uses
    the normalized B-spline density and therefore remains stable when knots
    repeat.  An all-equal knot vector is treated as a point mass.
    """
    ordered = np.sort(np.asarray(knots, dtype=float))
    if ordered.ndim != 1 or len(ordered) < 2:
        raise ValueError("at least two finite knots are required")
    if np.any(~np.isfinite(ordered)):
        raise ValueError("knots must be finite")
    lower, upper = float(ordered[0]), float(ordered[-1])
    if value < lower:
        return 0.0
    if value >= upper:
        return 1.0
    if lower == upper:
        return float(value >= lower)

    degree_plus_one = len(ordered) - 1
    density_spline = BSpline.basis_element(ordered, extrapolate=False)
    antiderivative = density_spline.antiderivative()
    normalizer = degree_plus_one / (upper - lower)
    cdf = normalizer * (antiderivative(value) - antiderivative(lower))
    # Roundoff can put an endpoint a few ulps outside the probability range.
    return float(np.clip(cdf, 0.0, 1.0))


def dirichlet_average_quantile(
        knots: Sequence[float], probability: float,
        *, iterations: int = 80) -> float:
    """Quantile of a uniform-Dirichlet average, rounded upward.

    The returned value is the upper endpoint of a fixed-iteration bisection
    bracket, followed by one ``nextafter`` step toward ``+infinity``.  This
    avoids accidentally reporting a value just below the numerical root.
    """
    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must lie in [0, 1]")
    if iterations < 1:
        raise ValueError("iterations must be positive")
    ordered = np.sort(np.asarray(knots, dtype=float))
    if ordered.ndim != 1 or len(ordered) < 2:
        raise ValueError("at least two knots are required")
    if np.any(~np.isfinite(ordered)):
        raise ValueError("knots must be finite")
    lower, upper = float(ordered[0]), float(ordered[-1])
    if lower == upper or probability <= 0.0:
        return lower
    if probability >= 1.0:
        return upper

    # Construct the spline once rather than rebuilding it on every CDF call.
    density_spline = BSpline.basis_element(ordered, extrapolate=False)
    antiderivative = density_spline.antiderivative()
    base = antiderivative(lower)
    normalizer = (len(ordered) - 1) / (upper - lower)

    for _ in range(iterations):
        midpoint = (lower + upper) / 2.0
        cdf = normalizer * (antiderivative(midpoint) - base)
        if cdf < probability:
            lower = midpoint
        else:
            upper = midpoint
    return float(np.nextafter(upper, np.inf))


def dirichlet_average_tail_exact(
        value: object, knots: Sequence[object]) -> Fraction:
    r"""Return an exact Dirichlet-average upper tail at a rational point.

    For ``m=n+1`` knots ``a_0,...,a_n`` and uniform Dirichlet weights,

    .. math::

       \Pr\{\textstyle\sum_i a_iD_i>x\}
       =[a_0,\ldots,a_n](t-x)_+^n.

    The right side is a divided difference.  Repeated knots are evaluated
    confluently, using the exact derivative divided by its factorial.  All
    operations below are therefore rational when the supplied values are
    decimal strings, integers, or :class:`fractions.Fraction` objects.
    """
    x = _as_fraction(value)
    ordered = tuple(sorted(_as_fraction(knot) for knot in knots))
    if len(ordered) < 2:
        raise ValueError("at least two knots are required")

    if ordered[0] == ordered[-1]:
        # The Dirichlet average is a point mass.  This is P(Y > x), which is
        # the tail convention used by the continuous nondegenerate formula.
        return Fraction(int(ordered[0] > x))
    if x < ordered[0]:
        return Fraction(1)
    if x >= ordered[-1]:
        return Fraction(0)

    degree = len(ordered) - 1
    multiplicities = tuple(sorted(Counter(ordered).items()))

    # A confluent divided difference is the sum of the residues of
    #
    #     (z-x)_+^degree / prod_b (z-b)^multiplicity[b].
    #
    # At a knot a>x of multiplicity m, its residue is the coefficient of
    # y^(m-1) in
    #
    #   (a+y-x)^degree * prod_{b!=a}(a-b+y)^(-multiplicity[b]).
    #
    # Computing only that truncated coefficient is dramatically faster than
    # a full O(n^2) divided-difference table for the zero-heavy samples usual
    # in audit work, while remaining exact and handling repeated positive
    # taints.  Knots at or below x contribute zero because the truncated
    # power has all required derivatives equal to zero there.
    tail = Fraction(0)
    for knot, multiplicity in multiplicities:
        if knot <= x:
            continue
        target_degree = multiplicity - 1
        series = [
            Fraction(comb(degree, order))
            * (knot - x) ** (degree - order)
            for order in range(target_degree + 1)
        ]
        for other, other_multiplicity in multiplicities:
            if other == knot:
                continue
            separation = knot - other
            factor = [
                Fraction(
                    (-1) ** order
                    * comb(other_multiplicity + order - 1, order),
                    1,
                )
                / separation ** (other_multiplicity + order)
                for order in range(target_degree + 1)
            ]
            product: list[Fraction] = []
            for total_degree in range(target_degree + 1):
                product.append(sum(
                    series[left_degree]
                    * factor[total_degree - left_degree]
                    for left_degree in range(total_degree + 1)
                ))
            series = product
        tail += series[target_degree]

    if not 0 <= tail <= 1:
        raise ArithmeticError("exact Dirichlet tail left [0, 1]")
    return tail


@dataclass(frozen=True)
class GaffkeQuantileCertificate:
    """Exact-sign dyadic enclosure of a one-sided Gaffke quantile."""

    alpha: Fraction
    lower: Fraction
    upper: Fraction
    tail_at_lower: Fraction
    tail_at_upper: Fraction
    grid_bits: int
    degenerate: bool = False

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower


def gaffke_quantile_certificate(
        taints: Iterable[object], n: int, alpha: object,
        *, grid_bits: int = DEFAULT_CERTIFICATE_BITS
) -> GaffkeQuantileCertificate:
    """Certify a dyadic bracket for the Gaffke upper endpoint.

    A numerical B-spline quantile is used only to propose a narrow dyadic
    bracket.  Both tail signs are then evaluated exactly.  If the proposal
    is inaccurate, the bracket expands until the exact signs establish

    ``tail(lower) >= alpha >= tail(upper)``.

    Thus ``upper`` is always a conservative enclosure of the target
    ``(1-alpha)`` quantile.  The default grid has width on the order of
    ``2**-48``, negligible for audit reporting but inexpensive to certify.
    """
    if grid_bits < 8:
        raise ValueError("grid_bits must be at least 8")
    tail_probability = _as_fraction(alpha)
    if not 0 < tail_probability < 1:
        raise ValueError("alpha must lie strictly between zero and one")

    sample = _complete_sample_fraction(taints, n)
    knots = sample + (Fraction(1),)
    if knots[0] == knots[-1]:
        point = knots[0]
        return GaffkeQuantileCertificate(
            alpha=tail_probability,
            lower=point,
            upper=point,
            tail_at_lower=Fraction(0),
            tail_at_upper=Fraction(0),
            grid_bits=grid_bits,
            degenerate=True,
        )

    numerical = dirichlet_average_quantile(
        [float(knot) for knot in knots],
        1.0 - float(tail_probability),
    )
    denominator = 1 << grid_bits
    center = int(Fraction.from_float(numerical) * denominator)
    center = min(max(center, 0), denominator)

    radius = 2
    while True:
        lower_numerator = max(0, center - radius)
        upper_numerator = min(denominator, center + radius + 1)
        lower = Fraction(lower_numerator, denominator)
        upper = Fraction(upper_numerator, denominator)
        tail_at_lower = dirichlet_average_tail_exact(lower, knots)
        tail_at_upper = dirichlet_average_tail_exact(upper, knots)
        if (tail_at_lower >= tail_probability
                and tail_at_upper <= tail_probability):
            return GaffkeQuantileCertificate(
                alpha=tail_probability,
                lower=lower,
                upper=upper,
                tail_at_lower=tail_at_lower,
                tail_at_upper=tail_at_upper,
                grid_bits=grid_bits,
            )
        if lower_numerator == 0 and upper_numerator == denominator:
            raise ArithmeticError("failed to bracket the Gaffke quantile")
        radius *= 2


def gaffke_upper_bound(
        taints: Iterable[object], n: int, alpha: object,
        *, grid_bits: int = DEFAULT_CERTIFICATE_BITS) -> float:
    """Certified upper enclosure of the one-sided Gaffke endpoint."""
    certificate = gaffke_quantile_certificate(
        taints, n, alpha, grid_bits=grid_bits)
    return _fraction_to_float_up(certificate.upper)


@dataclass(frozen=True)
class SafeguardedResult:
    """Components of the safeguarded report."""

    n: int
    alpha: float
    method: str
    stringer: float
    gaffke_lower: float
    gaffke: float
    gaffke_bracket_width: float
    gaffke_certificate_bits: int
    gaffke_upper_dyadic: str
    safeguarded: float
    uplift: float
    governing_bound: str


def safeguarded_stringer_bound(
        taints: Iterable[object], n: int, alpha: object,
        method: str = "poisson", dps: int = 80,
        *, certificate_bits: int = DEFAULT_CERTIFICATE_BITS
) -> SafeguardedResult:
    """Return ``max(Stringer, Gaffke)`` and its two components.

    The default is the Poisson-factor Stringer bound used in audit practice.
    The safeguard's coverage follows from the Gaffke component for every
    sample size; no unresolved general-``n`` Stringer claim is assumed.
    """
    if method not in METHODS:
        raise ValueError("method must be one of %r" % (METHODS,))
    # Materialize once because callers may pass a generator.
    observed = tuple(taints)
    _complete_sample_fraction(observed, n)
    sb = float(stringer_bound(observed, n, str(alpha), method, dps))
    if not isfinite(sb):
        raise ArithmeticError("Stringer calculation was not finite")
    certificate = gaffke_quantile_certificate(
        observed, n, alpha, grid_bits=certificate_bits)
    gaffke_lower = _fraction_to_float_down(certificate.lower)
    gaffke = _fraction_to_float_up(certificate.upper)
    safeguarded = max(sb, gaffke)
    return SafeguardedResult(
        n=n,
        alpha=float(alpha),
        method=method,
        stringer=sb,
        gaffke_lower=gaffke_lower,
        gaffke=gaffke,
        gaffke_bracket_width=float(certificate.width),
        gaffke_certificate_bits=certificate.grid_bits,
        gaffke_upper_dyadic=(
            f"{certificate.upper.numerator}/{certificate.upper.denominator}"
        ),
        safeguarded=safeguarded,
        uplift=safeguarded - sb,
        governing_bound="stringer" if sb >= gaffke else "gaffke",
    )


def _parse_taints(raw: str) -> list[str]:
    if not raw.strip():
        return []
    return [piece.strip() for piece in raw.split(",")]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute Stringer, Gaffke, and their valid maximum.")
    parser.add_argument("--n", type=int, required=True,
                        help="sample size; zero taints may be omitted")
    parser.add_argument("--alpha", default="0.05",
                        help="tail probability, e.g. 0.05 for 95%%")
    parser.add_argument("--method", choices=METHODS, default="poisson")
    parser.add_argument("--certificate-bits", type=int,
                        default=DEFAULT_CERTIFICATE_BITS,
                        help="dyadic grid bits for exact Gaffke sign checks")
    parser.add_argument("--taints", default="",
                        help="comma-separated observed taints; omit zeros")
    args = parser.parse_args(argv)
    result = safeguarded_stringer_bound(
        _parse_taints(args.taints), args.n, args.alpha, args.method,
        certificate_bits=args.certificate_bits)
    print(json.dumps(asdict(result), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DEFAULT_CERTIFICATE_BITS",
    "GaffkeQuantileCertificate",
    "SafeguardedResult",
    "dirichlet_average_cdf",
    "dirichlet_average_quantile",
    "dirichlet_average_tail_exact",
    "gaffke_quantile_certificate",
    "gaffke_upper_bound",
    "safeguarded_stringer_bound",
]
