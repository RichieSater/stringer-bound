"""Directed scalar lemmas used by the ``n=7`` Poissonization proof."""

from __future__ import annotations

import hashlib
import itertools
import math
from collections.abc import Callable
from fractions import Fraction

from flint import arb, ctx, fmpq

ARB_PRECISION_BITS = 160
SCALAR_DOMAIN_END = Fraction(7)
SMALL_U_CUTOFF = Fraction(1, 8)

FACE6_Q4_EVENTS = (
    (Fraction(3, 8), Fraction(7), Fraction(1, 20)),
    (Fraction(1, 2), Fraction(43, 20), Fraction(1, 20)),
    (Fraction(14, 25), Fraction(2), Fraction(1, 10)),
    (Fraction(16, 25), Fraction(9, 5), Fraction(1, 10)),
    (Fraction(7, 10), Fraction(17, 10), Fraction(1, 10)),
    (Fraction(77, 100), Fraction(8, 5), Fraction(1, 10)),
    (Fraction(169, 200), Fraction(3, 2), Fraction(1, 10)),
    (Fraction(187, 200), Fraction(21, 20), Fraction(1, 25)),
)
FACE7_Q5_EVENTS = (
    (Fraction(1, 3), Fraction(6), Fraction(1, 20)),
    (Fraction(11, 25), Fraction(3), Fraction(1, 20)),
    (Fraction(49, 100), Fraction(5, 3), Fraction(1, 10)),
    (Fraction(11, 20), Fraction(8, 5), Fraction(1, 10)),
    (Fraction(59, 100), Fraction(3, 2), Fraction(1, 5)),
    (Fraction(67, 100), Fraction(10, 7), Fraction(3, 10)),
    (Fraction(83, 100), Fraction(4, 3), Fraction(1, 5)),
    (Fraction(27, 25), Fraction(63, 50), Fraction(1, 10)),
    (Fraction(109, 100), Fraction(49, 40), Fraction(1, 10)),
    (Fraction(28, 25), Fraction(6, 5), Fraction(1, 40)),
)


def _arb_exact(value: Fraction) -> arb:
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def _arb_interval(lower: Fraction, upper: Fraction) -> arb:
    lower = Fraction(lower)
    upper = Fraction(upper)
    lower_q = fmpq(lower.numerator, lower.denominator)
    upper_q = fmpq(upper.numerator, upper.denominator)
    return arb((lower_q + upper_q) / 2, (upper_q - lower_q) / 2)


def _q_arb(value: arb, order: int, upper_piece: bool) -> arb:
    lam = 7 / value
    result = (-lam).exp() * sum(
        (lam**index / math.factorial(index) for index in range(order + 1)),
        arb(0),
    )
    if upper_piece:
        corrections = {
            3: -1 + 35 / value**4 - 84 / value**5 + 70 / value**6 - 20 / value**7,
            4: -1 + 21 / value**5 - 35 / value**6 + 15 / value**7,
            5: -1 + 7 / value**6 - 6 / value**7,
            6: -1 + value**-7,
        }
        result += corrections[order]
    return result


def _q_derivative_arb(value: arb, order: int, upper_piece: bool) -> arb:
    lam = 7 / value
    result = (-lam).exp() * lam ** (order + 2) / (7 * math.factorial(order))
    if upper_piece:
        corrections = {
            3: -140 / value**5 + 420 / value**6 - 420 / value**7 + 140 / value**8,
            4: -105 / value**6 + 210 / value**7 - 105 / value**8,
            5: -42 / value**7 + 42 / value**8,
            6: -7 / value**8,
        }
        result += corrections[order]
    return result


def _symbolic_checks() -> None:
    from sympy import diff, exp, factorial, simplify, symbols

    u = symbols("u", positive=True)
    lam = 7 / u
    for order in range(3, 7):
        lower_f = u**order * exp(-lam)
        upper_f = lower_f - (u - 1) ** 7 / u ** (7 - order)
        lower_q = exp(-lam) * sum(
            lam**index / factorial(index) for index in range(order + 1)
        )
        corrections = {
            3: -1 + 35 / u**4 - 84 / u**5 + 70 / u**6 - 20 / u**7,
            4: -1 + 21 / u**5 - 35 / u**6 + 15 / u**7,
            5: -1 + 7 / u**6 - 6 / u**7,
            6: -1 + u**-7,
        }
        upper_q = lower_q + corrections[order]
        lower_derivative = exp(-lam) * lam ** (order + 2) / (7 * factorial(order))
        derivative_corrections = {
            3: -140 / u**5 + 420 / u**6 - 420 / u**7 + 140 / u**8,
            4: -105 / u**6 + 210 / u**7 - 105 / u**8,
            5: -42 / u**7 + 42 / u**8,
            6: -7 / u**8,
        }
        assert simplify(diff(lower_f, u, order) / factorial(order) - lower_q) == 0
        assert simplify(diff(upper_f, u, order) / factorial(order) - upper_q) == 0
        assert simplify(diff(lower_q, u) - lower_derivative) == 0
        assert (
            simplify(
                diff(upper_q, u) - lower_derivative - derivative_corrections[order]
            )
            == 0
        )


def verify_n7_scalar_bounds() -> dict[str, object]:
    """Regenerate every directed scalar sign used by the four faces."""
    _symbolic_checks()
    previous_precision = ctx.prec
    ctx.prec = ARB_PRECISION_BITS
    records: dict[str, dict[str, object]] = {}
    digest = hashlib.sha256()

    def certify(
        label: str,
        start: Fraction,
        stop: Fraction,
        evaluator: Callable[[arb], arb],
        maximum_depth: int = 50,
    ) -> None:
        stack = [(Fraction(start), Fraction(stop), 0)]
        leaves = 0
        depth_seen = 0
        while stack:
            lower, upper, depth = stack.pop()
            if evaluator(_arb_interval(lower, upper)) > 0:
                leaves += 1
                depth_seen = max(depth_seen, depth)
                digest.update(
                    (
                        f"{label}|{lower.numerator}/{lower.denominator}|"
                        f"{upper.numerator}/{upper.denominator}\n"
                    ).encode("ascii")
                )
                continue
            if depth >= maximum_depth:
                raise AssertionError(
                    f"n=7 scalar subdivision did not close: {label}, [{lower},{upper}]"
                )
            midpoint = (lower + upper) / 2
            stack.append((lower, midpoint, depth + 1))
            stack.append((midpoint, upper, depth + 1))
        records[label] = {
            "interval": f"[{start},{stop}]",
            "terminal_intervals": leaves,
            "maximum_bisection_depth": depth_seen,
        }

    try:
        core_ends = {
            3: Fraction(9, 4),
            4: Fraction(27, 16),
            5: Fraction(173, 128),
            6: Fraction(9, 8),
        }
        minorants = {
            3: (Fraction(7, 8), Fraction(35, 8), Fraction(1, 25), Fraction(7, 4)),
            4: (Fraction(4, 5), Fraction(19, 5), Fraction(9, 100), Fraction(7, 5)),
            5: (Fraction(47, 64), Fraction(213, 64), Fraction(1, 5), Fraction(7, 6)),
            6: (Fraction(91, 128), Fraction(175, 64), Fraction(481, 1000), Fraction(1)),
        }
        for order, stop in core_ends.items():
            certify(
                f"q{order}-nonnegative-core",
                Fraction(1),
                stop,
                lambda value, order=order: _q_arb(value, order, True),
            )
        for order, (start, stop, slope, center) in minorants.items():
            certify(
                f"q{order}-minorant-lower-piece",
                start,
                Fraction(1),
                lambda value, order=order, slope=slope, center=center: (
                    _q_arb(value, order, False)
                    - _arb_exact(slope) * (_arb_exact(center) - value)
                ),
            )
            certify(
                f"q{order}-minorant-upper-piece",
                Fraction(1),
                stop,
                lambda value, order=order, slope=slope, center=center: (
                    _q_arb(value, order, True)
                    - _arb_exact(slope) * (_arb_exact(center) - value)
                ),
            )

        derivative_bounds = {
            3: Fraction(39, 100),
            4: Fraction(16, 25),
            5: Fraction(49, 40),
        }
        caps = {
            3: Fraction(81, 500),
            4: Fraction(223, 1000),
            5: Fraction(499, 1600),
        }
        for order, bound in derivative_bounds.items():
            for suffix, upper_piece, start, stop in (
                ("lower", False, SMALL_U_CUTOFF, Fraction(1)),
                ("upper", True, Fraction(1), SCALAR_DOMAIN_END),
            ):
                certify(
                    f"q{order}-derivative-{suffix}-positive-side",
                    start,
                    stop,
                    lambda value, order=order, bound=bound, upper_piece=upper_piece: (
                        _arb_exact(bound) + _q_derivative_arb(value, order, upper_piece)
                    ),
                )
                certify(
                    f"q{order}-derivative-{suffix}-negative-side",
                    start,
                    stop,
                    lambda value, order=order, bound=bound, upper_piece=upper_piece: (
                        _arb_exact(bound) - _q_derivative_arb(value, order, upper_piece)
                    ),
                )
        for order, cap in caps.items():
            certify(
                f"q{order}-cap-lower-piece",
                SMALL_U_CUTOFF,
                Fraction(1),
                lambda value, order=order, cap=cap: (
                    _arb_exact(cap) - _q_arb(value, order, False)
                ),
            )
            certify(
                f"q{order}-cap-upper-piece",
                Fraction(1),
                SCALAR_DOMAIN_END,
                lambda value, order=order, cap=cap: (
                    _arb_exact(cap) - _q_arb(value, order, True)
                ),
            )

        for order, base, events in (
            (4, Fraction(1, 500), FACE6_Q4_EVENTS),
            (5, Fraction(1, 600), FACE7_Q5_EVENTS),
        ):
            breakpoints = sorted(
                {
                    SMALL_U_CUTOFF,
                    Fraction(1),
                    SCALAR_DOMAIN_END,
                    *(endpoint for event in events for endpoint in event[:2]),
                }
            )
            for index, (start, stop) in enumerate(itertools.pairwise(breakpoints)):
                midpoint = (start + stop) / 2
                bound = base + sum(
                    step for lower, upper, step in events if lower <= midpoint <= upper
                )
                upper_piece = start >= 1
                certify(
                    f"q{order}-localized-{index}-positive-side",
                    start,
                    stop,
                    lambda value, order=order, bound=bound, upper_piece=upper_piece: (
                        _arb_exact(bound) + _q_derivative_arb(value, order, upper_piece)
                    ),
                )
                certify(
                    f"q{order}-localized-{index}-negative-side",
                    start,
                    stop,
                    lambda value, order=order, bound=bound, upper_piece=upper_piece: (
                        _arb_exact(bound) - _q_derivative_arb(value, order, upper_piece)
                    ),
                )

        # For 0<u<=1/8, each lower-piece q_k is increasing.  Its derivative
        # is also increasing because lambda=7/u>=56>k+2.  Directed checks at
        # the right endpoint therefore close the omitted half-open tail.
        small_u_records = {}
        for order, derivative_bound in {
            3: Fraction(39, 100),
            4: Fraction(1, 500),
            5: Fraction(1, 600),
        }.items():
            derivative_margin = _arb_exact(derivative_bound) - _q_derivative_arb(
                _arb_exact(SMALL_U_CUTOFF), order, False
            )
            cap_margin = _arb_exact(caps[order]) - _q_arb(
                _arb_exact(SMALL_U_CUTOFF), order, False
            )
            assert derivative_margin > 0
            assert cap_margin > 0
            small_u_records[f"q{order}"] = {
                "interval": "(0,1/8]",
                "derivative_bound": str(derivative_bound),
                "cap": str(caps[order]),
                "directed_endpoint_checks": "passed",
            }
    finally:
        ctx.prec = previous_precision

    return {
        "arb_precision_bits": ARB_PRECISION_BITS,
        "symbolic_integrand_and_derivative_identities": "passed",
        "proved_bounds": {
            "q3": {
                "nonnegative_core": "0<=u<=9/4",
                "affine_minorant": "q3(u)>=(1/25)(7/4-u) on [7/8,35/8]",
                "global_derivative": "abs(q3'(u))<39/100 on [0,7]",
                "pointwise_cap": "q3(u)<81/500 on [0,7]",
            },
            "q4": {
                "nonnegative_core": "0<=u<=27/16",
                "affine_minorant": "q4(u)>=(9/100)(7/5-u) on [4/5,19/5]",
                "global_derivative": "abs(q4'(u))<16/25 on [0,7]",
                "pointwise_cap": "q4(u)<223/1000 on [0,7]",
            },
            "q5": {
                "nonnegative_core": "0<=u<=173/128",
                "affine_minorant": "q5(u)>=(1/5)(7/6-u) on [47/64,213/64]",
                "global_derivative": "abs(q5'(u))<49/40 on [0,7]",
                "pointwise_cap": "q5(u)<499/1600 on [0,7]",
            },
            "q6": {
                "nonnegative_core": "0<=u<=9/8",
                "affine_minorant": ("q6(u)>=(481/1000)(1-u) on [91/128,175/64]"),
            },
        },
        "compact_interval_records": records,
        "compact_terminal_intervals": sum(
            int(record["terminal_intervals"]) for record in records.values()
        ),
        "compact_maximum_bisection_depth": max(
            int(record["maximum_bisection_depth"]) for record in records.values()
        ),
        "compact_transcript_sha256": digest.hexdigest(),
        "small_u_analytic_tails": small_u_records,
        "q4_localized_events": [tuple(map(str, event)) for event in FACE6_Q4_EVENTS],
        "q5_localized_events": [tuple(map(str, event)) for event in FACE7_Q5_EVENTS],
    }
