"""Generate a reproducible, decision-oriented audit-method benchmark.

The benchmark compares six upper-taint-rate calculations on fixed,
zero-heavy sample profiles at practical sample sizes:

* ordinary binomial and Poisson Stringer (reference calculations whose
  general coverage is unresolved at these sample sizes);
* the finite-sample-valid Gaffke bounded-mean endpoint;
* the pre-specified maximum of Poisson Stringer and Gaffke;
* the proved full-scale calibrated Poisson rule; and
* the proved zero-taint-preserving calibrated Poisson rule.

Every displayed rate is an exact rational upper enclosure.  Gaffke brackets
have exact tail-sign checks; Stringer factors have exact endpoint signs; and
the calibration multipliers and Poisson-factor brackets are read from the
committed exact certificate.  The output translates rate differences to
basis points and to an illustrative amount per $10 million solely to expose
decision scale.  It is not a coverage simulation, a sampling-design bridge,
or professional guidance.

Usage::

    python audit_decision_benchmark.py
    python audit_decision_benchmark.py --out audit-decision-benchmark.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from gaffke import gaffke_quantile_certificate
from poisson_band_calibration import exact_calibrated_report
from stringer import exact_binomial_factor_brackets


HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[2]
CERTIFICATE_DIR = HERE.parent / "certificates"
CALIBRATION_PATH = CERTIFICATE_DIR / "poisson-band-calibration-certificate.json"
OUTPUT_NAME = "audit-decision-benchmark.json"
FACTOR_BITS = 96
GAFFKE_BITS = 64
ILLUSTRATIVE_RECORDED_AMOUNT = 10_000_000

# Zero taints are omitted but counted through n.  All profiles are fixed in
# source so regeneration cannot silently select favorable examples.
CASES = (
    ("all-zero-n25-95", 25, "0.05", ()),
    ("mixed-three-n25-95", 25, "0.05", ("1", "0.4", "0.1")),
    ("all-zero-n100-95", 100, "0.05", ()),
    ("one-small-n100-95", 100, "0.05", ("0.1",)),
    ("one-full-n100-95", 100, "0.05", ("1",)),
    ("mixed-three-n100-95", 100, "0.05", ("1", "0.4", "0.1")),
    (
        "six-small-n100-95",
        100,
        "0.05",
        ("0.2", "0.1", "0.05", "0.02", "0.01", "0.005"),
    ),
    ("mixed-three-n200-95", 200, "0.05", ("1", "0.4", "0.1")),
    ("mixed-three-n100-90", 100, "0.10", ("1", "0.4", "0.1")),
    ("mixed-three-n100-99", 100, "0.01", ("1", "0.4", "0.1")),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_from_record(record) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def _ceil_div(numerator: int, denominator: int) -> int:
    return -(-numerator // denominator)


def _fixed_decimal(integer: int, places: int) -> str:
    if places == 0:
        return str(integer)
    scale = 10 ** places
    whole, fractional = divmod(integer, scale)
    return f"{whole}.{fractional:0{places}d}"


def _decimal_ceiling(value: Fraction, places: int) -> str:
    scaled = _ceil_div(value.numerator * (10 ** places), value.denominator)
    return _fixed_decimal(scaled, places)


def _decimal_floor(value: Fraction, places: int) -> str:
    scaled = value.numerator * (10 ** places) // value.denominator
    return _fixed_decimal(scaled, places)


def _fraction_record(value: Fraction, places: int = 12):
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        f"decimal_ceiling_{places}": _decimal_ceiling(value, places),
    }


def _interval_record(lower: Fraction, upper: Fraction, places: int = 12):
    lower, upper = Fraction(lower), Fraction(upper)
    if lower > upper:
        raise ValueError("interval endpoints are reversed")
    return {
        "lower": {
            "numerator": str(lower.numerator),
            "denominator": str(lower.denominator),
            f"decimal_floor_{places}": _decimal_floor(lower, places),
        },
        "upper": _fraction_record(upper, places),
        "width": _fraction_record(upper - lower, places),
    }


def _money_ceiling_record(rate: Fraction):
    cents = _ceil_div(
        rate.numerator * ILLUSTRATIVE_RECORDED_AMOUNT * 100,
        rate.denominator,
    )
    return {
        "exact_ceiling_cents": str(cents),
        "decimal_dollars": _fixed_decimal(cents, 2),
    }


def _rate_record(rate: Fraction):
    return {
        "rate_upper": _fraction_record(rate),
        "basis_points_upper": _decimal_ceiling(rate * 10_000, 6),
        "amount_upper_per_10m": _money_ceiling_record(rate),
    }


def _stringer_interval(factor_brackets, n: int, taints, *, poisson=False):
    """Enclose Stringer via its nonnegative summation-by-parts form."""
    ordered = sorted((Fraction(value) for value in taints if Fraction(value)),
                     reverse=True)
    if any(value < 0 or value > 1 for value in ordered):
        raise ValueError("taints must lie in [0,1]")
    if len(ordered) > n:
        raise ValueError("more taints than sample units")
    if len(factor_brackets) <= len(ordered):
        raise ValueError("factor prefix is too short")

    if ordered:
        coefficients = [1 - ordered[0]]
        coefficients.extend(
            ordered[j - 1] - ordered[j] for j in range(1, len(ordered)))
        coefficients.append(ordered[-1])
    else:
        coefficients = [Fraction(1)]
    if any(coefficient < 0 for coefficient in coefficients):
        raise AssertionError("summation-by-parts coefficient is negative")
    divisor = n if poisson else 1
    lower = sum(
        coefficient * factor_brackets[j][0] / divisor
        for j, coefficient in enumerate(coefficients)
    )
    upper = sum(
        coefficient * factor_brackets[j][1] / divisor
        for j, coefficient in enumerate(coefficients)
    )
    return lower, upper


def _calibration_level(data, alpha: str):
    return next(level for level in data["levels"] if level["alpha"] == alpha)


def _lambda_brackets(level):
    output = []
    for record in level["poisson_lambda_brackets"]:
        denominator = int(record["dyadic_denominator"])
        output.append((
            Fraction(int(record["lower_numerator"]), denominator),
            Fraction(int(record["upper_numerator"]), denominator),
        ))
    return tuple(output)


def _case_multiplier(level, n: int, key: str, field: str) -> Fraction:
    case = next(item for item in level[key] if item["n"] == n)
    return _fraction_from_record(case[field])


def _switch_band(left_name: str, left: Fraction,
                 right_name: str, right: Fraction):
    """Describe thresholds at which two conservative upper endpoints differ."""
    if left <= right:
        lower_name, lower, upper_name, upper = left_name, left, right_name, right
    else:
        lower_name, lower, upper_name, upper = right_name, right, left_name, left
    width = upper - lower
    return {
        "lower_endpoint_method": lower_name,
        "upper_endpoint_method": upper_name,
        "threshold_interval": {
            "lower_inclusive": _fraction_record(lower),
            "upper_exclusive": _fraction_record(upper),
        },
        "has_nonempty_switch_band": width > 0,
        "width_rate": _fraction_record(width),
        "width_basis_points_upper": _decimal_ceiling(width * 10_000, 6),
        "width_amount_per_10m_upper": _money_ceiling_record(width),
        "interpretation": (
            "For a rate threshold in [lower,upper), the method at the lower "
            "endpoint is at or below the threshold while the method at the "
            "upper endpoint is above it, using the certified upper "
            "enclosures shown."
        ),
    }


def _build_case(case, calibration_data):
    case_id, n, alpha, taint_text = case
    taints = tuple(Fraction(value) for value in taint_text)
    level = _calibration_level(calibration_data, alpha)
    poisson_brackets = _lambda_brackets(level)
    kappa = _case_multiplier(level, n, "cases", "kappa_upper")
    eta = _case_multiplier(level, n, "zero_anchor_cases", "eta_upper")

    binomial_brackets = exact_binomial_factor_brackets(
        n, alpha, FACTOR_BITS)
    binomial_lower, binomial_upper = _stringer_interval(
        binomial_brackets, n, taints)
    poisson_lower, poisson_upper = _stringer_interval(
        poisson_brackets, n, taints, poisson=True)

    gaffke = gaffke_quantile_certificate(
        taints, n, alpha, grid_bits=GAFFKE_BITS)
    if not (gaffke.tail_at_lower >= Fraction(alpha)
            >= gaffke.tail_at_upper):
        raise AssertionError("Gaffke tail signs do not bracket alpha")

    calibrated = exact_calibrated_report(
        n, poisson_brackets, kappa, taints, eta)
    ordinary_from_calibration = _fraction_from_record(
        calibrated["ordinary_poisson_stringer_upper"])
    if ordinary_from_calibration != poisson_upper:
        raise AssertionError("ordinary Poisson paths disagree")
    full_upper = _fraction_from_record(
        calibrated["factorwise_capped_calibrated_poisson_upper"])
    anchor_upper = _fraction_from_record(
        calibrated["zero_anchor_factorwise_capped_poisson_upper"])

    safeguard_lower = max(poisson_lower, gaffke.lower)
    safeguard_upper = max(poisson_upper, gaffke.upper)
    if poisson_lower >= gaffke.upper:
        governing = "ordinary_poisson_stringer"
    elif gaffke.lower >= poisson_upper:
        governing = "gaffke"
    else:
        governing = "not_resolved_by_brackets"

    rates = {
        "ordinary_binomial_stringer_reference": binomial_upper,
        "ordinary_poisson_stringer_reference": poisson_upper,
        "gaffke_valid_endpoint": gaffke.upper,
        "poisson_stringer_gaffke_max": safeguard_upper,
        "full_scale_calibrated_factorwise_cap": full_upper,
        "zero_anchor_calibrated_factorwise_cap": anchor_upper,
    }
    methods = {
        "ordinary_binomial_stringer_reference": {
            "mathematical_status": (
                "ordinary coverage unresolved at this sample size; exact "
                "reference calculation only"
            ),
            "rate_interval": _interval_record(
                binomial_lower, binomial_upper),
            **_rate_record(binomial_upper),
        },
        "ordinary_poisson_stringer_reference": {
            "mathematical_status": (
                "ordinary coverage unresolved at this sample size; exact "
                "reference calculation only"
            ),
            "rate_interval": _interval_record(poisson_lower, poisson_upper),
            **_rate_record(poisson_upper),
        },
        "gaffke_valid_endpoint": {
            "mathematical_status": (
                "proved finite-sample valid under the independent bounded-"
                "mean model, using the cited external validity theorem"
            ),
            "rate_interval": _interval_record(gaffke.lower, gaffke.upper),
            "tail_at_lower_at_least_alpha": True,
            "tail_at_upper_at_most_alpha": True,
            "certificate_bits": GAFFKE_BITS,
            **_rate_record(gaffke.upper),
        },
        "poisson_stringer_gaffke_max": {
            "mathematical_status": (
                "proved finite-sample valid under the independent bounded-"
                "mean model"
            ),
            "rate_interval": _interval_record(
                safeguard_lower, safeguard_upper),
            "governing_component_certified": governing,
            **_rate_record(safeguard_upper),
        },
        "full_scale_calibrated_factorwise_cap": {
            "mathematical_status": (
                "proved finite-sample valid modified Poisson rule under the "
                "independent bounded-mean model"
            ),
            "kappa_upper": _fraction_record(kappa),
            **_rate_record(full_upper),
        },
        "zero_anchor_calibrated_factorwise_cap": {
            "mathematical_status": (
                "proved finite-sample valid modified Poisson rule under the "
                "independent bounded-mean model"
            ),
            "eta_upper": _fraction_record(eta),
            **_rate_record(anchor_upper),
        },
    }

    comparisons = {
        "gaffke_vs_ordinary_poisson": _switch_band(
            "gaffke_valid_endpoint", rates["gaffke_valid_endpoint"],
            "ordinary_poisson_stringer_reference",
            rates["ordinary_poisson_stringer_reference"],
        ),
        "full_scale_vs_ordinary_poisson": _switch_band(
            "full_scale_calibrated_factorwise_cap",
            rates["full_scale_calibrated_factorwise_cap"],
            "ordinary_poisson_stringer_reference",
            rates["ordinary_poisson_stringer_reference"],
        ),
        "zero_anchor_vs_ordinary_poisson": _switch_band(
            "zero_anchor_calibrated_factorwise_cap",
            rates["zero_anchor_calibrated_factorwise_cap"],
            "ordinary_poisson_stringer_reference",
            rates["ordinary_poisson_stringer_reference"],
        ),
        "full_scale_vs_zero_anchor": _switch_band(
            "full_scale_calibrated_factorwise_cap",
            rates["full_scale_calibrated_factorwise_cap"],
            "zero_anchor_calibrated_factorwise_cap",
            rates["zero_anchor_calibrated_factorwise_cap"],
        ),
    }
    return {
        "id": case_id,
        "sample_size": n,
        "alpha": alpha,
        "nominal_confidence": str(1 - Fraction(alpha)),
        "nonzero_taints_descending": [str(value) for value in sorted(
            taints, reverse=True)],
        "nonzero_taint_count": len(taints),
        "zero_taint_count": n - len(taints),
        "methods": methods,
        "threshold_switch_bands": comparisons,
    }


def build_benchmark():
    calibration_raw = CALIBRATION_PATH.read_bytes()
    calibration_data = json.loads(calibration_raw)
    cases = [_build_case(case, calibration_data) for case in CASES]
    return {
        "schema_version": 1,
        "generated_by": (
            "supporting-materials/computations/python/"
            "audit_decision_benchmark.py"
        ),
        "claim": (
            "A fixed suite of practical-size, zero-heavy samples has exact "
            "method-comparison outputs and exact threshold-switch bands."
        ),
        "mathematical_status": (
            "Reproducible exact benchmark. It illustrates the output scale "
            "of proved valid alternatives; it does not establish coverage "
            "of ordinary Stringer at these sample sizes."
        ),
        "model_scope": (
            "The validity labels assume independent [0,1]-valued "
            "observations with a common mean. No without-replacement, "
            "systematic-PPS, stratification, or audit-unit-cap bridge is "
            "asserted."
        ),
        "decision_semantics": (
            "For a pre-specified upper-rate method and a separately given "
            "rate threshold tau, the numerical comparison is bound<=tau. "
            "The benchmark does not choose tau or turn that comparison into "
            "an audit conclusion."
        ),
        "method_selection_semantics": (
            "The full-scale and zero-anchor calibration paths are separate "
            "proved procedures. Selecting the smaller path after observing "
            "the sample is not certified."
        ),
        "illustrative_recorded_amount_dollars": str(
            ILLUSTRATIVE_RECORDED_AMOUNT),
        "factor_bits": FACTOR_BITS,
        "gaffke_certificate_bits": GAFFKE_BITS,
        "source_sha256": {
            str(CALIBRATION_PATH.relative_to(REPOSITORY_ROOT)):
                hashlib.sha256(calibration_raw).hexdigest(),
            str((HERE / "gaffke.py").relative_to(REPOSITORY_ROOT)):
                _sha256(HERE / "gaffke.py"),
            str((HERE / "poisson_band_calibration.py").relative_to(
                REPOSITORY_ROOT)):
                _sha256(HERE / "poisson_band_calibration.py"),
            str((HERE / "stringer.py").relative_to(REPOSITORY_ROOT)):
                _sha256(HERE / "stringer.py"),
        },
        "case_count": len(cases),
        "cases": cases,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        help="write the machine-readable benchmark")
    args = parser.parse_args(argv)
    benchmark = build_benchmark()
    rendered = json.dumps(benchmark, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    print(
        f"audit decision benchmark: {benchmark['case_count']} exact cases; "
        "all exact bracket and relationship checks passed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
