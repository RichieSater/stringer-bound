"""Tests for the all-n scalar Poisson-band calibration."""

from __future__ import annotations

import json
import re
import sys
import unittest
from fractions import Fraction
from pathlib import Path

from mpmath import mp


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    Path(__file__).resolve().parents[1]
    / "certificates"
    / "poisson-band-calibration-certificate.json"
)
PAPER = Path(__file__).resolve().parents[2] / "paper" / "stringer.tex"
THEORY = Path(__file__).resolve().parents[2] / "theory" / "POISSON-BAND-CALIBRATION.md"
sys.path.insert(0, str(PYTHON_DIR))

from poisson_band_calibration import (  # noqa: E402
    CASES,
    bolshev_probability_numeric,
    certify_case,
    exact_calibrated_report,
)
from poisson_band_certificate import bolshev_probability  # noqa: E402
from stringer import exact_poisson_lambda_brackets  # noqa: E402


def _record_fraction(record):
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def _upward_decimal(value, digits):
    value = Fraction(value)
    scale = 10 ** digits
    numerator = (value.numerator * scale + value.denominator - 1) // value.denominator
    return f"{numerator // scale}.{numerator % scale:0{digits}d}"


class ScalarCalibrationTests(unittest.TestCase):
    def test_factor_capping_preserves_calibrated_band_boundaries(self):
        for kappa in (Fraction(1), Fraction(5, 4), Fraction(3)):
            for factor in (Fraction(0), Fraction(1, 2), Fraction(1),
                           Fraction(3, 2)):
                self.assertEqual(
                    min(Fraction(1), kappa * min(Fraction(1), factor)),
                    min(Fraction(1), kappa * factor),
                )

    def test_exact_report_uses_nonnegative_factor_form(self):
        n = 5
        brackets = exact_poisson_lambda_brackets("0.05", n, 40)
        taints = (Fraction(1), Fraction(2, 5), Fraction(1, 10))
        kappa = Fraction(5, 4)
        report = exact_calibrated_report(n, brackets, kappa, taints)
        factors = [upper / n for _, upper in brackets]
        direct = (
            factors[0]
            + (factors[1] - factors[0]) * taints[0]
            + (factors[2] - factors[1]) * taints[1]
            + (factors[3] - factors[2]) * taints[2]
        )
        ordinary = _record_fraction(
            report["ordinary_poisson_stringer_upper"])
        calibrated = _record_fraction(report["calibrated_poisson_upper"])
        self.assertEqual(
            report["factor_convention"],
            "untruncated_poisson_factors_then_final_cap",
        )
        self.assertEqual(report["sample_size"], n)
        self.assertEqual(report["nonzero_taint_count"], 3)
        self.assertEqual(report["zero_taint_count"], 2)
        self.assertEqual(
            [_record_fraction(record)
             for record in report["sorted_nonzero_taints"]],
            list(taints),
        )
        self.assertEqual(ordinary, direct)
        self.assertEqual(calibrated, min(Fraction(1), kappa * direct))

    def test_exact_report_validates_taints(self):
        brackets = exact_poisson_lambda_brackets("0.05", 1, 32)
        with self.assertRaises(ValueError):
            exact_calibrated_report(
                1, brackets, Fraction(1), (Fraction(6, 5),))

    def test_numeric_bolshev_locator_matches_exact_special_case(self):
        boundaries = [Fraction(1, 4), Fraction(2, 3), Fraction(1)]
        exact = bolshev_probability(boundaries)
        with mp.workdps(80):
            numeric = bolshev_probability_numeric([
                mp.mpf(value.numerator) / value.denominator
                for value in boundaries
            ])
            self.assertAlmostEqual(
                float(numeric), float(exact), places=14)

    def test_unadjusted_certified_frontier_returns_one(self):
        brackets = exact_poisson_lambda_brackets("0.10", 7, 48)
        case = certify_case(8, Fraction("0.10"), brackets)
        self.assertEqual(_record_fraction(case["kappa_lower"]), 1)
        self.assertEqual(_record_fraction(case["kappa_upper"]), 1)
        self.assertGreaterEqual(
            _record_fraction(
                case["event_probability_lower_at_kappa_upper"]),
            Fraction(9, 10),
        )

    def test_certificate_routine_rejects_out_of_scope_tail_level(self):
        brackets = exact_poisson_lambda_brackets("0.50", 1, 32)
        with self.assertRaises(ValueError):
            certify_case(2, Fraction(1, 2), brackets)

    def test_first_post_frontier_case_has_adjacent_dyadic_bracket(self):
        brackets = exact_poisson_lambda_brackets("0.10", 8, 48)
        case = certify_case(9, Fraction("0.10"), brackets)
        lower = _record_fraction(case["kappa_lower"])
        upper = _record_fraction(case["kappa_upper"])
        self.assertEqual(upper - lower, Fraction(1, 1 << 28))
        self.assertLess(
            _record_fraction(
                case["event_probability_upper_at_kappa_lower"]),
            Fraction(9, 10),
        )
        self.assertGreaterEqual(
            _record_fraction(
                case["event_probability_lower_at_kappa_upper"]),
            Fraction(9, 10),
        )

    def test_committed_certificate_has_opposite_exact_signs(self):
        payload = json.loads(CERTIFICATE.read_text())
        self.assertEqual(payload["schema_version"], 1)
        self.assertGreaterEqual(payload["exponential_series_pairs"], 1)
        self.assertIn(
            "kappa_upper.valid_decimal_ceiling_12",
            payload["decimal_policy"],
        )
        self.assertEqual(
            {
                level["alpha"]: tuple(level["sample_sizes"])
                for level in payload["levels"]
            },
            CASES,
        )
        for level in payload["levels"]:
            nominal = 1 - Fraction(level["alpha"])
            factor_records = level["poisson_lambda_brackets"]
            self.assertEqual(len(factor_records), max(level["sample_sizes"]))
            factor_denominator = 1 << payload["factor_bits"]
            previous_lower = None
            for j, record in enumerate(factor_records):
                self.assertEqual(record["j"], j)
                self.assertEqual(
                    int(record["dyadic_denominator"]), factor_denominator)
                lower_numerator = int(record["lower_numerator"])
                upper_numerator = int(record["upper_numerator"])
                self.assertEqual(upper_numerator - lower_numerator, 1)
                if previous_lower is not None:
                    self.assertLess(previous_lower, lower_numerator)
                previous_lower = lower_numerator
                self.assertEqual(
                    record["endpoint_signs"],
                    {"cdf_at_lower": "> alpha", "cdf_at_upper": "< alpha"},
                )
            for case in level["cases"]:
                lower = _record_fraction(case["kappa_lower"])
                upper = _record_fraction(case["kappa_upper"])
                self.assertEqual(upper - lower, Fraction(1, 1 << 28))
                self.assertEqual(
                    case["kappa_upper"]["valid_decimal_ceiling_12"],
                    _upward_decimal(upper, 12),
                )
                self.assertLess(
                    _record_fraction(
                        case["event_probability_upper_at_kappa_lower"]),
                    nominal,
                )
                self.assertGreaterEqual(
                    _record_fraction(
                        case["event_probability_lower_at_kappa_upper"]),
                    nominal,
                )

    def test_manuscript_table_rounds_every_certified_upper_upward(self):
        payload = json.loads(CERTIFICATE.read_text())
        paper = PAPER.read_text()
        theory = THEORY.read_text()
        confidence = {"0.10": "90", "0.05": "95", "0.01": "99"}
        for level in payload["levels"]:
            displayed = [
                _upward_decimal(
                    _record_fraction(case["kappa_upper"]), 6)
                for case in level["cases"]
            ]
            row = "$%s\\%%$ & %s\\\\" % (
                confidence[level["alpha"]], " & ".join(f"${x}$" for x in displayed))
            self.assertIn(row, paper)
            for case in level["cases"]:
                displayed_12 = _upward_decimal(
                    _record_fraction(case["kappa_upper"]), 12)
                pattern = r"\|\s*%s%%\s*\|\s*%d\s*\|\s*%s\s*\|" % (
                    confidence[level["alpha"]], case["n"],
                    re.escape(displayed_12))
                self.assertRegex(theory, pattern)


if __name__ == "__main__":
    unittest.main()
