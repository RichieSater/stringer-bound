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
PRACTICE = (
    Path(__file__).resolve().parents[2]
    / "audit"
    / "PRACTICE-SAFEGUARD.md"
)
sys.path.insert(0, str(PYTHON_DIR))

from poisson_band_calibration import (  # noqa: E402
    CASES,
    _exact_zero_anchor_boundaries,
    bolshev_probability_numeric,
    certify_case,
    certify_zero_anchor_case,
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
    def test_zero_anchor_mixed_endpoint_enclosures_have_correct_direction(self):
        n = 3
        brackets = (
            (Fraction(1, 5), Fraction(1, 4)),
            (Fraction(1, 2), Fraction(3, 5)),
            (Fraction(4, 5), Fraction(9, 10)),
        )
        eta = Fraction(2)
        lambdas = [(lower + upper) / 2 for lower, upper in brackets]
        p0 = lambdas[0] / n
        actual = tuple(
            min(Fraction(1), p0 + eta * (lambdas[j] / n - p0))
            for j in range(n)
        )
        lower = _exact_zero_anchor_boundaries(
            n, brackets, eta, "lower")
        upper = _exact_zero_anchor_boundaries(
            n, brackets, eta, "upper")
        self.assertTrue(all(lo <= value <= hi for lo, value, hi in zip(
            lower, actual, upper)))
        self.assertEqual(tuple(sorted(lower)), lower)
        self.assertEqual(tuple(sorted(upper)), upper)

    def test_factor_capping_preserves_calibrated_band_boundaries(self):
        for kappa in (Fraction(1), Fraction(5, 4), Fraction(3)):
            for factor in (Fraction(0), Fraction(1, 2), Fraction(1),
                           Fraction(3, 2)):
                self.assertEqual(
                    min(Fraction(1), kappa * min(Fraction(1), factor)),
                    min(Fraction(1), kappa * factor),
                )
        for eta in (Fraction(1), Fraction(3, 2), Fraction(4)):
            for p0 in (Fraction(0), Fraction(1, 5), Fraction(1)):
                for factor in (p0, Fraction(1, 2), Fraction(1),
                               Fraction(3, 2)):
                    if factor < p0:
                        continue
                    self.assertEqual(
                        min(Fraction(1), p0 + eta *
                            (min(Fraction(1), factor) - p0)),
                        min(Fraction(1), p0 + eta * (factor - p0)),
                    )

    def test_exact_report_uses_nonnegative_factor_form(self):
        n = 5
        brackets = exact_poisson_lambda_brackets("0.05", n, 40)
        taints = (Fraction(1), Fraction(2, 5), Fraction(1, 10))
        kappa = Fraction(5, 4)
        eta = Fraction(3, 2)
        report = exact_calibrated_report(
            n, brackets, kappa, taints, eta)
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
        factorwise_calibrated = _record_fraction(
            report["factorwise_capped_calibrated_poisson_upper"])
        self.assertEqual(
            report["factor_convention"],
            "untruncated_poisson_factors_then_final_cap",
        )
        self.assertEqual(
            report["reported_calibration_variants"],
            ["final_cap", "calibrated_factorwise_cap"],
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
        full_factorwise = [min(Fraction(1), kappa * factor)
                           for factor in factors]
        full_factorwise_direct = (
            full_factorwise[0] * (1 - taints[0])
            + full_factorwise[1] * (taints[0] - taints[1])
            + full_factorwise[2] * (taints[1] - taints[2])
            + full_factorwise[3] * taints[2]
        )
        self.assertEqual(factorwise_calibrated, full_factorwise_direct)
        self.assertLessEqual(factorwise_calibrated, calibrated)
        coefficient_p0 = 1 - eta * taints[0]
        anchored_direct = coefficient_p0 * brackets[0][0] / n
        anchored_direct += eta * (
            factors[1] * (taints[0] - taints[1])
            + factors[2] * (taints[1] - taints[2])
            + factors[3] * taints[2]
        )
        self.assertEqual(
            _record_fraction(
                report["zero_anchor_calibrated_poisson_upper"]),
            min(Fraction(1), anchored_direct),
        )
        p0_upper = brackets[0][1] / n
        p0_lower = brackets[0][0] / n
        anchored_factorwise = [p0_upper]
        anchored_factorwise.extend(
            min(Fraction(1), eta * factors[j] - (eta - 1) * p0_lower)
            for j in range(1, 4)
        )
        anchored_factorwise_direct = (
            anchored_factorwise[0] * (1 - taints[0])
            + anchored_factorwise[1] * (taints[0] - taints[1])
            + anchored_factorwise[2] * (taints[1] - taints[2])
            + anchored_factorwise[3] * taints[2]
        )
        self.assertEqual(
            _record_fraction(
                report["zero_anchor_factorwise_capped_poisson_upper"]),
            anchored_factorwise_direct,
        )
        self.assertLessEqual(
            anchored_factorwise_direct,
            _record_fraction(
                report["zero_anchor_calibrated_poisson_upper"]),
        )
        self.assertIn(
            "selected before observing",
            report["calibration_selection_warning"],
        )

    def test_zero_anchor_report_leaves_no_error_factor_unscaled(self):
        n = 25
        brackets = exact_poisson_lambda_brackets("0.05", 1, 40)
        report = exact_calibrated_report(
            n, brackets, Fraction(5, 4), (), Fraction(2))
        self.assertEqual(
            _record_fraction(
                report["zero_anchor_calibrated_poisson_upper"]),
            brackets[0][1] / n,
        )
        self.assertEqual(
            _record_fraction(
                report["zero_anchor_factorwise_capped_poisson_upper"]),
            brackets[0][1] / n,
        )

    def test_calibrated_factor_capping_can_strictly_reduce_report(self):
        n = 25
        brackets = exact_poisson_lambda_brackets("0.05", 20, 40)
        report = exact_calibrated_report(
            n,
            brackets,
            Fraction(5, 4),
            (Fraction(1, 10),) * 20,
            Fraction(2),
        )
        self.assertTrue(report["factorwise_capped_upper_is_lower"])
        self.assertTrue(
            report["zero_anchor_factorwise_capped_upper_is_lower"])
        self.assertLess(
            _record_fraction(
                report["factorwise_capped_calibrated_poisson_upper"]),
            _record_fraction(report["calibrated_poisson_upper"]),
        )
        self.assertLess(
            _record_fraction(
                report["zero_anchor_factorwise_capped_poisson_upper"]),
            _record_fraction(
                report["zero_anchor_calibrated_poisson_upper"]),
        )

    def test_anchor_factor_enclosure_preserves_pointwise_ordering(self):
        # Independent upper endpoints for the capped effective factors can
        # be microscopically larger than the sign-aware affine enclosure.
        # Intersecting the two valid enclosures preserves the theorem's
        # pointwise ordering at the reporting layer.
        brackets = exact_poisson_lambda_brackets("0.05", 5, 32)
        report = exact_calibrated_report(
            5,
            brackets,
            Fraction(5, 4),
            (Fraction(22, 25),),
            Fraction(11, 10),
        )
        self.assertEqual(
            _record_fraction(
                report["zero_anchor_factorwise_capped_poisson_upper"]),
            _record_fraction(
                report["zero_anchor_calibrated_poisson_upper"]),
        )

    def test_zero_anchor_trivial_capped_case(self):
        brackets = exact_poisson_lambda_brackets("0.05", 0, 40)
        case = certify_zero_anchor_case(1, Fraction("0.05"), brackets)
        self.assertEqual(_record_fraction(case["eta_lower"]), 1)
        self.assertEqual(_record_fraction(case["eta_upper"]), 1)
        report = exact_calibrated_report(
            1, brackets, Fraction(1), (), Fraction(1))
        self.assertEqual(
            _record_fraction(
                report["zero_anchor_calibrated_poisson_upper"]),
            1,
        )

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
        anchor_case = certify_zero_anchor_case(
            8, Fraction("0.10"), brackets)
        self.assertEqual(_record_fraction(anchor_case["eta_lower"]), 1)
        self.assertEqual(_record_fraction(anchor_case["eta_upper"]), 1)

    def test_certificate_routine_rejects_out_of_scope_tail_level(self):
        brackets = exact_poisson_lambda_brackets("0.50", 1, 32)
        with self.assertRaises(ValueError):
            certify_case(2, Fraction(1, 2), brackets)
        with self.assertRaises(ValueError):
            certify_zero_anchor_case(2, Fraction(1, 2), brackets)

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
        anchor_case = certify_zero_anchor_case(
            9, Fraction("0.10"), brackets)
        anchor_lower = _record_fraction(anchor_case["eta_lower"])
        anchor_upper = _record_fraction(anchor_case["eta_upper"])
        self.assertEqual(
            anchor_upper - anchor_lower, Fraction(1, 1 << 28))
        self.assertGreater(anchor_lower, 1)

    def test_committed_certificate_has_opposite_exact_signs(self):
        payload = json.loads(CERTIFICATE.read_text())
        self.assertEqual(payload["schema_version"], 2)
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
            self.assertEqual(
                [case["n"] for case in level["zero_anchor_cases"]],
                level["sample_sizes"],
            )
            kappa_by_n = {
                case["n"]: _record_fraction(case["kappa_upper"])
                for case in level["cases"]
            }
            for case in level["zero_anchor_cases"]:
                lower = _record_fraction(case["eta_lower"])
                upper = _record_fraction(case["eta_upper"])
                self.assertGreaterEqual(upper, kappa_by_n[case["n"]])
                self.assertEqual(upper - lower, Fraction(1, 1 << 28))
                self.assertEqual(
                    case["eta_upper"]["valid_decimal_ceiling_12"],
                    _upward_decimal(upper, 12),
                )
                self.assertLess(
                    _record_fraction(
                        case["event_probability_upper_at_eta_lower"]),
                    nominal,
                )
                self.assertGreaterEqual(
                    _record_fraction(
                        case["event_probability_lower_at_eta_upper"]),
                    nominal,
                )

    def test_manuscript_table_rounds_every_certified_upper_upward(self):
        payload = json.loads(CERTIFICATE.read_text())
        paper = PAPER.read_text()
        theory = THEORY.read_text()
        practice = PRACTICE.read_text()
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
            anchored_displayed = [
                _upward_decimal(
                    _record_fraction(case["eta_upper"]), 6)
                for case in level["zero_anchor_cases"]
            ]
            anchored_row = "| %s%% | %s |" % (
                confidence[level["alpha"]],
                " | ".join(anchored_displayed),
            )
            self.assertIn(anchored_row, practice)
            for case in level["zero_anchor_cases"]:
                displayed_12 = _upward_decimal(
                    _record_fraction(case["eta_upper"]), 12)
                pattern = r"\|\s*%s%%\s*\|\s*%d\s*\|\s*%s\s*\|" % (
                    confidence[level["alpha"]], case["n"],
                    re.escape(displayed_12))
                self.assertRegex(theory, pattern)


if __name__ == "__main__":
    unittest.main()
