"""Regression tests for the exact binomial-certificate path."""

from __future__ import annotations

import hashlib
import json
import sys
import unittest
from fractions import Fraction
from math import comb
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

import stringer  # noqa: E402
from coverage import coverage_exact  # noqa: E402
from n3_gaffke_certificate import build_certificate  # noqa: E402
from n5_gaffke_certificate import Interval as N5Interval  # noqa: E402
from summarize_certificates import default_paths, summarize  # noqa: E402
from two_point_lemma import check_poisson_dominates  # noqa: E402


def direct_cdf(p: Fraction, n: int, j: int) -> Fraction:
    return sum(
        Fraction(comb(n, i)) * p ** i * (1 - p) ** (n - i)
        for i in range(j + 1)
    )


class ExactFactorTests(unittest.TestCase):
    def test_dyadic_sign_matches_direct_fraction_polynomial(self):
        alpha = Fraction(7, 20)
        bits = 5
        denominator = 1 << bits
        for n in range(1, 7):
            for j in range(n):
                for k in range(denominator + 1):
                    expected_delta = direct_cdf(
                        Fraction(k, denominator), n, j) - alpha
                    expected = ((expected_delta > 0)
                                - (expected_delta < 0))
                    actual = stringer._binomial_cdf_sign_dyadic(
                        k, bits, n, j, alpha)
                    self.assertEqual(actual, expected, (n, j, k))

    def test_exact_brackets_have_certified_endpoint_signs(self):
        n, bits = 12, 36
        alpha = Fraction(13, 20)
        brackets = stringer.exact_binomial_factor_brackets(
            n, str(alpha), bits)
        self.assertEqual(len(brackets), n + 1)
        for j, (lower, upper) in enumerate(brackets[:-1]):
            self.assertLessEqual(upper - lower, Fraction(1, 1 << bits))
            lower_k = lower.numerator << (bits -
                                           (lower.denominator.bit_length() - 1))
            upper_k = upper.numerator << (bits -
                                           (upper.denominator.bit_length() - 1))
            self.assertGreaterEqual(
                stringer._binomial_cdf_sign_dyadic(
                    lower_k, bits, n, j, alpha), 0)
            self.assertLessEqual(
                stringer._binomial_cdf_sign_dyadic(
                    upper_k, bits, n, j, alpha), 0)
        self.assertEqual(brackets[-1], (Fraction(1), Fraction(1)))
        self.assertTrue(all(
            brackets[j][0] <= brackets[j + 1][0]
            for j in range(n)
        ))

    def test_poisson_domination_is_confidence_level_dependent(self):
        self.assertTrue(check_poisson_dominates(2, "0.05", dps=80))
        self.assertFalse(check_poisson_dominates(2, "0.70", dps=80))

    def test_n5_directed_dyadic_intervals_enclose_exact_arithmetic(self):
        values = (
            Fraction(-7, 13), Fraction(-2, 7),
            Fraction(2, 7), Fraction(7, 13), Fraction(4, 9),
        )
        for left in values:
            x = N5Interval(left)
            self.assertLessEqual(x.lower, left)
            self.assertGreaterEqual(x.upper, left)
            if left:
                reciprocal = x.reciprocal()
                self.assertLessEqual(reciprocal.lower, 1 / left)
                self.assertGreaterEqual(reciprocal.upper, 1 / left)
            for right in values:
                y = N5Interval(right)
                for enclosure, exact in (
                    (x + y, left + right),
                    (x - y, left - right),
                    (x * y, left * right),
                ):
                    self.assertLessEqual(enclosure.lower, exact)
                    self.assertGreaterEqual(enclosure.upper, exact)


class ExactCoverageTests(unittest.TestCase):
    def test_one_nonzero_value_matches_direct_binomial_sum(self):
        n = 7
        alpha = "0.05"
        value, probability = Fraction(3, 5), Fraction(2, 7)
        coverage, theta, min_gap = coverage_exact(
            [value], [probability], n, alpha, factor_bits=64)

        numerical_factors = [
            float(p) for p, _ in stringer.factors(n, alpha, dps=80)
        ]
        expected = Fraction(0)
        for k in range(n + 1):
            sb = ((1 - float(value)) * numerical_factors[0]
                  + float(value) * numerical_factors[k])
            if sb >= float(theta):
                expected += (Fraction(comb(n, k)) * probability ** k
                             * (1 - probability) ** (n - k))

        self.assertEqual(coverage, expected)
        self.assertGreaterEqual(coverage, 1 - Fraction(alpha))
        self.assertGreaterEqual(min_gap, 0)

    def test_poisson_is_not_labeled_formally_exact(self):
        with self.assertRaises(ValueError):
            coverage_exact(
                [Fraction(1, 2)], [Fraction(1, 3)], 4, "0.05",
                method="poisson")


class CertificateSummaryTests(unittest.TestCase):
    def test_n3_conventional_level_certificate(self):
        certificate = build_certificate()
        self.assertEqual(
            [level["alpha"] for level in certificate["levels"]],
            ["0.01", "0.05", "0.10"],
        )
        for level in certificate["levels"]:
            self.assertIn("pointwise dominates", level["conclusion"])
            self.assertTrue(all(
                Fraction(int(value["lower"]["numerator"]),
                         int(value["lower"]["denominator"])) > 0
                for value in level["region_c_am_gm_margins"].values()
            ))
            for region in (
                level["region_a"],
                level["region_b_y_equals_1_boundary"],
            ):
                enclosure = region["structural_zero_enclosure"]
                lower = Fraction(int(enclosure["lower"]["numerator"]),
                                 int(enclosure["lower"]["denominator"]))
                upper = Fraction(int(enclosure["upper"]["numerator"]),
                                 int(enclosure["upper"]["denominator"]))
                self.assertLessEqual(lower, 0)
                self.assertGreaterEqual(upper, 0)

    def test_n4_conventional_level_certificate_artifact(self):
        certificate_path = (PYTHON_DIR.parent / "certificates"
                            / "n4-gaffke-certificate.json")
        certificate = json.loads(certificate_path.read_text())
        self.assertEqual(
            [level["alpha"] for level in certificate["levels"]],
            ["0.01", "0.05", "0.10"],
        )
        for level in certificate["levels"]:
            self.assertIn("pointwise dominates", level["conclusion"])
            minimum = level["minimum_positive_bernstein_coefficient"]
            self.assertGreater(
                Fraction(int(minimum["lower"]["numerator"]),
                         int(minimum["lower"]["denominator"])),
                0,
            )
        expected_region_minima = {
            "0.01": ("3.10e-10", "4.78e-14", "5.84e-12"),
            "0.05": ("3.47e-07", "2.66e-10", "4.60e-09"),
            "0.10": ("7.74e-06", "1.18e-08", "7.93e-08"),
        }
        for level in certificate["levels"]:
            observed = []
            for region_name in ("A", "B", "C"):
                tetrahedra = level["polynomial_regions"][region_name][
                    "tetrahedra"]
                minima = [
                    Fraction(
                        int(item["minimum_positive_coefficient"]["lower"]
                            ["numerator"]),
                        int(item["minimum_positive_coefficient"]["lower"]
                            ["denominator"]),
                    )
                    for item in tetrahedra
                ]
                observed.append(f"{float(min(minima)):.2e}")
            self.assertEqual(
                tuple(observed), expected_region_minima[level["alpha"]])

    def test_n5_conventional_level_certificate_artifacts(self):
        certificate_dir = PYTHON_DIR.parent / "certificates"
        structure_path = (certificate_dir
                          / "n5-gaffke-bernstein-structure.json")
        certificate_path = certificate_dir / "n5-gaffke-certificate.json"
        structure_bytes = structure_path.read_bytes()
        structure = json.loads(structure_bytes)
        certificate = json.loads(certificate_path.read_text())

        self.assertEqual(
            structure["face_order_verification"],
            "exact_polynomial_ideal_membership",
        )
        self.assertEqual(
            certificate["structure_sha256"],
            hashlib.sha256(structure_bytes).hexdigest(),
        )
        self.assertEqual(
            [level["alpha"] for level in certificate["levels"]],
            ["0.01", "0.05", "0.10"],
        )
        self.assertEqual(
            {
                region: [simplex["structural_zero_count"]
                         for simplex in record["simplices"]]
                for region, record in structure["regions"].items()
            },
            {
                "A": [1],
                "B": [66, 15, 39],
                "C": [59, 35, 92, 59, 92],
                "D": [15, 39, 66],
                "E": [1],
            },
        )
        expected_region_minima = {
            "0.01": ("8.16e-13", "7.45e-20", "3.23e-21",
                     "7.81e-16", "3.47e-04"),
            "0.05": ("6.70e-09", "1.89e-14", "6.24e-16",
                     "3.57e-12", "1.38e-03"),
            "0.10": ("3.81e-07", "4.93e-12", "1.36e-13",
                     "1.34e-10", "2.20e-03"),
        }
        for level in certificate["levels"]:
            self.assertEqual(level["factor_bits"], 240)
            self.assertEqual(level["interval_bits"], 256)
            observed = []
            for region_name in ("A", "B", "C", "D", "E"):
                simplices = level["polynomial_regions"][region_name][
                    "simplices"]
                for item in simplices:
                    determinant = item["affine_determinant"]
                    lower = Fraction(
                        int(determinant["lower"]["numerator"]),
                        int(determinant["lower"]["denominator"]),
                    )
                    upper = Fraction(
                        int(determinant["upper"]["numerator"]),
                        int(determinant["upper"]["denominator"]),
                    )
                    self.assertFalse(lower <= 0 <= upper)
                minima = [
                    Fraction(
                        int(item["minimum_positive_coefficient"]["lower"]
                            ["numerator"]),
                        int(item["minimum_positive_coefficient"]["lower"]
                            ["denominator"]),
                    )
                    for item in simplices
                ]
                observed.append(f"{float(min(minima)):.2e}")
            self.assertEqual(
                tuple(observed), expected_region_minima[level["alpha"]])

    def test_generated_rows_match_the_manuscript_table(self):
        rows = summarize(default_paths())
        observed = [
            (row["nominal_percent"], row["n"], row["table_display"],
             row["certified_examples"])
            for row in rows
        ]
        expected = [
            ("30", 50, "0.29815", 3),
            ("32", 100, "0.30956", 9),
            ("35", 200, "0.33729", 1),
            ("35", 400, "0.32630", 1),
            ("37", 200, "0.36079", 6),
            ("37", 400, "0.36777", 13),
        ]
        self.assertEqual(observed, expected)

        manuscript = (PYTHON_DIR.parents[1] / "paper" / "stringer.tex"
                      ).read_text()
        for percent, n, display, count in expected:
            row = (f"${percent}\\%$ & ${n}$ & ${display}$ & "
                   f"${count}$\\\\")
            self.assertIn(row, manuscript)


if __name__ == "__main__":
    unittest.main()
