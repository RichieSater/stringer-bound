"""Regression tests for the exact binomial-certificate path."""

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from math import comb
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

import stringer  # noqa: E402
from coverage import coverage_exact  # noqa: E402
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
