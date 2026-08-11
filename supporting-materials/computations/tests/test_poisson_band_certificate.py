"""Tests for the exact Poisson simultaneous-band certificate."""

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path

from mpmath import mp


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

import stringer  # noqa: E402
from poisson_band_certificate import (  # noqa: E402
    LEVELS,
    bolshev_probability,
    build_certificate,
)


class ExactPoissonFactorTests(unittest.TestCase):
    def test_rational_exponential_bounds_enclose_high_precision_value(self):
        for value in (Fraction(0), Fraction(1, 3), Fraction(7, 4),
                      Fraction(31, 2)):
            lower, upper = stringer.exact_exp_neg_bounds(value)
            # The rational enclosure is much tighter than 100 decimal
            # places, so the independent numerical value needs ample guard
            # digits to avoid mistaking its own rounding for a failed bound.
            with mp.workdps(220):
                exact = mp.exp(-mp.mpf(value.numerator) / value.denominator)
                self.assertLessEqual(mp.mpf(lower.numerator)
                                     / lower.denominator, exact)
                self.assertGreaterEqual(mp.mpf(upper.numerator)
                                        / upper.denominator, exact)

    def test_poisson_brackets_have_exact_endpoint_signs(self):
        alpha = Fraction(1, 20)
        bits = 48
        brackets = stringer.exact_poisson_lambda_brackets(
            "0.05", 7, bits)
        self.assertEqual(len(brackets), 8)
        for j, (lower, upper) in enumerate(brackets):
            self.assertEqual(upper - lower, Fraction(1, 1 << bits))
            lower_cdf, _ = stringer.exact_poisson_cdf_bounds(lower, j)
            _, upper_cdf = stringer.exact_poisson_cdf_bounds(upper, j)
            self.assertGreater(lower_cdf, alpha)
            self.assertLess(upper_cdf, alpha)
        self.assertTrue(all(
            brackets[j][0] < brackets[j + 1][0]
            for j in range(len(brackets) - 1)
        ))

    def test_poisson_input_validation(self):
        with self.assertRaises(ValueError):
            stringer.exact_exp_neg_bounds(Fraction(-1))
        with self.assertRaises(ValueError):
            stringer.exact_poisson_cdf_bounds(Fraction(1), -1)
        with self.assertRaises(ValueError):
            stringer.exact_poisson_lambda_brackets("1", 2)


class BolshevCertificateTests(unittest.TestCase):
    def test_bolshev_special_cases(self):
        a = Fraction(2, 5)
        self.assertEqual(bolshev_probability([a]), a)
        self.assertEqual(bolshev_probability([a, a]), a**2)
        self.assertEqual(
            bolshev_probability([a, Fraction(1)]),
            1 - (1 - a)**2,
        )

    def test_boundary_validation(self):
        with self.assertRaises(ValueError):
            bolshev_probability([Fraction(2, 3), Fraction(1, 3)])
        with self.assertRaises(ValueError):
            bolshev_probability([Fraction(4, 3)])

    def test_conventional_level_frontiers_are_exactly_certified(self):
        certificate = build_certificate()
        actual = {
            level["alpha"]: level["certified_n_range"][1]
            for level in certificate["levels"]
        }
        self.assertEqual(actual, LEVELS)
        for level in certificate["levels"]:
            self.assertGreaterEqual(
                Fraction(int(level["frontier_margin_over_nominal"]["numerator"]),
                         int(level["frontier_margin_over_nominal"]["denominator"])),
                0,
            )
            limitation = level["next_size_limitation"]
            self.assertLess(
                Fraction(int(limitation["event_probability_upper_bound"]["numerator"]),
                         int(limitation["event_probability_upper_bound"]["denominator"])),
                1 - Fraction(level["alpha"]),
            )


if __name__ == "__main__":
    unittest.main()
