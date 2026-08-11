"""Tests for the all-sample-size Gaffke safeguard."""

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.stats import beta


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from gaffke import (  # noqa: E402
    dirichlet_average_cdf,
    dirichlet_average_quantile,
    dirichlet_average_tail_exact,
    gaffke_quantile_certificate,
    gaffke_upper_bound,
    safeguarded_stringer_bound,
)
from stringer import factor_prefix, factor_values, stringer_bound  # noqa: E402


class DirichletAverageTests(unittest.TestCase):
    def test_uniform_two_knot_case(self):
        knots = [0.2, 0.8]
        self.assertAlmostEqual(dirichlet_average_cdf(0.5, knots), 0.5)
        self.assertAlmostEqual(
            dirichlet_average_quantile(knots, 0.95), 0.77, places=14)

    def test_repeated_knots_reduce_to_beta_distribution(self):
        # The sum of three out of five exchangeable Dirichlet coordinates is
        # Beta(3, 2).  Repeated knots exercise the confluent B-spline path.
        knots = [0.0, 0.0, 1.0, 1.0, 1.0]
        for value in (0.1, 0.4, 0.8):
            self.assertAlmostEqual(
                dirichlet_average_cdf(value, knots),
                beta.cdf(value, 3, 2), places=13)
        self.assertAlmostEqual(
            dirichlet_average_quantile(knots, 0.95),
            beta.ppf(0.95, 3, 2), places=13)

    def test_all_equal_knots_are_degenerate(self):
        self.assertEqual(dirichlet_average_cdf(0.3, [0.3] * 5), 1.0)
        self.assertEqual(
            dirichlet_average_quantile([0.3] * 5, 0.95), 0.3)

    def test_exact_confluent_tail_matches_beta_cases(self):
        x = Fraction(2, 5)
        # With two zeros and three ones the average is Beta(3, 2).
        actual = dirichlet_average_tail_exact(x, [0, 0, 1, 1, 1])
        # Its CDF is 4*x^3-3*x^4.
        expected = Fraction(1) - (4 * x**3 - 3 * x**4)
        self.assertEqual(actual, expected)

        # With n zeros and a terminal one the average is Beta(1, n).
        self.assertEqual(
            dirichlet_average_tail_exact(Fraction(1, 4), [0] * 4 + [1]),
            Fraction(3, 4) ** 4,
        )

    def test_exact_quantile_certificate_checks_both_tail_signs(self):
        certificate = gaffke_quantile_certificate(
            ["1", "0.4", "0.1"], 20, "0.05")
        self.assertGreaterEqual(
            certificate.tail_at_lower, certificate.alpha)
        self.assertLessEqual(
            certificate.tail_at_upper, certificate.alpha)
        self.assertLessEqual(certificate.width, Fraction(6, 1 << 48))


class GaffkeSafeguardTests(unittest.TestCase):
    def test_stringer_factor_prefix_matches_full_factor_table(self):
        for method in ("binomial", "poisson"):
            prefix = factor_prefix(12, "0.05", 3, method, dps=60)
            full = factor_values(12, "0.05", method, dps=60)
            self.assertEqual(
                [float(value) for value, _ in prefix],
                [float(value) for value in full[:4]],
            )
        omitted = stringer_bound(
            [1.0, 0.4, 0.1], 100, "0.05", "poisson", dps=60)
        explicit = stringer_bound(
            [1.0, 0.4, 0.1] + [0.0] * 97,
            100, "0.05", "poisson", dps=60)
        self.assertEqual(float(omitted), float(explicit))

    def test_bernoulli_samples_equal_clopper_pearson(self):
        alpha = "0.05"
        for n in (2, 5, 12):
            factors = factor_values(n, alpha, "binomial", dps=80)
            for errors in range(n + 1):
                observed = [1.0] * errors  # zeros are intentionally omitted
                actual = gaffke_upper_bound(observed, n, alpha)
                self.assertAlmostEqual(
                    actual, float(factors[errors]), places=12,
                    msg=(n, errors))

    def test_safeguard_is_the_maximum_of_its_components(self):
        result = safeguarded_stringer_bound(
            [1.0, 0.4, 0.1], 100, "0.05", "poisson")
        self.assertEqual(result.safeguarded,
                         max(result.stringer, result.gaffke))
        self.assertGreaterEqual(result.uplift, 0.0)
        self.assertIn(result.governing_bound, ("stringer", "gaffke"))

    def test_certified_small_n_cases_need_no_uplift(self):
        samples = (
            [],
            [0.1],
            [1.0],
            [0.9, 0.4],
            [1.0, 0.8, 0.2],
        )
        for n in (3, 4, 5):
            for alpha in ("0.01", "0.05", "0.10"):
                for observed in samples:
                    if len(observed) > n:
                        continue
                    result = safeguarded_stringer_bound(
                        observed, n, alpha, "binomial")
                    # The exact theorem says Stringer >= Gaffke.  The
                    # implementation deliberately reports a dyadic upper
                    # enclosure of Gaffke, so equality cases can acquire a
                    # sub-1e-13 numerical uplift without changing the
                    # mathematical pointwise comparison.
                    self.assertLessEqual(
                        result.gaffke - result.stringer, 6e-14,
                        msg=(n, alpha, observed, result))
                    self.assertLessEqual(result.uplift, 6e-14)
                    self.assertAlmostEqual(
                        result.stringer,
                        float(stringer_bound(
                            observed, n, alpha, "binomial", dps=80)),
                        places=14)

    def test_all_n_one_cap_region_needs_no_uplift(self):
        # These samples have their largest taint below binomial Stringer, so
        # the analytic one-cap theorem applies without a sample-size cutoff.
        for n in (6, 20, 100):
            for alpha in ("0.01", "0.05", "0.10"):
                observed = ["0.002", "0.001"]
                result = safeguarded_stringer_bound(
                    observed, n, alpha, "binomial")
                self.assertGreaterEqual(
                    result.stringer, max(map(float, observed)))
                # Gaffke is reported as a dyadic upper enclosure, so allow
                # only the certificate-scale rounding seen at equality.
                self.assertLessEqual(
                    result.gaffke - result.stringer, 6e-14,
                    msg=(n, alpha, result))

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            gaffke_upper_bound([1.1], 4, "0.05")
        with self.assertRaises(ValueError):
            gaffke_upper_bound([0.2] * 5, 4, "0.05")
        with self.assertRaises(ValueError):
            safeguarded_stringer_bound([0.2], 4, "0.05", "unknown")


if __name__ == "__main__":
    unittest.main()
