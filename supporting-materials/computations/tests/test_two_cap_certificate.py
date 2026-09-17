from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = PYTHON_DIR.parent / "certificates" / "two-cap-certificate.json"
sys.path.insert(0, str(PYTHON_DIR))

from gaffke import dirichlet_average_tail_exact  # noqa: E402
from two_cap_certificate import (  # noqa: E402
    build_certificate,
    derivative_coefficients,
    two_cap_tail_exact,
    two_cap_uniform_bound,
)


class TwoCapCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.committed = json.loads(CERTIFICATE.read_text())

    def test_certificate_regenerates(self) -> None:
        # Use the committed full result as the authoritative slow regression;
        # a small regeneration here keeps the unit suite fast.  The Make target
        # performs the byte-for-byte n<=200 regeneration.
        small = build_certificate(n_max=8, factor_bits=64)
        self.assertEqual([r["alpha"] for r in small["results"]],
                         ["0.10", "0.05", "0.01"])
        self.assertEqual([r["q_min"]["decimal"] for r in self.committed["results"]],
                         ["0.416666666666666667", "0.333333333333333333",
                          "0.250000000000000000"])
        self.assertTrue(all(r["endpoint_inequalities_checked"] == 19900
                            for r in self.committed["results"]))

    def test_two_cap_formula_matches_exact_divided_difference(self) -> None:
        knots = (Fraction(0), Fraction(1, 4), Fraction(1, 2),
                 Fraction(3, 4), Fraction(1))
        s = Fraction(5, 8)
        self.assertEqual(two_cap_tail_exact(knots, s),
                         dirichlet_average_tail_exact(s, knots))

    def test_derivative_coefficients_have_only_negative_to_positive_change(self) -> None:
        for n in range(2, 15):
            for r in range(1, n):
                coeffs = derivative_coefficients(
                    n, r, Fraction(7, 13), Fraction(2, 17))
                signs = [1 if value > 0 else -1 for value in coeffs if value]
                changes = [(a, b) for a, b in zip(signs, signs[1:]) if a != b]
                self.assertLessEqual(len(changes), 1)
                self.assertNotIn((1, -1), changes)

    def test_uniform_bound_dominates_constructed_exact_caps(self) -> None:
        for n in range(2, 8):
            total = sum(range(1, n + 2))
            weights = tuple(Fraction(i, total) for i in range(1, n + 2))
            q = Fraction(1, 2)
            K = weights[n] + weights[n - 1] * (1 - q)
            # An equal mixture of all prefix-step budget vertices.
            lambdas = [Fraction(1, n - 1)] * (n - 1)
            prefixes = []
            running = Fraction(0)
            for value in weights[:n - 1]:
                running += value
                prefixes.append(running)
            increments = [lambdas[r] * K / prefixes[r] for r in range(n - 1)]
            z = [1 + sum(increments[j:]) for j in range(n - 1)]
            S = Fraction(1, 2) / z[0]
            knots = tuple([1 - S * value for value in z]
                          + [1 - S * q, Fraction(1)])
            threshold = 1 - S
            # The constructed budget makes threshold the c-barycenter.
            self.assertEqual(sum(c * x for c, x in zip(weights, knots)), threshold)
            exact = two_cap_tail_exact(knots, threshold)
            self.assertLessEqual(exact, two_cap_uniform_bound(weights, q))

    def test_reported_worst_cases_and_positive_margins(self) -> None:
        expected = {"0.10": (7, 6), "0.05": (9, 8), "0.01": (12, 11)}
        for result in self.committed["results"]:
            worst = result["worst_endpoint_case"]
            self.assertEqual((worst["n"], worst["r"]), expected[result["alpha"]])
            self.assertGreater(int(worst["strict_margin_below_alpha"]["numerator"]), 0)


if __name__ == "__main__":
    unittest.main()
