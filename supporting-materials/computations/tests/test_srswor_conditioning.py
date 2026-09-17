from __future__ import annotations

import itertools
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = ROOT / "computations" / "certificates" / "srswor-conditioning-certificate.json"
sys.path.insert(0, str(PYTHON_DIR))

from srswor_conditioning import (  # noqa: E402
    adjusted_tail,
    build_certificate,
    distinct_probability,
    minimum_population_size,
)


class SrsworConditioningTests(unittest.TestCase):
    def test_exact_collision_probabilities(self) -> None:
        self.assertEqual(distinct_probability(3, 2), Fraction(2, 3))
        self.assertEqual(distinct_probability(3, 3), Fraction(2, 9))
        self.assertEqual(distinct_probability(10, 1), Fraction(1, 1))
        self.assertEqual(adjusted_tail(Fraction(1, 20), 3, 3), Fraction(1, 90))

    def test_reported_minimum_population_sizes_are_exact(self) -> None:
        expected = {
            4: (10, 5, 6),
            5: (17, 7, 8),
            6: (24, 9, 12),
            7: (33, 12, 16),
        }
        for n, values in expected.items():
            for required, wanted in zip((Fraction(1, 2), Fraction(1, 10), Fraction(1, 5)), values):
                N = minimum_population_size(n, required)
                self.assertEqual(N, wanted)
                self.assertGreaterEqual(distinct_probability(N, n), required)
                self.assertLess(distinct_probability(N - 1, n), required)

    def test_conditioning_identity_and_bound_by_exhaustion(self) -> None:
        # Three different deterministic event shapes.  The theorem applies to
        # any statistic/event; these finite checks guard the implementation of
        # the conditioning factor, not the written probability proof.
        N, n = 4, 3
        tuples = list(itertools.product(range(N), repeat=n))
        distinct = [x for x in tuples if len(set(x)) == n]
        c = distinct_probability(N, n)
        events = (
            lambda x: sum(x) >= 5,
            lambda x: max(x) - min(x) >= 2,
            lambda x: x[0] <= x[1] <= x[2],
        )
        for event in events:
            p_wr = Fraction(sum(event(x) for x in tuples), len(tuples))
            p_joint = Fraction(sum(event(x) for x in distinct), len(tuples))
            p_wor = Fraction(sum(event(x) for x in distinct), len(distinct))
            self.assertEqual(p_wor, p_joint / c)
            self.assertLessEqual(p_wor, p_wr / c)

    def test_certificate_matches_generator(self) -> None:
        committed = json.loads(CERTIFICATE.read_text())
        self.assertEqual(committed, build_certificate())


if __name__ == "__main__":
    unittest.main()
