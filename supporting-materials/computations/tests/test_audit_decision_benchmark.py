"""Regression tests for the exact decision-oriented benchmark."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE_PATH = (
    PYTHON_DIR.parent / "certificates" / "audit-decision-benchmark.json")
sys.path.insert(0, str(PYTHON_DIR))

from audit_decision_benchmark import build_benchmark  # noqa: E402


def fraction(record):
    return Fraction(int(record["numerator"]), int(record["denominator"]))


class AuditDecisionBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.generated = build_benchmark()
        cls.committed = json.loads(CERTIFICATE_PATH.read_text())

    def test_regenerates_committed_artifact(self):
        self.assertEqual(self.generated, self.committed)

    def test_fixed_case_suite_and_exact_relationships(self):
        self.assertEqual(self.generated["case_count"], 10)
        self.assertEqual(len(self.generated["cases"]), 10)
        for case in self.generated["cases"]:
            methods = case["methods"]
            poisson = fraction(methods[
                "ordinary_poisson_stringer_reference"]["rate_upper"])
            gaffke = fraction(methods[
                "gaffke_valid_endpoint"]["rate_upper"])
            safeguard = fraction(methods[
                "poisson_stringer_gaffke_max"]["rate_upper"])
            self.assertGreater(poisson, gaffke)
            self.assertEqual(safeguard, poisson)
            self.assertEqual(
                methods["poisson_stringer_gaffke_max"]
                ["governing_component_certified"],
                "ordinary_poisson_stringer",
            )
            for comparison in case["threshold_switch_bands"].values():
                width = fraction(comparison["width_rate"])
                self.assertEqual(
                    comparison["has_nonempty_switch_band"], width > 0)

    def test_zero_anchor_preserves_all_zero_result(self):
        zero_cases = [case for case in self.generated["cases"]
                      if case["nonzero_taint_count"] == 0]
        self.assertEqual(len(zero_cases), 2)
        for case in zero_cases:
            methods = case["methods"]
            ordinary = fraction(methods[
                "ordinary_poisson_stringer_reference"]["rate_upper"])
            anchored = fraction(methods[
                "zero_anchor_calibrated_factorwise_cap"]["rate_upper"])
            self.assertEqual(ordinary, anchored)


if __name__ == "__main__":
    unittest.main()
