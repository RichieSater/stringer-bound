"""Targeted tests for the rigorous three-exponential interior certificate."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent / "certificates"
    / "three-exponential-interior-certificate.json")
sys.path.insert(0, str(PYTHON_DIR))

from three_exponential_interior_certificate import certify_cell  # noqa: E402


class ThreeExponentialInteriorTests(unittest.TestCase):
    def test_certificate_summary_has_complete_partition(self):
        data = json.loads(CERTIFICATE.read_text())
        partition = data["partition"]
        self.assertEqual(partition["coarse_grid_shape"], [150, 150])
        self.assertEqual(partition["coarse_leaf_count"], 15587)
        self.assertEqual(partition["fine_leaf_count"], 27652)
        self.assertEqual(partition["total_leaf_count"], 43239)
        self.assertEqual(
            data["leaf_trace_sha256"],
            "ee5f5131586001d5ede0f8595da9eeba0d0ba4fb1a3673f7806b7a06e65749e4",
        )

    def test_representative_and_weak_margin_cells_certify(self):
        boxes = (
            (Fraction(1), Fraction(101, 100),
             Fraction(1), Fraction(101, 100)),
            (Fraction(32, 25), Fraction(129, 100),
             Fraction(199, 50), Fraction(399, 100)),
            (Fraction(69, 50), Fraction(139, 100),
             Fraction(153, 50), Fraction(307, 100)),
            (Fraction(399, 100), Fraction(4),
             Fraction(1), Fraction(101, 100)),
        )
        for box in boxes:
            certified, margins, _ = certify_cell(box)
            self.assertTrue(certified, box)
            self.assertTrue(all(margin > 0 for margin in margins), box)

    def test_a_failed_coarse_cell_is_closed_by_four_refinements(self):
        coarse = (
            Fraction(1), Fraction(51, 50),
            Fraction(1), Fraction(51, 50),
        )
        certified, _, _ = certify_cell(coarse)
        self.assertFalse(certified)
        step = Fraction(1, 100)
        for i in range(2):
            for j in range(2):
                subbox = (
                    Fraction(1) + i * step,
                    Fraction(1) + (i + 1) * step,
                    Fraction(1) + j * step,
                    Fraction(1) + (j + 1) * step,
                )
                subcertified, _, _ = certify_cell(subbox)
                self.assertTrue(subcertified, subbox)


if __name__ == "__main__":
    unittest.main()
