"""Tests for the exact five-coordinate repeated-knot face certificate."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent
    / "certificates"
    / "five-coordinate-block-faces-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from five_coordinate_block_faces import (  # noqa: E402
    build_certificate,
    derive_block_face_polynomials,
)


class FiveCoordinateBlockFaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected = json.loads(CERTIFICATE.read_text())
        cls.regenerated = build_certificate()
        cls.derived = derive_block_face_polynomials()["cases"]

    def test_checked_in_certificate_is_exactly_regenerated(self):
        self.assertEqual(self.regenerated, self.expected)

    def test_every_certificate_coefficient_has_exact_sign(self):
        pair32 = self.expected["three_below_top_pair"]
        self.assertEqual(
            pair32["compactness_sign_audit"]["positive_coefficients"], 89
        )
        self.assertEqual(
            pair32["compactness_sign_audit"]["negative_coefficients"], 0
        )
        multiplier32 = pair32["bernstein_basis"]["multiplier_S"]
        residual32 = pair32["bernstein_basis"]["residual_H"]
        self.assertEqual(multiplier32["positive_coefficients"], 111)
        self.assertEqual(multiplier32["negative_coefficients"], 0)
        self.assertEqual(residual32["positive_coefficients"], 727)
        self.assertEqual(residual32["zero_coefficients"], 729)
        self.assertEqual(residual32["negative_coefficients"], 0)
        self.assertEqual(
            multiplier32["sha256_all_coefficients"],
            "4dc342db208460a61249a85d8248d03c741c58371231bf848829e96fa0572900",
        )
        self.assertEqual(
            residual32["sha256_all_coefficients"],
            "4a7eb401dea31627cab53b245a1217905b24d6b1d77214844641c2a5d2e88403",
        )

        pair23 = self.expected["two_below_top_pair"]["bernstein_basis"]
        quotient23 = pair23["quotient_S"]
        remainder23 = pair23["remainder_H"]
        self.assertEqual(quotient23["positive_coefficients"], 8733)
        self.assertEqual(quotient23["zero_coefficients"], 651)
        self.assertEqual(quotient23["negative_coefficients"], 0)
        self.assertEqual(remainder23["positive_coefficients"], 4424)
        self.assertEqual(remainder23["zero_coefficients"], 1176)
        self.assertEqual(remainder23["negative_coefficients"], 0)
        self.assertEqual(
            quotient23["sha256_all_coefficients"],
            "22000cfb50238ca919f3a3026ee9935d544f02101074e7ed8ef8d41a4ff3a942",
        )
        self.assertEqual(
            remainder23["sha256_all_coefficients"],
            "d9ee92094145995a706f1263d96dd105415a3dc8a1afa84af1202cb3a6118571",
        )

        triple23 = self.expected["two_below_top_triple"]
        self.assertEqual(
            triple23["compactness_sign_audit"]["positive_coefficients"], 19
        )
        self.assertEqual(
            triple23["compactness_sign_audit"]["negative_coefficients"], 0
        )
        multiplier_triple = triple23["bernstein_basis"]["multiplier_S"]
        residual_triple = triple23["bernstein_basis"]["residual_H"]
        self.assertEqual(multiplier_triple["positive_coefficients"], 38)
        self.assertEqual(multiplier_triple["negative_coefficients"], 0)
        self.assertEqual(residual_triple["positive_coefficients"], 73)
        self.assertEqual(residual_triple["zero_coefficients"], 107)
        self.assertEqual(residual_triple["negative_coefficients"], 0)
        self.assertEqual(
            multiplier_triple["sha256_all_coefficients"],
            "e1fe9343b0233f274d8817f5534e005783cf37bbc922dc943f8249b7b4ddbf37",
        )
        self.assertEqual(
            residual_triple["sha256_all_coefficients"],
            "f17649607e13bb21dbe75b0197bdc952607fd36520acad828c06cb935fbb8a67",
        )

    def _check_fixture(self, name, values, expected_cap, expected_gamma0):
        case = self.derived[name]
        point = dict(zip(case["variables"], values, strict=True))
        cap = Fraction(case["cap"].subs(point))
        gamma0 = Fraction(case["centroid"][0].subs(point))
        if name == "two_below_top_triple":
            gap = Fraction(
                case["centroid"][2].subs(point)
                - case["centroid"][1].subs(point)
            )
        else:
            gap = Fraction(
                case["centroid"][3].subs(point)
                - case["centroid"][2].subs(point)
            )
        self.assertEqual(cap, expected_cap)
        self.assertEqual(gamma0, expected_gamma0)
        self.assertGreater(gap, 0)
        self.assertLess(cap, 1 - (1 - gamma0) ** 4)

    def test_exact_rational_three_below_top_pair_fixture(self):
        self._check_fixture(
            "three_below_top_pair",
            (Fraction(1, 2), Fraction(1, 2), Fraction(1, 10)),
            Fraction(1597, 106722),
            Fraction(311, 9988),
        )

    def test_exact_rational_two_below_top_pair_fixture(self):
        self._check_fixture(
            "two_below_top_pair",
            (Fraction(1, 2), Fraction(1, 2), Fraction(1, 10)),
            Fraction(871, 15246),
            Fraction(233, 5852),
        )

    def test_exact_rational_two_below_top_triple_fixture(self):
        self._check_fixture(
            "two_below_top_triple",
            (Fraction(1, 2), Fraction(1, 10)),
            Fraction(21871, 287496),
            Fraction(817, 17468),
        )


if __name__ == "__main__":
    unittest.main()
