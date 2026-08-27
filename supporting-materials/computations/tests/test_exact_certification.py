"""Regression tests for the exact binomial-certificate path."""

from __future__ import annotations

import collections
import gzip
import hashlib
import json
import sys
import unittest
from fractions import Fraction
from math import comb
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
REPOSITORY_ROOT = PYTHON_DIR.parents[2]
sys.path.insert(0, str(PYTHON_DIR))

import stringer  # noqa: E402
import derive_n6_bernstein_structure as n6_structure  # noqa: E402
import derive_n7_bernstein_structure as n7_structure  # noqa: E402
import n6_gaffke_certificate as n6_certificate  # noqa: E402
import n7_gaffke_certificate as n7_certificate  # noqa: E402
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


class SubmissionPolicyTests(unittest.TestCase):
    def test_ai_disclosure_occurs_once_and_only_in_principal_tex(self):
        principal = (
            REPOSITORY_ROOT / "supporting-materials" / "paper" / "stringer.tex"
        )
        circulation_documents = [
            REPOSITORY_ROOT / "README.md",
            REPOSITORY_ROOT / "FINDINGS.md",
            *(
                REPOSITORY_ROOT / "supporting-materials"
            ).rglob("*.md"),
            *(
                REPOSITORY_ROOT / "supporting-materials" / "paper"
            ).rglob("*.tex"),
        ]
        principal_text = principal.read_text(encoding="utf-8")
        self.assertEqual(principal_text.count(r"\section*{AI disclosure}"), 1)
        self.assertEqual(principal_text.count("Anthropic Claude"), 1)
        self.assertEqual(principal_text.count("OpenAI Codex"), 1)

        prohibited_status_phrases = (
            "AI-assisted",
            "AI-assisted review",
            "".join(("needs", " human", " validation")),
            "claimed complete solution",
        )
        for path in circulation_documents:
            text = path.read_text(encoding="utf-8")
            for phrase in prohibited_status_phrases:
                self.assertNotIn(phrase, text, path)
            if path != principal:
                self.assertNotIn("AI disclosure", text, path)
                self.assertNotIn("Anthropic Claude", text, path)
                self.assertNotIn("OpenAI Codex", text, path)


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
    def test_n7_empty_generator_product_is_one(self):
        self.assertEqual(n7_structure._product(value for value in ()), "1")

    def test_n7_sparse_power_to_bernstein_transform(self):
        n7_certificate.ctx.prec = n7_certificate.BALL_BITS
        polynomial = {
            (0, 0, 0, 0, 0, 0): n7_certificate.arb(2),
            (1, 0, 0, 0, 0, 0): n7_certificate.arb(3),
            (0, 2, 0, 0, 0, 0): n7_certificate.arb(5),
            (1, 0, 0, 0, 0, 1): n7_certificate.arb(-7),
        }
        coefficients = n7_certificate._bernstein_coefficients(
            polynomial, 2)
        for powers, value in coefficients.items():
            expected = (
                Fraction(2)
                + Fraction(3 * powers[0], 2)
                + Fraction(5 * comb(powers[1], 2))
                - Fraction(7 * powers[0] * powers[5], 2)
            )
            lower = value.lower().fmpq()
            upper = value.upper().fmpq()
            self.assertLessEqual(
                Fraction(int(lower.p), int(lower.q)), expected)
            self.assertGreaterEqual(
                Fraction(int(upper.p), int(upper.q)), expected)

    def test_n6_generic_derivative_chunks_are_exhaustive(self):
        for order, chunk_count in ((1, 6), (2, 6), (3, 6),
                                   (4, 6), (5, 6), (6, 3)):
            with self.subTest(order=order, chunk_count=chunk_count):
                full_script, full_count = (
                    n6_structure._singular_differential_checks(order))
                chunks = [n6_structure._singular_differential_checks(
                    order, (index, chunk_count))
                    for index in range(chunk_count)]
                self.assertEqual(
                    full_count, comb(5 + order - 1, 5))
                self.assertEqual(
                    sum(count for _, count in chunks), full_count)
                full_checks = {
                    line for line in full_script.splitlines()
                    if line.startswith("poly D")
                }
                chunk_checks = [
                    {line for line in script.splitlines()
                     if line.startswith("poly D")}
                    for script, _ in chunks
                ]
                self.assertEqual(set().union(*chunk_checks), full_checks)
                for left in range(chunk_count):
                    for right in range(left + 1, chunk_count):
                        self.assertFalse(
                            chunk_checks[left] & chunk_checks[right])

    def test_n6_sparse_power_to_bernstein_transform(self):
        polynomial = n6_certificate.SimplexPolynomial({
            (0, 0, 0, 0, 0): 2,
            (1, 0, 0, 0, 0): 3,
            (0, 2, 0, 0, 0): 5,
            (1, 0, 0, 0, 1): -7,
        })
        coefficients = n6_certificate._bernstein_coefficients(
            polynomial, 2)
        for index, value in coefficients.items():
            alpha = index[:5]
            expected = (
                Fraction(2)
                + Fraction(3 * alpha[0], 2)
                + Fraction(5 * comb(alpha[1], 2))
                - Fraction(7 * alpha[0] * alpha[4], 2)
            )
            self.assertEqual(value.lower, expected)
            self.assertEqual(value.upper, expected)

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

    def test_n6_conventional_level_certificate_artifacts(self):
        certificate_dir = PYTHON_DIR.parent / "certificates"
        structure_path = (certificate_dir
                          / "n6-gaffke-bernstein-structure.json")
        certificate_path = certificate_dir / "n6-gaffke-certificate.json"
        structure_bytes = structure_path.read_bytes()
        structure = json.loads(structure_bytes)
        certificate = json.loads(certificate_path.read_text())

        self.assertEqual(
            structure["face_order_verification"],
            "exact_generic_differential_ideal_membership",
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
                "B": [456, 456, 302, 162, 57],
                "C": [666, 666, 911, 911, 911, 911,
                      813, 358, 554, 813],
                "D": [813, 358, 554, 813, 666, 666,
                      911, 911, 911, 911],
                "E": [456, 456, 57, 302, 162],
                "F": [1],
            },
        )
        for record in structure["regions"].values():
            for simplex in record["simplices"]:
                for proof in simplex["face_order_proofs"]:
                    self.assertIn(
                        "over_QQ(b,c,d,e,f,g)", proof["verification"])
                    self.assertGreaterEqual(
                        proof["derivative_reductions_checked"], 1)

        expected_region_minima = {
            "0.01": ("2.02e-15", "9.63e-26", "1.22e-30",
                     "6.05e-29", "6.29e-20", "2.41e-04"),
            "0.05": ("1.23e-10", "1.17e-18", "6.57e-23",
                     "4.01e-22", "1.66e-15", "9.57e-04"),
            "0.10": ("1.80e-08", "1.85e-15", "1.93e-19",
                     "4.24e-19", "1.35e-13", "1.53e-03"),
        }
        for level in certificate["levels"]:
            self.assertEqual(level["factor_bits"], 320)
            self.assertEqual(level["interval_bits"], 384)
            rank_certificate = level["face_normal_rank_certificate"]
            self.assertEqual(
                rank_certificate["generator_sets_certified"], 16)
            self.assertEqual(
                rank_certificate[
                    "generic_region_face_order_checks_linked"], 26)
            self.assertGreater(
                Fraction(
                    int(rank_certificate["minimum_absolute_normal_minor"]
                        ["numerator"]),
                    int(rank_certificate["minimum_absolute_normal_minor"]
                        ["denominator"]),
                ),
                0,
            )
            for record in rank_certificate["sets"]:
                determinant = record["normal_minor_determinant"]
                lower = Fraction(int(determinant["lower"]["numerator"]),
                                 int(determinant["lower"]["denominator"]))
                upper = Fraction(int(determinant["upper"]["numerator"]),
                                 int(determinant["upper"]["denominator"]))
                self.assertFalse(lower <= 0 <= upper)
            observed = []
            for region_name in ("A", "B", "C", "D", "E", "F"):
                region = level["polynomial_regions"][region_name]
                vertex_certificate = region["vertex_region_certificate"]
                self.assertGreater(vertex_certificate["vertex_count"], 0)
                self.assertGreater(
                    vertex_certificate[
                        "ordered_domain_and_region_inequalities_checked"],
                    0,
                )
                slack = vertex_certificate["minimum_strict_vertex_slack"]
                self.assertIsNotNone(slack)
                self.assertGreater(
                    Fraction(int(slack["lower"]["numerator"]),
                             int(slack["lower"]["denominator"])),
                    0,
                )
                simplices = region["simplices"]
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
                tuple(observed),
                expected_region_minima[level["alpha"]],
            )

            chain = level["triangulation_chain_certificate"]
            self.assertEqual(chain["simplex_count"], 32)
            self.assertEqual(chain["distinct_vertex_pairs_certified"], 210)
            self.assertEqual(
                chain[
                    "internal_facets_paired_with_opposite_orientation"], 48)
            self.assertEqual(chain["unpaired_outer_boundary_facets"], 96)
            self.assertEqual(
                chain["outer_boundary_facets_by_hyperplane"],
                {label: 16 for label in (
                    "x=0", "x=y", "y=z", "z=w", "w=u", "u=1")},
            )
            volume = chain["normalized_oriented_volume"]
            volume_lower = Fraction(int(volume["lower"]["numerator"]),
                                    int(volume["lower"]["denominator"]))
            volume_upper = Fraction(int(volume["upper"]["numerator"]),
                                    int(volume["upper"]["denominator"]))
            self.assertLess(volume_lower, 1)
            self.assertGreater(volume_upper, 1)
            self.assertGreater(volume_lower, Fraction(1, 2))
            self.assertLess(volume_upper, Fraction(3, 2))
            self.assertEqual(chain["integer_relative_chain_degree"], 1)

    def test_n7_conventional_level_certificate_artifacts(self):
        certificate_dir = PYTHON_DIR.parent / "certificates"
        structure_path = (
            certificate_dir / "n7-gaffke-bernstein-structure.json.gz")
        certificate_path = certificate_dir / "n7-gaffke-certificate.json"
        structure_bytes = structure_path.read_bytes()
        with gzip.open(structure_path, "rt") as handle:
            structure = json.load(handle)
        certificate = json.loads(certificate_path.read_text())

        self.assertEqual(structure["schema_version"], 2)
        self.assertEqual(
            structure["face_order_verification"],
            "exact_generic_I_adic_factor_order_or_Singular_quotient_Horner_"
            "ideal_power_membership",
        )
        generic_face_certificates = structure[
            "generic_face_order_certificates"]
        self.assertEqual(len(generic_face_certificates), 26)
        self.assertEqual(
            sum(record["verification"].startswith("exact_I_adic_")
                for record in generic_face_certificates),
            22,
        )
        self.assertEqual(
            sum("Singular" in record["verification"]
                for record in generic_face_certificates),
            4,
        )
        self.assertEqual(
            certificate["structure_sha256"],
            hashlib.sha256(structure_bytes).hexdigest(),
        )
        self.assertEqual(
            [level["alpha"] for level in certificate["levels"]],
            ["0.01", "0.05", "0.10"],
        )
        expected_region_minima = {
            "0.01": ("4.80e-18", "1.11e-31", "3.65e-40",
                     "3.21e-42", "4.35e-37", "3.41e-24",
                     "1.76e-04"),
            "0.05": ("2.21e-12", "6.73e-23", "5.89e-30",
                     "3.54e-32", "9.75e-29", "5.17e-19",
                     "7.01e-04"),
            "0.10": ("8.29e-10", "6.54e-19", "2.41e-25",
                     "1.13e-27", "4.97e-25", "9.07e-17",
                     "1.12e-03"),
        }
        self.assertEqual(
            {
                region: (
                    record["degree"], record["simplex_count"],
                    record["total_structural_zero_count"])
                for region, record in structure["regions"].items()
            },
            {
                "A": (7, 1, 1),
                "B": (12, 6, 10462),
                "C": (15, 15, 99507),
                "D": (16, 20, 193352),
                "E": (15, 15, 99507),
                "F": (12, 6, 10462),
                "G": (7, 1, 1),
            },
        )
        generic_keys = set()
        stored_pivots = {}
        face_conditions = 0
        verification_counts = collections.Counter()
        for record in structure["regions"].values():
            for simplex in record["simplices"]:
                for proof in simplex["face_order_proofs"]:
                    face_conditions += 1
                    verification_counts[proof["verification"]] += 1
                    self.assertEqual(
                        len(proof["proof_source_inverse_pivot_columns"]),
                        len(proof["proof_source_ideal_generators"]),
                    )
                    self.assertEqual(
                        len(set(proof[
                            "proof_source_inverse_pivot_columns"])),
                        len(proof["proof_source_inverse_pivot_columns"]),
                    )
                    key = (
                        proof["proof_source_region"],
                        tuple(proof["proof_source_ideal_generators"]),
                        proof["vanishing_order"],
                    )
                    generic_keys.add(key)
                    stored_pivots.setdefault(key, set()).add(tuple(
                        proof["proof_source_inverse_pivot_columns"]))
        self.assertEqual(face_conditions, 322)
        self.assertEqual(len(generic_keys), 26)
        self.assertEqual(
            verification_counts,
            collections.Counter({
                "exact_I_adic_factor_order_using_N_equals_F_times_P_"
                "and_domain_associated_graded_ring": 258,
                "exact_generic_Singular_quotient_Horner_"
                "ideal_power_membership_over_QQ(b,c,d,e,f,g,h)": 64,
            }),
        )
        for key in generic_keys:
            _, pivot = n7_structure._inverse_face_map(key[1])
            expected = tuple(
                n7_structure.COORDINATE_NAMES[index] for index in pivot)
            self.assertEqual(stored_pivots[key], {expected})

        for level in certificate["levels"]:
            self.assertEqual(level["factor_bits"], 512)
            self.assertEqual(level["arb_precision_bits"], 768)
            global_minimum = level[
                "minimum_positive_bernstein_lower_bound"]
            self.assertGreater(
                Fraction(int(global_minimum["numerator"]),
                         int(global_minimum["denominator"])),
                0,
            )
            observed_region_minima = []
            for region in "ABCDEFG":
                minimum = level["polynomial_regions"][region][
                    "minimum_positive_lower_bound"]
                value = Fraction(
                    int(minimum["numerator"]),
                    int(minimum["denominator"]),
                )
                observed_region_minima.append(f"{float(value):.2e}")
            self.assertEqual(
                tuple(observed_region_minima),
                expected_region_minima[level["alpha"]],
            )
            rank = level["face_normal_rank_certificate"]
            self.assertEqual(rank["generator_sets_certified"], 22)
            self.assertEqual(
                rank["proof_source_generator_sets_certified"], 42)
            self.assertEqual(
                rank["simplex_face_order_conditions_linked"], 322)
            self.assertEqual(
                rank["distinct_generic_ideal_power_proofs_linked"], 26)
            rank_minimum = rank[
                "minimum_absolute_normal_minor_lower_bound"]
            self.assertGreater(
                Fraction(int(rank_minimum["numerator"]),
                         int(rank_minimum["denominator"])),
                0,
            )
            for record in rank["sets"] + rank["proof_source_sets"]:
                determinant = record["normal_minor_determinant"]
                lower = determinant["lower"]
                upper = determinant["upper"]
                lower = Fraction(int(lower["numerator"]),
                                 int(lower["denominator"]))
                upper = Fraction(int(upper["numerator"]),
                                 int(upper["denominator"]))
                self.assertFalse(lower <= 0 <= upper)
            for record in rank["proof_source_sets"]:
                self.assertEqual(
                    record["minor_coordinate_columns"],
                    record["exact_proof_inverse_pivot_columns"],
                )
            for region_name, structure_region in structure[
                    "regions"].items():
                region = level["polynomial_regions"][region_name]
                self.assertEqual(
                    len(region["simplices"]),
                    structure_region["simplex_count"],
                )
                vertex_slack = region["vertex_region_certificate"][
                    "minimum_strict_vertex_slack_lower_bound"]
                self.assertGreater(
                    Fraction(int(vertex_slack["numerator"]),
                             int(vertex_slack["denominator"])),
                    0,
                )
                self.assertEqual(
                    sum(item["structural_zero_count"]
                        for item in region["simplices"]),
                    structure_region["total_structural_zero_count"],
                )
                for item in region["simplices"]:
                    determinant = item["affine_determinant"]
                    lower = determinant["lower"]
                    upper = determinant["upper"]
                    lower = Fraction(int(lower["numerator"]),
                                     int(lower["denominator"]))
                    upper = Fraction(int(upper["numerator"]),
                                     int(upper["denominator"]))
                    self.assertFalse(lower <= 0 <= upper)
                    minimum = item["minimum_positive_lower_bound"]
                    self.assertGreater(
                        Fraction(int(minimum["numerator"]),
                                 int(minimum["denominator"])),
                        0,
                    )

            chain = level["triangulation_chain_certificate"]
            self.assertEqual(chain["simplex_count"], 64)
            self.assertEqual(chain["distinct_vertex_pairs_certified"], 378)
            self.assertEqual(
                chain["internal_facets_paired_with_opposite_orientation"],
                112,
            )
            self.assertEqual(chain["unpaired_outer_boundary_facets"], 224)
            self.assertEqual(
                chain["outer_boundary_facets_by_hyperplane"],
                {label: 32 for label in (
                    "x=0", "x=y", "y=z", "z=w", "w=u", "u=v", "v=1")},
            )
            volume = chain["normalized_oriented_volume"]
            lower = volume["lower"]
            upper = volume["upper"]
            lower = Fraction(int(lower["numerator"]),
                             int(lower["denominator"]))
            upper = Fraction(int(upper["numerator"]),
                             int(upper["denominator"]))
            self.assertLessEqual(lower, 1)
            self.assertGreaterEqual(upper, 1)
            self.assertGreater(lower, Fraction(1, 2))
            self.assertLess(upper, Fraction(3, 2))
            self.assertEqual(chain["integer_relative_chain_degree"], 1)

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
