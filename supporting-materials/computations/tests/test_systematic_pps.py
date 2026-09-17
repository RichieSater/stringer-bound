from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = (
    ROOT / "computations" / "certificates" / "systematic-pps-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from systematic_pps import (  # noqa: E402
    build_certificate,
    completion_bound,
    completion_bound_from_observed,
    exhaustive_binary_census,
    finite_population_target,
    hybrid_design_bound,
    markov_design_bound,
    periodic_population,
    phase_sample_mean,
    systematic_mean_expectation,
    systematic_phases,
)


class SystematicPpsTests(unittest.TestCase):
    def test_phase_enumerator_is_exactly_design_unbiased(self) -> None:
        fixtures = (
            ((1, 2, 1), (0, Fraction(1, 2), 1), 2),
            ((1, 3, 2, 4), (Fraction(1, 5), 0, 1, Fraction(2, 3)), 3),
            ((7, 1, 2), (Fraction(1, 7), 1, 0), 4),
        )
        for weights, taints, sample_size in fixtures:
            with self.subTest(weights=weights, sample_size=sample_size):
                phases = systematic_phases(weights, taints, sample_size)
                self.assertEqual(
                    systematic_mean_expectation(phases),
                    finite_population_target(weights, taints),
                )
                self.assertEqual(
                    sum((phase.probability for phase in phases), Fraction(0)),
                    1,
                )

    def test_periodic_family_has_only_two_sample_vectors(self) -> None:
        weights, taints = periodic_population(2, Fraction(4, 5))
        phases = systematic_phases(weights, taints, 2)
        self.assertEqual(len(phases), 2)
        by_taints = {phase.taints: phase.probability for phase in phases}
        self.assertEqual(by_taints[(Fraction(1), Fraction(1))], Fraction(4, 5))
        self.assertEqual(by_taints[(Fraction(0), Fraction(0))], Fraction(1, 5))
        self.assertEqual(finite_population_target(weights, taints), Fraction(4, 5))

    def test_markov_and_hybrid_bounds_cover_each_fixture(self) -> None:
        fixtures = (
            ((1, 2, 1), (0, Fraction(1, 2), 1), 2),
            ((1, 3, 2, 4), (Fraction(1, 5), 0, 1, Fraction(2, 3)), 3),
            ((7, 1, 2), (Fraction(1, 7), 1, 0), 4),
        )
        for weights, taints, sample_size in fixtures:
            phases = systematic_phases(weights, taints, sample_size)
            target = finite_population_target(weights, taints)
            for alpha in (Fraction(1, 10), Fraction(1, 20), Fraction(1, 100)):
                markov_failure = sum(
                    (
                        phase.probability
                        for phase in phases
                        if markov_design_bound(phase.taints, alpha) < target
                    ),
                    Fraction(0),
                )
                hybrid_failure = sum(
                    (
                        phase.probability
                        for phase in phases
                        if hybrid_design_bound(weights, taints, phase, alpha) < target
                    ),
                    Fraction(0),
                )
                self.assertLessEqual(markov_failure, alpha)
                self.assertLessEqual(hybrid_failure, alpha)
                for phase in phases:
                    self.assertGreaterEqual(
                        completion_bound(weights, taints, phase.item_indices),
                        target,
                    )

    def test_sample_mean_can_be_highly_dependent_but_remains_unbiased(self) -> None:
        weights, taints = periodic_population(100, Fraction(1, 10))
        phases = systematic_phases(weights, taints, 100)
        distribution = {
            phase_sample_mean(phase): phase.probability for phase in phases
        }
        self.assertEqual(
            distribution,
            {Fraction(1): Fraction(1, 10), Fraction(0): Fraction(9, 10)},
        )
        self.assertEqual(systematic_mean_expectation(phases), Fraction(1, 10))

        low_phase = next(phase for phase in phases if not any(phase.taints))
        observed = {index: Fraction(0) for index in set(low_phase.item_indices)}
        self.assertEqual(
            completion_bound_from_observed(weights, observed), Fraction(1, 10)
        )

    def test_small_binary_census_finds_the_nine_unit_witness(self) -> None:
        census = exhaustive_binary_census()
        witness = census["first_failure"]
        self.assertEqual(census["population_design_pairs_checked"], 2188)
        self.assertEqual(census["phase_comparisons_checked"], 4444)
        self.assertEqual(witness["total_units"], 9)
        self.assertEqual(witness["sample_size"], 3)
        self.assertEqual(witness["systematic_interval_units"], 3)
        self.assertEqual(witness["ordered_taints"], [1, 1, 0] * 3)
        self.assertEqual(witness["sample_one_counts_by_start"], [3, 3, 0])
        self.assertEqual(witness["coverage"]["numerator"], "2")
        self.assertEqual(witness["coverage"]["denominator"], "3")

    def test_certificate_matches_generator(self) -> None:
        self.assertEqual(json.loads(CERTIFICATE.read_text()), build_certificate())


if __name__ == "__main__":
    unittest.main()
