"""Build the solution index for the reusable simplex-cap challenge.

The challenge JSON contains only the theorem template, instance scopes, and
verification obligations.  This module keeps solution paths, commands, and
digests in a separate index.  The index is an integrity/composition layer;
the instance-specific commands remain responsible for the mathematical sign
checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CHALLENGE = (
    ROOT / "supporting-materials" / "computations" / "challenges"
    / "simplex-cap-challenge.json"
)

INSTANCE_SOURCES = {
    "binomial-n3-confidence-range": {
        "status": "proved",
        "sources": [
            ("supporting-materials/theory/N3-CONVENTIONAL.md",
             "written cap reduction and uniform-alpha proof"),
            ("supporting-materials/computations/certificates/n3-gaffke-bernstein-formulas.json",
             "derived middle-region Bernstein formulas"),
            ("supporting-materials/computations/certificates/n3-gaffke-certificate.json",
             "exact 414-cell range sign data"),
            ("supporting-materials/computations/python/derive_n3_bernstein_formulas.py",
             "symbolic formula regeneration"),
            ("supporting-materials/computations/python/n3_gaffke_certificate.py",
             "exact range-certificate regeneration"),
        ],
        "commands": [
            "make n3-formula-check",
            "make n3-certificate-check",
        ],
    },
    "binomial-n4-confidence-range": {
        "status": "unresolved",
        "proof_obligation": 'The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.',
        "sources": [
            ("supporting-materials/theory/FIVE-COORDINATE-VERTEX-BARRIER.md",
             "local inequalities and incomplete global face implication"),
            ("supporting-materials/theory/FIVE-COORDINATE-MIXED-SECTION.md",
             "distinct-knot and unique-top local reductions"),
            ("supporting-materials/computations/certificates/five-coordinate-mixed-section-certificate.json",
             "exact mixed-section coefficient-sign data"),
            ("supporting-materials/computations/python/five_coordinate_mixed_section.py",
             "exact mixed-section identity and sign regeneration"),
            ("supporting-materials/computations/python/five_coordinate_mixed_witness.py",
             "index-only sparse exact-system witness"),
            ("supporting-materials/computations/certificates/five-coordinate-block-faces-certificate.json",
             "exact repeated-top coefficient-sign data"),
            ("supporting-materials/computations/python/five_coordinate_block_faces.py",
             "exact repeated-top identity and sign regeneration"),
            ("supporting-materials/computations/python/five_coordinate_block_witness.py",
             "index-only repeated-top exact-system witnesses"),
            ("supporting-materials/theory/N4-MONOTONE-RANGE.md",
             "proved weight-order analysis and conditional coverage implication"),
            ("supporting-materials/computations/certificates/n4-weight-monotonicity-certificate.json",
             "exact resultant, Sturm, and endpoint data"),
            ("supporting-materials/computations/python/n4_weight_monotonicity.py",
             "exact resultant and weight-range regeneration"),
        ],
        "commands": [
            "make five-coordinate-mixed-section-check",
            "make five-coordinate-block-faces-check",
            "make n4-weight-monotonicity-check",
        ],
    },
    "binomial-n5-confidence-range": {
        "status": "unresolved",
        "proof_obligation": 'The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.',
        "sources": [
            ("supporting-materials/theory/ADJACENT-BLOCK-HIERARCHY.md",
             "ordered-composition reduction and exhaustive family count"),
            ("supporting-materials/theory/SIX-COORDINATE-REPEATED-TOP.md",
             "all six repeated-top local barriers"),
            ("supporting-materials/computations/certificates/six-coordinate-repeated-top-certificate.json",
             "exact repeated-top coefficient-sign data"),
            ("supporting-materials/computations/python/six_coordinate_repeated_top.py",
             "exact repeated-top identity and sign regeneration"),
            ("supporting-materials/computations/python/six_coordinate_face32.py",
             "exact final repeated-top projective certificate"),
            ("supporting-materials/theory/SIX-COORDINATE-UNIQUE-TOP.md",
             "unique-top L=4 projective reduction"),
            ("supporting-materials/computations/certificates/six-coordinate-unique-top-certificate.json",
             "exact unique-top L=4 coefficient-sign data"),
            ("supporting-materials/computations/python/six_coordinate_unique_top.py",
             "exact unique-top L=4 sign regeneration"),
            ("supporting-materials/theory/SIX-COORDINATE-UNIQUE-TOP-L2-L3.md",
             "sharp minorant and final unique-top projective reductions"),
            ("supporting-materials/computations/certificates/six-coordinate-unique-top-l2-l3-certificate.json",
             "exact and rigorously enclosed L=2/L=3 sign data"),
            ("supporting-materials/computations/python/six_coordinate_unique_top_l2_l3.py",
             "L=2/L=3 symbolic and certificate regeneration"),
            ("supporting-materials/computations/python/certified_binary64_bernstein.py",
             "exact forward-error verification of the large L=2 chart"),
            ("supporting-materials/theory/N5-MONOTONE-RANGE.md",
             "proved weight-order analysis and conditional coverage implication"),
            ("supporting-materials/computations/certificates/n5-weight-monotonicity-certificate.json",
             "exact resultant, Sturm, and endpoint data"),
            ("supporting-materials/computations/python/n5_weight_monotonicity.py",
             "exact resultant and weight-range regeneration"),
        ],
        "commands": [
            "make six-coordinate-repeated-top-check",
            "make six-coordinate-unique-top-check",
            "make six-coordinate-unique-top-l2-l3-check",
            "make n5-weight-monotonicity-check",
        ],
    },
    "binomial-n5-conventional-levels": {
        "status": "proved",
        "sources": [
            ("supporting-materials/theory/N5-CONVENTIONAL.md",
             "written five-region cap reduction"),
            ("supporting-materials/computations/certificates/n5-gaffke-bernstein-structure.json",
             "exact residual, face, and triangulation structure"),
            ("supporting-materials/computations/certificates/n5-gaffke-certificate.json",
             "directed-dyadic coefficient-sign data"),
            ("supporting-materials/computations/python/derive_n5_bernstein_structure.py",
             "symbolic structure regeneration"),
            ("supporting-materials/computations/python/n5_gaffke_certificate.py",
             "directed-dyadic sign regeneration"),
        ],
        "commands": [
            "make n5-structure-check",
            "make n5-certificate-check",
        ],
    },
    "binomial-n6-conventional-levels": {
        "status": "proved",
        "sources": [
            ("supporting-materials/theory/N6-CONVENTIONAL.md",
             "written six-region cap reduction"),
            ("supporting-materials/computations/certificates/n6-gaffke-bernstein-structure.json",
             "exact residual, face, and triangulation structure"),
            ("supporting-materials/computations/certificates/n6-gaffke-certificate.json",
             "directed-dyadic coefficient-sign data"),
            ("supporting-materials/computations/python/derive_n6_bernstein_structure.py",
             "symbolic structure and face-ideal regeneration"),
            ("supporting-materials/computations/python/n6_gaffke_certificate.py",
             "directed-dyadic sign regeneration"),
        ],
        "commands": [
            "make n6-structure-check",
            "make n6-certificate-check",
        ],
    },
    "binomial-n7-conventional-levels": {
        "status": "proved",
        "sources": [
            ("supporting-materials/theory/N7-CONVENTIONAL.md",
             "written seven-region cap reduction"),
            ("supporting-materials/computations/certificates/n7-gaffke-bernstein-structure.json.gz",
             "exact compressed residual and face structure"),
            ("supporting-materials/computations/certificates/n7-gaffke-certificate.json",
             "rigorous Arb sign summary"),
            ("supporting-materials/computations/python/n7_gaffke_structure_data.py",
             "fixed triangulation and face conditions"),
            ("supporting-materials/computations/python/derive_n7_bernstein_structure.py",
             "symbolic structure and face-proof regeneration"),
            ("supporting-materials/computations/python/n7_gaffke_certificate.py",
             "outward-rounded Arb sign regeneration"),
        ],
        "commands": [
            "make n7-structure-data-check",
            "make n7-certificate-check",
        ],
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _walk_keys(value):
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from _walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_keys(child)


def _validate_challenge(challenge):
    if challenge.get("schema_version") != 1:
        raise ValueError("unexpected simplex-cap challenge schema")
    if challenge.get("interface") != "ordered-simplex-cap-comparison":
        raise ValueError("unexpected simplex-cap interface name")
    forbidden_solution_keys = {
        "solution", "solutions", "evidence", "source", "sources",
        "command", "commands", "sha256", "digest", "digests",
    }
    present = forbidden_solution_keys.intersection(_walk_keys(challenge))
    if present:
        raise ValueError(
            "challenge statement contains solution-layer keys: "
            + ", ".join(sorted(present)))
    obligations = {
        item["id"] for item in challenge.get("verification_obligations", [])
    }
    expected_obligations = {
        "factor_scope", "region_partition", "cap_identity",
        "triangulation", "structural_zeros", "coefficient_signs",
        "degenerate_knots", "artifact_closure",
    }
    if obligations != expected_obligations:
        raise ValueError("challenge verification obligations changed")
    ids = [item["id"] for item in challenge.get("instances", [])]
    if len(ids) != len(set(ids)):
        raise ValueError("challenge instance ids are not unique")
    if set(ids) != set(INSTANCE_SOURCES):
        raise ValueError("challenge and solution-index instances differ")


def build_solution_index():
    challenge_raw = CHALLENGE.read_bytes()
    challenge = json.loads(challenge_raw)
    _validate_challenge(challenge)
    instances = []
    for challenge_instance in challenge["instances"]:
        instance_id = challenge_instance["id"]
        specification = INSTANCE_SOURCES[instance_id]
        sources = []
        for relative_path, role in specification["sources"]:
            path = ROOT / relative_path
            if not path.is_file():
                raise FileNotFoundError(path)
            sources.append({
                "path": relative_path,
                "role": role,
                "sha256": _sha256(path),
            })
        instances.append({
            "id": instance_id,
            "mathematical_status": specification["status"],
            "proof_obligation": specification.get("proof_obligation"),
            "scope_copy": {
                "sample_size": challenge_instance["sample_size"],
                "alpha_domain": challenge_instance["alpha_domain"],
            },
            "sources": sources,
            "verification_commands": specification["commands"],
        })
    return {
        "schema_version": 1,
        "interface": challenge["interface"],
        "challenge_path": str(CHALLENGE.relative_to(ROOT)),
        "challenge_sha256": hashlib.sha256(challenge_raw).hexdigest(),
        "composition_status": (
            "The challenge is solution-free; this index binds each declared "
            "instance to its written reduction, deterministic data, checker, "
            "and regeneration commands."
        ),
        "limitation": (
            "Index regeneration checks scope and artifact integrity. The "
            "listed commands perform symbolic and numerical sign checks. "
            "An unresolved global implication is not proved by artifact integrity."
        ),
        "instances": instances,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        help="write the machine-readable solution index")
    args = parser.parse_args(argv)
    index = build_solution_index()
    rendered = json.dumps(index, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    print(
        "simplex-cap interface: solution-free challenge and "
        f"{len(index['instances'])} instance bindings checked (status explicit)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
