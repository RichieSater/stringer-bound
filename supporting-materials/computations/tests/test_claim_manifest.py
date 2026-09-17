"""Document-anchor and mathematical-status mutations for the claim manifest."""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from validate_claim_manifest import MANIFEST, validate  # noqa: E402


class ClaimManifestTests(unittest.TestCase):
    def test_repository_links_and_unresolved_scopes(self):
        self.assertEqual(validate(), [])
        claims = {
            claim["id"]: claim
            for claim in json.loads(MANIFEST.read_text())["claims"]
        }
        for name in (
            "five-coordinate-monotone-cap-theorem",
            "six-coordinate-monotone-cap-theorem",
            "binomial-n4-expanded-conservatism",
            "binomial-n5-expanded-conservatism",
        ):
            self.assertEqual(claims[name]["classification"], "unresolved_claim")
            self.assertIn("active-prefix", claims[name]["proof_obligation"])
            self.assertIn("repeated-lowest", claims[name]["proof_obligation"])
        self.assertEqual(
            claims["binomial-n5-monotone-weight-range"]["classification"],
            "proved_theorem",
        )

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        (self.root / "main.tex").write_text(r"\label{thm:main}")
        (self.root / "note.md").write_text("# Local inequality\n## Exact signs\n")
        self.payload = {
            "manuscript": "main.tex",
            "claims": [{
                "id": "local-signs", "classification": "exact_computational_certificate",
                "statement": "The listed polynomial coefficients are nonnegative.",
                "manuscript": "note.md", "manuscript_labels": [],
                "source_markers": ["## Exact signs"],
                "evidence": [{"path": "note.md", "role": "exact identities",
                              "trust": "proof_essential"}],
                "verification_commands": [],
                "limitations": "Only the listed polynomials are covered.",
            }],
        }

    def test_note_anchor_passes(self):
        self.assertEqual(validate(self.payload, self.root), [])

    def test_deleted_and_partial_heading_mutations_are_rejected(self):
        for marker in ("## Missing signs", "Exact signs", ""):
            mutated = copy.deepcopy(self.payload)
            mutated["claims"][0]["source_markers"] = [marker]
            self.assertTrue(validate(mutated, self.root))

    def test_unanchored_claim_is_rejected(self):
        self.payload["claims"][0]["source_markers"] = []
        self.assertTrue(validate(self.payload, self.root))

    def test_unresolved_status_requires_exact_obligation(self):
        claim = self.payload["claims"][0]
        for status in ("unresolved_claim", "conjectural_claim"):
            claim["classification"] = status
            claim.pop("proof_obligation", None)
            self.assertTrue(validate(self.payload, self.root))
            claim["proof_obligation"] = (
                "Establish the active-prefix inequality on repeated-lowest faces."
            )
            self.assertEqual(validate(self.payload, self.root), [])

    def test_missing_evidence_and_dead_tex_label_are_rejected(self):
        for field, value in (("manuscript_labels", ["thm:removed"]),
                             ("evidence", [{"path": "missing.md", "role": "signs",
                                            "trust": "proof_essential"}])):
            mutated = copy.deepcopy(self.payload)
            mutated["claims"][0][field] = value
            self.assertTrue(validate(mutated, self.root))


if __name__ == "__main__":
    unittest.main()
