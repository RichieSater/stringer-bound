"""Mutation tests for the fail-closed public-corpus policy check."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from check_public_corpus import (  # noqa: E402
    PRINCIPAL_MANUSCRIPT,
    _disclosure_violations,
    decode_public_text,
    scan_path_name,
    scan_repository,
    scan_text,
)


class PublicCorpusPolicyTests(unittest.TestCase):
    def test_current_tracked_corpus_passes(self):
        self.assertEqual(scan_repository(), [])

    def test_process_language_mutations_are_rejected(self):
        words = " ".join
        mutations = (
            words(("An", "independent", "human", "proof", "review", "is", "required")),
            words(("Use", "this", "review", "protocol", "before", "release")),
            "Suggested " + "sign" + "-off block",
            words(("The", "draft", "has", "not", "undergone", "peer", "review")),
            words(("Add", "a", "reader", "test", "with", "named", "participants")),
            words(("Reviewer", "assignment", "table")),
            words(("Publish", "a", "role", "matrix")),
            words(("Publish", "a", "reviewer", "roster")),
            words(("Status:", "pending", "evaluation")),
            words(("Status:", "awaiting", "validation")),
            words(("This", "is", "an", "approval", "gate")),
            words(("A", "statistician", "must", "approve", "the", "release")),
            words(("The", "draft", "must", "be", "validated", "by", "a", "reader")),
            words(
                ("Submission", "is", "contingent", "on", "the", "editor's", "approval")
            ),
            words(("Records", "of", "two", "adversarial", "review", "passes")),
            "This is an " + "AI" + "-assisted review",
            words(("The", "draft", "needs", "human", "validation")),
        )
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                self.assertTrue(scan_text(Path("README.md"), mutation))

    def test_internal_process_artifact_name_mutations_are_rejected(self):
        names = (
            "supporting-materials/audit/HUMAN-" + "REVIEW-PACKET.md",
            "supporting-materials/audit/REVIEW-" + "N2.md",
            "supporting-materials/audit/REVIEW-" + "PROTOCOL.md",
            "supporting-materials/audit/SIGN-" + "OFF.md",
            "supporting-materials/audit/ROLE-" + "MATRIX.md",
            "supporting-materials/audit/REVIEWER-" + "ROSTER.md",
            "supporting-materials/audit/REVIEWER_" + "ASSIGNMENT.md",
            "supporting-materials/audit/PENDING-" + "EVALUATION.md",
        )
        for name in names:
            with self.subTest(name=name):
                self.assertTrue(scan_path_name(Path(name)))

    def test_duplicate_disclosure_mutations_are_rejected(self):
        principal = "AI " + "disclosure Anthropic " + "Claude OpenAI " + "Codex"
        mutations = (
            "Generative " + "AI was used to draft this note.",
            "A large language " + "model prepared this appendix.",
            "Chat" + "GPT generated this README.",
            "Google " + "Gemini assisted with this supplement.",
        )
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                violations = _disclosure_violations(
                    {
                        PRINCIPAL_MANUSCRIPT: principal,
                        Path("README.md"): mutation,
                    }
                )
                self.assertTrue(violations)

    def test_unknown_binary_format_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "public-note.md"
            with self.assertRaises(UnicodeError):
                decode_public_text(path, b"\xff\xfe\x00")


if __name__ == "__main__":
    unittest.main()
