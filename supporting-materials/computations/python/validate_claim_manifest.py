"""Validate the machine-readable map from manuscript claims to evidence."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST = ROOT / "supporting-materials" / "claim-evidence.json"

CLASSIFICATIONS = {
    "proved_theorem",
    "proved_literature_correction",
    "exact_computational_certificate",
    "numerical_search_evidence",
    "date_stamped_literature_status",
}
TRUST_LEVELS = {
    "proof_essential",
    "corroboration_only",
    "review_record",
    "heuristic_search",
    "literature_review",
}


def validate() -> list[str]:
    payload = json.loads(MANIFEST.read_text())
    errors: list[str] = []
    claims = payload.get("claims", [])
    ids = [claim.get("id") for claim in claims]
    if not claims:
        errors.append("manifest contains no claims")
    if len(ids) != len(set(ids)):
        errors.append("claim ids are not unique")

    manuscript_path = ROOT / payload.get("manuscript", "")
    if not manuscript_path.is_file():
        errors.append(f"missing manuscript: {manuscript_path}")
        manuscript = ""
    else:
        manuscript = manuscript_path.read_text()

    for claim in claims:
        claim_id = claim.get("id", "<missing-id>")
        if claim.get("classification") not in CLASSIFICATIONS:
            errors.append(f"{claim_id}: invalid classification")
        if not claim.get("statement"):
            errors.append(f"{claim_id}: missing statement")
        if "limitations" not in claim or not claim["limitations"]:
            errors.append(f"{claim_id}: missing limitations")
        if "verification_commands" not in claim:
            errors.append(f"{claim_id}: missing verification_commands")

        for label in claim.get("manuscript_labels", []):
            if (f"\\label{{{label}}}" not in manuscript
                    and f"\\tag{{{label}}}" not in manuscript):
                errors.append(f"{claim_id}: manuscript label not found: {label}")

        evidence = claim.get("evidence", [])
        if not evidence:
            errors.append(f"{claim_id}: no evidence entries")
        for item in evidence:
            path = ROOT / item.get("path", "")
            if not path.is_file():
                errors.append(f"{claim_id}: missing evidence file: {path}")
            if item.get("trust") not in TRUST_LEVELS:
                errors.append(f"{claim_id}: invalid trust level")
            if not item.get("role"):
                errors.append(f"{claim_id}: evidence role is empty")
    return errors


def main() -> int:
    errors = validate()
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    payload = json.loads(MANIFEST.read_text())
    print(f"claim manifest: {len(payload['claims'])} claims, all links valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
