"""Fail-closed policy check for text circulated from this repository."""

from __future__ import annotations

import argparse
import dataclasses
import pathlib
import re
import subprocess
import sys
from collections.abc import Iterable

REPOSITORY_ROOT = pathlib.Path(__file__).resolve().parents[3]
PRINCIPAL_MANUSCRIPT = pathlib.Path("supporting-materials/paper/stringer.tex")
BINARY_SUFFIXES = {".gz", ".jpeg", ".jpg", ".pdf", ".png", ".zip"}


@dataclasses.dataclass(frozen=True)
class Violation:
    path: pathlib.Path
    line: int
    rule: str
    excerpt: str


def _pattern(*parts: str) -> re.Pattern[str]:
    return re.compile("".join(parts), re.IGNORECASE | re.MULTILINE)


# Patterns use regular-expression separators so this checker does not itself
# place the disallowed prose into the circulated text corpus.
PROSE_RULES: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "person-gated-check",
        _pattern(
            r"\b(?:independent\s+)?human\s+(?:proof\s+)?",
            r"(?:review|validation)\b",
        ),
    ),
    (
        "publication-status-process",
        _pattern(
            r"\b(?:peer[-\s]?review(?:ed)?|not\s+independently\s+refereed)",
            r"\b",
        ),
    ),
    ("process-authorization-mark", _pattern(r"\bsign[-\s]?off\b")),
    (
        "process-disposition",
        _pattern(r"\b(?:requested|pending)\s+(?:output:\s*)?dispositions?\b"),
    ),
    ("process-review-protocol", _pattern(r"\breview\s+protocol\b")),
    ("process-review-packet", _pattern(r"\breview\s+packet\b")),
    ("process-review-record", _pattern(r"\breview\s+record\b")),
    ("process-reader-test", _pattern(r"\breader\s+tests?\b")),
    ("process-named-participants", _pattern(r"\bnamed\s+participants?\b")),
    ("process-pending-evaluation", _pattern(r"\bpending\s+evaluation\b")),
    (
        "process-assignment",
        _pattern(
            r"\b(?:reviewer|reader|participant)\s+assignments?\b|",
            r"\brole\s+matri(?:x|ces)\b",
        ),
    ),
    (
        "process-roster",
        _pattern(r"\b(?:reviewer|reader|participant)\s+rosters?\b"),
    ),
    (
        "process-awaiting-evaluation",
        _pattern(r"\bawaiting\s+(?:review|validation|approval)\b"),
    ),
    (
        "process-approval-gate",
        _pattern(
            r"\b(?:pending|required|requires?)\s+(?:final\s+)?approval\b|",
            r"\bapproval\s+(?:gate|required)\b",
        ),
    ),
    (
        "person-prerequisite",
        _pattern(
            r"\b(?:person|researcher|statistician|reviewer|reader|participant|",
            r"author|editor|committee)\b",
            r"[^\n]{0,80}\b(?:must|required\s+to|should)\s+",
            r"(?:approve|validate|sign)\b|",
            r"\b(?:must|should|required\s+to)\s+be\s+",
            r"(?:approved|validated|signed)\s+by\s+(?:a|an|the)\s+",
            r"(?:person|researcher|statistician|reviewer|reader|participant|",
            r"author|editor|committee)\b|",
            r"\b(?:release|submission|readiness|completion)\b[^\n]{0,80}",
            r"\bcontingent\s+on\b[^\n]{0,80}\bapproval\b",
        ),
    ),
    (
        "vague-ai-status",
        _pattern(
            r"\bai[-\s]assisted(?:\s+review)?\b|",
            r"\bneeds\s+human\s+validation\b|",
            r"\bclaimed\s+complete\s+solution\b",
        ),
    ),
    (
        "internal-pass-logistics",
        _pattern(r"\badversarial\s+(?:review\s+)?passes\b"),
    ),
)

FORBIDDEN_PATH_PARTS = (
    "human-" + "review",
    "review-" + "packet",
    "review-" + "protocol",
    "review-" + "record",
    "review-" + "n2",
    "peer-" + "review",
    "sign-" + "off",
    "sign" + "off",
    "role-" + "matrix",
    "reader-" + "test",
    "reviewer-" + "roster",
    "reviewer-" + "assignment",
    "pending-" + "evaluation",
)
PRINCIPAL_AI_MARKERS = (
    "ai dis" + "closure",
    "anthropic " + "claude",
    "openai " + "codex",
)
DUPLICATE_AI_MARKERS = PRINCIPAL_AI_MARKERS + (
    "generative " + "ai",
    "large language " + "model",
    "chat" + "gpt",
    "google " + "gemini",
)


def tracked_paths(root: pathlib.Path = REPOSITORY_ROOT) -> list[pathlib.Path]:
    payload = subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=root, stderr=subprocess.STDOUT
    )
    return [pathlib.Path(item) for item in payload.decode("utf-8").split("\0") if item]


def decode_public_text(path: pathlib.Path, payload: bytes) -> str | None:
    """Decode a tracked public file, rejecting unknown non-UTF-8 formats."""
    if path.suffix.lower() in BINARY_SUFFIXES:
        return None
    if b"\0" in payload:
        raise UnicodeError(f"unclassified binary file: {path}")
    return payload.decode("utf-8")


def scan_text(path: pathlib.Path, text: str) -> list[Violation]:
    violations = []
    for rule, regex in PROSE_RULES:
        for match in regex.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            excerpt = text.splitlines()[line - 1].strip()
            violations.append(Violation(path, line, rule, excerpt))
    return violations


def scan_path_name(path: pathlib.Path) -> list[Violation]:
    normalized = path.as_posix().lower().replace("_", "-").replace(" ", "-")
    return [
        Violation(path, 1, "internal-process-artifact-path", normalized)
        for part in FORBIDDEN_PATH_PARTS
        if part in normalized
    ]


def _disclosure_violations(text_by_path: dict[pathlib.Path, str]) -> list[Violation]:
    violations = []
    principal = text_by_path.get(PRINCIPAL_MANUSCRIPT, "")
    for marker in PRINCIPAL_AI_MARKERS:
        count = principal.lower().count(marker)
        if count != 1:
            violations.append(
                Violation(
                    PRINCIPAL_MANUSCRIPT,
                    1,
                    "principal-disclosure-count",
                    f"{marker}: expected one occurrence, found {count}",
                )
            )
    for path, text in text_by_path.items():
        if path == PRINCIPAL_MANUSCRIPT:
            continue
        lowered = text.lower()
        for marker in DUPLICATE_AI_MARKERS:
            if marker in lowered:
                line = lowered.count("\n", 0, lowered.index(marker)) + 1
                violations.append(
                    Violation(path, line, "duplicate-ai-disclosure", marker)
                )
    return violations


def scan_repository(
    root: pathlib.Path = REPOSITORY_ROOT,
    paths: Iterable[pathlib.Path] | None = None,
) -> list[Violation]:
    """Scan every tracked public text file; unknown encodings fail closed."""
    violations = []
    text_by_path: dict[pathlib.Path, str] = {}
    for relative in tracked_paths(root) if paths is None else paths:
        absolute = root / relative
        if not absolute.is_file():
            continue
        violations.extend(scan_path_name(relative))
        try:
            text = decode_public_text(relative, absolute.read_bytes())
        except UnicodeError as error:
            violations.append(
                Violation(relative, 1, "unclassified-nontext-file", str(error))
            )
            continue
        if text is None:
            continue
        text_by_path[relative] = text
        violations.extend(scan_text(relative, text))
    violations.extend(_disclosure_violations(text_by_path))
    return sorted(
        violations, key=lambda item: (item.path.as_posix(), item.line, item.rule)
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.parse_args(argv)
    violations = scan_repository()
    for item in violations:
        print(
            f"{item.path}:{item.line}: {item.rule}: {item.excerpt}",
            file=sys.stderr,
        )
    if violations:
        return 1
    print("tracked public text satisfies the repository corpus policy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
