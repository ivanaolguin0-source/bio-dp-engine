"""Shared test utilities: repository paths and fixture loading."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from biodp.scoring import GapPenalty, SubstitutionMatrix

REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURES = REPO_ROOT / "fixtures"
SCHEMA_PATH = REPO_ROOT / "contracts" / "schemas" / "dptrace.schema.json"
NW_GOLDEN = sorted((FIXTURES / "golden" / "needleman_wunsch").glob("*.json"))


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def named_matrix(name: str) -> SubstitutionMatrix:
    return SubstitutionMatrix.from_ncbi_file(FIXTURES / "matrices" / f"{name}.txt")


def scoring_from_fixture(fixture: dict[str, Any]) -> tuple[SubstitutionMatrix, GapPenalty]:
    spec = fixture["input"]["substitution"]
    if spec["type"] == "match_mismatch":
        matrix = SubstitutionMatrix.match_mismatch(spec["match"], spec["mismatch"], spec["alphabet"])
    else:
        matrix = named_matrix(spec["name"])
    gap = GapPenalty(fixture["input"]["gap"]["open"], fixture["input"]["gap"]["extend"])
    return matrix, gap
