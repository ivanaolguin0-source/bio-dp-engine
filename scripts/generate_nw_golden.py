"""Regenerate the Needleman-Wunsch golden fixtures.

Inputs are chosen by hand to exercise one behaviour each. Expected values
come from the exhaustive oracle (engines/python/tests/oracles.py), never from
the engine under test, and are cross-checked against Biopython's
PairwiseAligner when it is installed. Any disagreement aborts the run.

Usage:
    python scripts/generate_nw_golden.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "engines" / "python" / "tests"))

from oracles import exhaustive_global, parse_ncbi_matrix  # noqa: E402

OUT = ROOT / "fixtures" / "golden" / "needleman_wunsch"
MAX_LISTED = 50

CASES: list[dict[str, Any]] = [
    {
        "id": "nw-001-both-empty",
        "description": "Empty against empty: score 0 and exactly one (empty) alignment.",
        "seq_a": "", "seq_b": "",
        "substitution": {"type": "match_mismatch", "match": 1, "mismatch": -1, "alphabet": "ACGT"},
        "gap": {"open": -1, "extend": -1},
    },
    {
        "id": "nw-002-one-empty-affine",
        "description": "Empty against ACG is one gap of length 3: open + 2 * extend.",
        "seq_a": "", "seq_b": "ACG",
        "substitution": {"type": "match_mismatch", "match": 1, "mismatch": -1, "alphabet": "ACGT"},
        "gap": {"open": -2, "extend": -1},
    },
    {
        "id": "nw-003-homopolymer-ties",
        "description": "AAA vs AA: the single gap can sit in any of three places, giving three co-optimal alignments.",
        "seq_a": "AAA", "seq_b": "AA",
        "substitution": {"type": "match_mismatch", "match": 1, "mismatch": -1, "alphabet": "ACGT"},
        "gap": {"open": -1, "extend": -1},
    },
    {
        "id": "nw-004-adjacent-opposite-gaps",
        "description": "With a harsh mismatch, A vs C is best as two adjacent gaps of different kinds, in either order.",
        "seq_a": "A", "seq_b": "C",
        "substitution": {"type": "match_mismatch", "match": 1, "mismatch": -10, "alphabet": "ACGT"},
        "gap": {"open": -1, "extend": -1},
    },
    {
        "id": "nw-005-textbook-gattaca",
        "description": "The GATTACA / GCATGCU teaching example with +1 / -1 / -1 scoring.",
        "seq_a": "GATTACA", "seq_b": "GCATGCU",
        "substitution": {"type": "match_mismatch", "match": 1, "mismatch": -1, "alphabet": "ACGTU"},
        "gap": {"open": -1, "extend": -1},
    },
    {
        "id": "nw-006-affine-merges-gap",
        "description": "Affine scoring should delete the TTTT block as one gap rather than scattering it.",
        "seq_a": "ACGTTTTACG", "seq_b": "ACGACG",
        "substitution": {"type": "match_mismatch", "match": 2, "mismatch": -3, "alphabet": "ACGT"},
        "gap": {"open": -5, "extend": -1},
    },
    {
        "id": "nw-007-durbin-blosum50-linear",
        "description": "HEAGAWGHEE / PAWHEAE, the worked example in Durbin et al. (1998), with BLOSUM50 and linear gap -8.",
        "seq_a": "HEAGAWGHEE", "seq_b": "PAWHEAE",
        "substitution": {"type": "named", "name": "BLOSUM50"},
        "gap": {"open": -8, "extend": -8},
    },
    {
        "id": "nw-008-blosum62-affine",
        "description": "Same peptides with BLOSUM62 and BLAST-style affine gaps (open -11, extend -1).",
        "seq_a": "HEAGAWGHEE", "seq_b": "PAWHEAE",
        "substitution": {"type": "named", "name": "BLOSUM62"},
        "gap": {"open": -11, "extend": -1},
    },
]


def score_function(spec: dict[str, Any]) -> Any:
    if spec["type"] == "match_mismatch":
        match, mismatch = spec["match"], spec["mismatch"]
        return lambda x, y: match if x == y else mismatch
    table = parse_ncbi_matrix((ROOT / "fixtures" / "matrices" / f"{spec['name']}.txt").read_text())
    return lambda x, y: table[(x, y)]


def biopython_check(case: dict[str, Any], score: int, count: int) -> str | None:
    """Return a provenance string if Biopython agrees, raise if it disagrees, None if unavailable."""
    try:
        import Bio
        from Bio.Align import PairwiseAligner, substitution_matrices
    except ImportError:
        return None
    if not case["seq_a"] or not case["seq_b"]:
        return None  # PairwiseAligner rejects empty sequences
    aligner = PairwiseAligner(mode="global")
    spec = case["substitution"]
    if spec["type"] == "match_mismatch":
        aligner.match_score, aligner.mismatch_score = spec["match"], spec["mismatch"]
    else:
        aligner.substitution_matrix = substitution_matrices.read(
            str(ROOT / "fixtures" / "matrices" / f"{spec['name']}.txt")
        )
    aligner.open_gap_score = case["gap"]["open"]
    aligner.extend_gap_score = case["gap"]["extend"]
    alignments = aligner.align(case["seq_a"], case["seq_b"])
    bp_score, bp_count = int(alignments.score), len(alignments)
    if (bp_score, bp_count) != (score, count):
        raise SystemExit(f"{case['id']}: oracle ({score}, {count}) != Biopython ({bp_score}, {bp_count})")
    return f"Biopython {Bio.__version__} PairwiseAligner"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        oracle = exhaustive_global(
            case["seq_a"], case["seq_b"], score_function(case["substitution"]),
            case["gap"]["open"], case["gap"]["extend"],
        )
        count = len(oracle.alignments)
        listed = sorted(oracle.alignments)
        fixture = {
            "id": case["id"],
            "algorithm": "needleman-wunsch",
            "description": case["description"],
            "input": {
                "seq_a": case["seq_a"],
                "seq_b": case["seq_b"],
                "substitution": case["substitution"],
                "gap": case["gap"],
            },
            "expected": {
                "score": oracle.score,
                "optimal_count": str(count),
                "alignments": [list(pair) for pair in listed] if count <= MAX_LISTED else None,
            },
            "provenance": {
                "generator": "exhaustive enumeration of all alignments (engines/python/tests/oracles.py)",
                "cross_checked_with": biopython_check(case, oracle.score, count),
            },
        }
        path = OUT / f"{case['id']}.json"
        path.write_text(json.dumps(fixture, indent=2) + "\n")
        print(f"{case['id']}: score={oracle.score} optimal={count}")


if __name__ == "__main__":
    main()
