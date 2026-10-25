"""Regenerate the Nussinov golden fixtures.

Expected values come from exhaustive enumeration of every secondary structure
(engines/python/tests/oracles.py), built from the first base of each interval.
The engine decomposes on the last base, so the two share no recurrence. The
classic four-case recurrence is also run as a second check on the score.

Usage:
    python scripts/generate_nussinov_golden.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "engines" / "python" / "tests"))

from oracles import classic_nussinov_score, max_pair_structures  # noqa: E402

OUT = ROOT / "fixtures" / "golden" / "nussinov"
MAX_LISTED = 50

CASES: list[dict[str, Any]] = [
    {"id": "nus-001-empty", "description": "Empty sequence: no pairs, one (empty) structure.",
     "sequence": "", "min_loop": 3, "wobble": True},
    {"id": "nus-002-too-short-to-pair", "description": "Four bases cannot enclose a 3-base loop, so nothing pairs.",
     "sequence": "GAUC", "min_loop": 3, "wobble": True},
    {"id": "nus-003-simple-hairpin", "description": "A three-pair stem closing a three-base loop.",
     "sequence": "GGGAAAUCC", "min_loop": 3, "wobble": True},
    {"id": "nus-004-no-loop-constraint", "description": "The same sequence with min_loop 0: more ways to place the same number of pairs.",
     "sequence": "GGGAAAUCC", "min_loop": 0, "wobble": True},
    {"id": "nus-005-wobble-on", "description": "G-U wobble pairs allowed.",
     "sequence": "GGGAAAUUU", "min_loop": 3, "wobble": True},
    {"id": "nus-006-wobble-off", "description": "The same sequence with Watson-Crick pairs only.",
     "sequence": "GGGAAAUUU", "min_loop": 3, "wobble": False},
    {"id": "nus-007-two-hairpins", "description": "Two separate stem-loops side by side (a bifurcation).",
     "sequence": "GGGAAACCCAGGGAAACCC", "min_loop": 3, "wobble": False},
    {"id": "nus-008-homopolymer-stem", "description": "A run of A's against a run of U's: the stem can register in several ways.",
     "sequence": "AAAAAAUUUUUU", "min_loop": 3, "wobble": True},
    {"id": "nus-009-many-ties", "description": "Four C's can pair with the G's on either side, so dozens of structures tie.",
     "sequence": "GGGGAAAACCCCAAAAGGGG", "min_loop": 3, "wobble": True},
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        best, structures = max_pair_structures(case["sequence"], case["min_loop"], case["wobble"])
        classic = classic_nussinov_score(case["sequence"], case["min_loop"], case["wobble"])
        if classic != best:
            raise SystemExit(f"{case['id']}: enumeration says {best}, classic recurrence says {classic}")
        listed = sorted([list(p) for p in s] for s in structures)
        fixture = {
            "id": case["id"],
            "algorithm": "nussinov",
            "description": case["description"],
            "input": {"sequence": case["sequence"], "min_loop": case["min_loop"], "wobble": case["wobble"]},
            "expected": {
                "score": best,
                "optimal_count": str(len(structures)),
                "structures": listed if len(structures) <= MAX_LISTED else None,
            },
            "provenance": {
                "generator": "exhaustive enumeration of all secondary structures, first-base decomposition (engines/python/tests/oracles.py)",
                "cross_checked_with": "classic four-case Nussinov recurrence (score)",
            },
        }
        (OUT / f"{case['id']}.json").write_text(json.dumps(fixture, indent=2) + "\n")
        print(f"{case['id']}: max pairs={best} optimal structures={len(structures)}")


if __name__ == "__main__":
    main()
