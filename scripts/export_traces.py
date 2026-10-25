"""Export DPTrace files for the web visualiser.

Runs the Python engine on every golden fixture plus a few longer demo inputs,
and writes one DPTrace JSON per case into web/public/traces/, together with an
index.json manifest the frontend reads to populate its picker.

The frontend never computes DP; it only replays these files. Regenerate after
any engine change:

    python scripts/export_traces.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "engines" / "python" / "src"))

from biodp.align import needleman_wunsch  # noqa: E402
from biodp.scoring import GapPenalty, SubstitutionMatrix  # noqa: E402
from biodp.trace import dumps, to_dptrace  # noqa: E402

GOLDEN = ROOT / "fixtures" / "golden" / "needleman_wunsch"
MATRICES = ROOT / "fixtures" / "matrices"
OUT = ROOT / "web" / "public" / "traces"

# Longer inputs that show the visualiser at a realistic density. Each one is
# chosen to display something the small fixtures cannot.
DEMOS: list[dict[str, Any]] = [
    {
        "id": "demo-dna-indel",
        "title": "DNA with one 3 bp deletion",
        "description": "Two 30-ish bp sequences differing by substitutions and a single 3 bp deletion; "
        "affine gaps keep the deletion in one block.",
        "seq_a": "ATGGCGTACGTTAGCCTAGGATCCGATTACAGG",
        "seq_b": "ATGGCTTACGTTAGCGGATCCGATAACAGG",
        "substitution": {"type": "match_mismatch", "match": 2, "mismatch": -3, "alphabet": "ACGT"},
        "gap": {"open": -5, "extend": -2},
    },
    {
        "id": "demo-homopolymer-ties",
        "title": "Homopolymer: many co-optimal alignments",
        "description": "A run of A's against a shorter run. Every placement of the gap scores the same, "
        "so the optimal-path band fills a whole region of the table.",
        "seq_a": "GCAAAAAAAATC",
        "seq_b": "GCAAAAATC",
        "substitution": {"type": "match_mismatch", "match": 1, "mismatch": -1, "alphabet": "ACGT"},
        "gap": {"open": -2, "extend": -2},
    },
    {
        "id": "demo-protein-blosum62",
        "title": "Protein fragments, BLOSUM62",
        "description": "A protein fragment against a variant with two substitutions (I/L, F/W) and a "
        "3-residue deletion, scored with BLOSUM62 and BLAST-style affine gaps (open -11, extend -1).",
        "seq_a": "MKTAYIAKQRQISFVKSHFSRQ",
        "seq_b": "MKTAYLAKQRQISWVKSRQ",
        "substitution": {"type": "named", "name": "BLOSUM62"},
        "gap": {"open": -11, "extend": -1},
    },
]


def matrix_from(spec: dict[str, Any]) -> SubstitutionMatrix:
    if spec["type"] == "match_mismatch":
        return SubstitutionMatrix.match_mismatch(spec["match"], spec["mismatch"], spec["alphabet"])
    return SubstitutionMatrix.from_ncbi_file(MATRICES / f"{spec['name']}.txt")


def export(case_id: str, title: str, description: str, seq_a: str, seq_b: str,
           substitution: dict[str, Any], gap: dict[str, int]) -> dict[str, Any]:
    result = needleman_wunsch(
        seq_a, seq_b, matrix=matrix_from(substitution), gap=GapPenalty(gap["open"], gap["extend"])
    )
    (OUT / f"{case_id}.json").write_text(dumps(to_dptrace(result)) + "\n")
    return {
        "id": case_id,
        "title": title,
        "description": description,
        "file": f"{case_id}.json",
        "score": result.score,
        "optimal_count": str(result.optimal_count),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for path in sorted(GOLDEN.glob("*.json")):
        fixture = json.loads(path.read_text())
        inp = fixture["input"]
        manifest.append(export(
            fixture["id"], fixture["id"].split("-", 2)[-1].replace("-", " ").capitalize(),
            fixture["description"], inp["seq_a"], inp["seq_b"], inp["substitution"], inp["gap"],
        ))
    for demo in DEMOS:
        manifest.append(export(
            demo["id"], demo["title"], demo["description"], demo["seq_a"], demo["seq_b"],
            demo["substitution"], demo["gap"],
        ))
    (OUT / "index.json").write_text(json.dumps({"traces": manifest}, indent=2) + "\n")
    for entry in manifest:
        print(f"{entry['id']}: score={entry['score']} optimal={entry['optimal_count']}")


if __name__ == "__main__":
    main()
