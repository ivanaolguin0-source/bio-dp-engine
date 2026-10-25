"""Serialise a filled DP into the language-neutral DPTrace contract.

The schema lives at ``contracts/schemas/dptrace.schema.json``. Every engine,
Python or Java, emits this shape, and the frontend only ever replays it.

Two views of the table are exported:

``layers``
    The raw Gotoh states (M, X, Y), each with scores and predecessor-state
    bitmasks. Needed to explain affine gaps faithfully.
``projection``
    The single-matrix view a reader expects from a textbook figure:
    ``H = max(M, X, Y)`` per cell, plus a *move* bitmask saying which
    directions reach that maximum (1 diagonal, 2 up, 4 left). For linear
    gaps this is exactly the classic Needleman-Wunsch table and arrows.
"""

from __future__ import annotations

import json
from itertools import islice
from typing import Any, Final

import numpy as np
import numpy.typing as npt

import biodp
from biodp.align.model import NEG_INF, State
from biodp.align.needleman_wunsch import GlobalAlignmentResult

SCHEMA_VERSION: Final = "1.0.0"

# Move bits for the projection. Deliberately distinct names from the FROM_*
# predecessor-state flags even though the numeric values coincide.
MOVE_DIAGONAL: Final = 1
MOVE_UP: Final = 2
MOVE_LEFT: Final = 4
_STATE_MOVE: Final = {State.M: MOVE_DIAGONAL, State.X: MOVE_UP, State.Y: MOVE_LEFT}

# About 300 x 300. Larger tables are slow to transfer as JSON and too dense
# to animate cell by cell; request compute-only output for those.
DEFAULT_MAX_TRACE_CELLS: Final = 90_000


def _grid(values: npt.NDArray[np.int64]) -> list[list[int | None]]:
    """Score grid with unreachable cells as ``null`` (JSON has no infinity)."""
    return [[None if v == NEG_INF else int(v) for v in row] for row in values.tolist()]


def to_dptrace(
    result: GlobalAlignmentResult,
    *,
    max_alignments: int = 16,
    max_cells: int = DEFAULT_MAX_TRACE_CELLS,
    include_layers: bool = True,
) -> dict[str, Any]:
    """Build a DPTrace document for one Needleman-Wunsch run.

    Args:
        result: Output of :func:`biodp.align.needleman_wunsch`.
        max_alignments: Cap on alignments listed; ``optimal_count`` is always
            exact and ``alignments_truncated`` says whether the list is partial.
        max_cells: Refuse tables larger than this (visualisation mode only).
        include_layers: Omit the raw M/X/Y layers to roughly quarter the size.

    The optimal count is serialised as a decimal string because JavaScript
    numbers lose integer precision above 2**53.
    """
    dp = result.dp
    rows, cols = dp.n + 1, dp.m + 1
    if rows * cols > max_cells:
        raise ValueError(
            f"table has {rows * cols:,} cells; visualisation traces are capped at {max_cells:,}. "
            "Request the score and alignment only for inputs this large."
        )

    projected = dp.scores.max(axis=0)
    moves = np.zeros((rows, cols), dtype=np.uint8)
    for state, bit in _STATE_MOVE.items():
        # A state contributes an arrow only if it attains the cell maximum AND
        # was itself reached from somewhere. The second condition keeps the
        # origin (score 0, no predecessor) and unreachable cells arrow-free.
        attains = (dp.scores[state] == projected) & (dp.pointers[state] != 0)
        moves |= attains.astype(np.uint8) * bit

    listed = list(islice(result.alignments, max_alignments))
    count = result.optimal_count

    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "algorithm": "needleman-wunsch",
        "engine": {"language": "python", "name": "biodp", "version": biodp.__version__},
        "inputs": {
            "seq_a": dp.seq_a,
            "seq_b": dp.seq_b,
            "scoring": {
                "substitution": {
                    "name": dp.matrix.name,
                    "alphabet": dp.matrix.alphabet,
                    "scores": dp.matrix.scores.tolist(),
                },
                "gap": {"open": dp.gap.open_score, "extend": dp.gap.extend_score},
            },
        },
        "shape": [rows, cols],
        "fill_order": {"kind": "row-major"},
        "projection": {
            "scores": _grid(projected),
            "moves": moves.tolist(),
            "on_optimal_path": result.alignments.on_path_cells().tolist(),
        },
        "result": {
            "score": dp.score,
            "optimal_count": str(count),
            "alignments": [
                {
                    "aligned_a": a.aligned_a,
                    "aligned_b": a.aligned_b,
                    "cigar": a.cigar(),
                    "path": [list(cell) for cell in a.path],
                }
                for a in listed
            ],
            "alignments_truncated": count > len(listed),
        },
    }
    if include_layers:
        document["layers"] = [
            {
                "name": state.name,
                "scores": _grid(dp.scores[state]),
                "predecessors": dp.pointers[state].tolist(),
            }
            for state in State
        ]
    return document


def dumps(document: dict[str, Any]) -> str:
    """Compact JSON encoding of a DPTrace document."""
    return json.dumps(document, separators=(",", ":"))
