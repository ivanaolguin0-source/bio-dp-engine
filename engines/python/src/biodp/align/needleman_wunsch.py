"""Needleman-Wunsch global alignment with affine gaps (Gotoh's recurrences).

This module only *fills* the tables. Reconstructing alignments from them
lives in :mod:`biodp.align.traceback`, so the two halves can be tested,
optimised, and replaced (for example by Hirschberg's linear-space variant)
independently.

Two interchangeable fill kernels are provided:

``"reference"``
    Plain nested loops that mirror the recurrence line by line. Slow, but
    it is the version to read alongside ``docs/algorithms/needleman-wunsch.md``.
``"vectorized"`` (default)
    Fills one row at a time with NumPy. The diagonal (M) and vertical (X)
    states depend only on the previous row, so they vectorise directly. The
    horizontal state (Y) depends on its left neighbour in the same row; that
    dependency is resolved exactly with a running maximum (see
    ``_fill_y_row``). The test suite asserts both kernels produce
    bit-identical score and pointer tensors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal

import numpy as np
import numpy.typing as npt

from biodp.align.model import (
    FROM_M,
    FROM_X,
    FROM_Y,
    NEG_INF,
    UNREACHABLE_BELOW,
    DPMatrices,
    PointerTensor,
    ScoreTensor,
    State,
)
from biodp.align.traceback import Alignment, OptimalAlignments
from biodp.scoring import GapPenalty, SubstitutionMatrix
from biodp.scoring.substitution import IndexArray, IntMatrix

Kernel = Literal["vectorized", "reference"]

# Three int64 score layers plus three uint8 pointer layers is 27 bytes per
# cell, so the default ceiling is roughly 675 MB. Above that, full-matrix
# traceback is the wrong tool; use a linear-space method instead.
DEFAULT_MAX_CELLS: Final = 25_000_000

_M, _X, _Y = int(State.M), int(State.X), int(State.Y)
Row = npt.NDArray[np.int64]


@dataclass(frozen=True, eq=False)
class GlobalAlignmentResult:
    """Everything produced by one Needleman-Wunsch run."""

    dp: DPMatrices
    alignments: OptimalAlignments

    @property
    def score(self) -> int:
        return self.dp.score

    @property
    def optimal_count(self) -> int:
        """Exact number of distinct co-optimal alignments (may be astronomically large)."""
        return self.alignments.count

    @property
    def best(self) -> Alignment:
        """The canonical optimal alignment under the documented tie-break."""
        return self.alignments.best()


def needleman_wunsch(
    seq_a: str,
    seq_b: str,
    *,
    matrix: SubstitutionMatrix,
    gap: GapPenalty,
    kernel: Kernel = "vectorized",
    max_cells: int = DEFAULT_MAX_CELLS,
) -> GlobalAlignmentResult:
    """Globally align two sequences and expose every optimal alignment.

    Scoring has no defaults on purpose: alignment results are meaningless
    without knowing the matrix and gap model, so callers must choose them.

    Args:
        seq_a: First sequence (rows of the DP table). Case-insensitive.
        seq_b: Second sequence (columns). Case-insensitive.
        matrix: Substitution scores; every residue must be in its alphabet.
        gap: Affine (or linear) gap scores.
        kernel: Fill implementation; results are identical either way.
        max_cells: Refuse problems whose table exceeds this many cells.

    Returns:
        A result holding the filled tables and a lazy view over all
        co-optimal alignments.

    Raises:
        InvalidResidueError: a sequence contains a symbol outside the alphabet.
        ValueError: the table would exceed ``max_cells``.
    """
    dp = fill_matrices(seq_a, seq_b, matrix=matrix, gap=gap, kernel=kernel, max_cells=max_cells)
    return GlobalAlignmentResult(dp=dp, alignments=OptimalAlignments(dp))


def fill_matrices(
    seq_a: str,
    seq_b: str,
    *,
    matrix: SubstitutionMatrix,
    gap: GapPenalty,
    kernel: Kernel = "vectorized",
    max_cells: int = DEFAULT_MAX_CELLS,
) -> DPMatrices:
    """Run only the forward pass and return the completed tables."""
    a_idx = matrix.encode(seq_a, "seq_a")
    b_idx = matrix.encode(seq_b, "seq_b")
    cells = (len(a_idx) + 1) * (len(b_idx) + 1)
    if cells > max_cells:
        raise ValueError(
            f"DP table would have {cells:,} cells (limit {max_cells:,}). "
            "Full-matrix traceback needs memory proportional to len(a) * len(b); "
            "use a linear-space method such as Hirschberg's algorithm for inputs this large."
        )

    if kernel == "vectorized":
        scores, pointers = _fill_vectorized(a_idx, b_idx, matrix.scores, gap)
    elif kernel == "reference":
        scores, pointers = _fill_reference(a_idx, b_idx, matrix.scores, gap)
    else:
        raise ValueError(f"unknown kernel {kernel!r}; expected 'vectorized' or 'reference'")

    scores.setflags(write=False)
    pointers.setflags(write=False)
    return DPMatrices(seq_a.upper(), seq_b.upper(), matrix, gap, scores, pointers)


# ---------------------------------------------------------------------- reference kernel


def _store(
    scores: ScoreTensor, pointers: PointerTensor, cell: tuple[int, int, int], value: int, mask: int
) -> None:
    """Write one cell, normalising anything derived from -inf to exactly NEG_INF."""
    if value < UNREACHABLE_BELOW:
        return  # tensors were initialised to (NEG_INF, 0)
    scores[cell] = value
    pointers[cell] = mask


def _argmax3(from_m: int, from_x: int, from_y: int) -> tuple[int, int]:
    """Best of three candidates plus a bitmask of every candidate that ties it."""
    best = max(from_m, from_x, from_y)
    mask = (FROM_M if from_m == best else 0) | (FROM_X if from_x == best else 0) | (FROM_Y if from_y == best else 0)
    return best, mask


def _fill_reference(
    a_idx: IndexArray, b_idx: IndexArray, sub: IntMatrix, gap: GapPenalty
) -> tuple[ScoreTensor, PointerTensor]:
    n, m = len(a_idx), len(b_idx)
    g_open, g_ext = gap.open_score, gap.extend_score
    S = np.full((3, n + 1, m + 1), NEG_INF, dtype=np.int64)
    P = np.zeros((3, n + 1, m + 1), dtype=np.uint8)

    # Base cases. The empty-vs-empty alignment scores 0 and is the start node.
    S[_M, 0, 0] = 0
    # A prefix of seq_a against nothing is one leading gap in seq_b...
    for i in range(1, n + 1):
        S[_X, i, 0] = g_open + (i - 1) * g_ext
        P[_X, i, 0] = FROM_M if i == 1 else FROM_X
    # ...and a prefix of seq_b against nothing is one leading gap in seq_a.
    for j in range(1, m + 1):
        S[_Y, 0, j] = g_open + (j - 1) * g_ext
        P[_Y, 0, j] = FROM_M if j == 1 else FROM_Y

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = int(sub[a_idx[i - 1], b_idx[j - 1]])
            # M: a_i opposite b_j, arriving from any state at (i-1, j-1).
            best, mask = _argmax3(int(S[_M, i - 1, j - 1]), int(S[_X, i - 1, j - 1]), int(S[_Y, i - 1, j - 1]))
            _store(S, P, (_M, i, j), best + s, mask)
            # X: a_i opposite a gap. Extending an X run costs g_ext; anything else opens.
            best, mask = _argmax3(int(S[_M, i - 1, j]) + g_open, int(S[_X, i - 1, j]) + g_ext, int(S[_Y, i - 1, j]) + g_open)
            _store(S, P, (_X, i, j), best, mask)
            # Y: b_j opposite a gap. Mirror image of X.
            best, mask = _argmax3(int(S[_M, i, j - 1]) + g_open, int(S[_X, i, j - 1]) + g_open, int(S[_Y, i, j - 1]) + g_ext)
            _store(S, P, (_Y, i, j), best, mask)
    return S, P


# ---------------------------------------------------------------------- vectorized kernel


def _tie_mask(from_m: Row, from_x: Row, from_y: Row, best: Row) -> PointerTensor:
    """Vectorised ``_argmax3`` bitmask: which candidates equal the row maximum."""
    mask = (
        (from_m == best).astype(np.uint8) * FROM_M
        | (from_x == best).astype(np.uint8) * FROM_X
        | (from_y == best).astype(np.uint8) * FROM_Y
    )
    return np.asarray(mask, dtype=np.uint8)


def _clamp_row(S: ScoreTensor, P: PointerTensor, state: int, i: int) -> None:
    dead = S[state, i] < UNREACHABLE_BELOW
    S[state, i][dead] = NEG_INF
    P[state, i][dead] = 0


def _fill_y_row(S: ScoreTensor, P: PointerTensor, i: int, g_open: int, g_ext: int, k_ext: Row) -> None:
    """Resolve the in-row dependency of the horizontal-gap state exactly.

    The recurrence ``Y[j] = max(M[j-1] + o, X[j-1] + o, Y[j-1] + e)`` unrolls to

        Y[j] = o + (j-1)*e + max_{k <= j-1} (B[k] - k*e),   B[k] = max(M[k], X[k])

    because a horizontal gap ending at column j must have opened right after
    some column k. The inner maximum is a running (prefix) maximum, which
    NumPy computes in one vectorised pass. Integers keep this exact.
    """
    best_open = np.maximum(S[_M, i], S[_X, i])
    running = np.maximum.accumulate(best_open - k_ext)
    y = np.full_like(best_open, NEG_INF)
    y[1:] = g_open + k_ext[:-1] + running[:-1]
    S[_Y, i] = y
    P[_Y, i, 1:] = _tie_mask(S[_M, i, :-1] + g_open, S[_X, i, :-1] + g_open, y[:-1] + g_ext, y[1:])
    P[_Y, i, 0] = 0
    _clamp_row(S, P, _Y, i)


def _fill_vectorized(
    a_idx: IndexArray, b_idx: IndexArray, sub: IntMatrix, gap: GapPenalty
) -> tuple[ScoreTensor, PointerTensor]:
    n, m = len(a_idx), len(b_idx)
    g_open, g_ext = gap.open_score, gap.extend_score
    S = np.full((3, n + 1, m + 1), NEG_INF, dtype=np.int64)
    P = np.zeros((3, n + 1, m + 1), dtype=np.uint8)
    k_ext = np.arange(m + 1, dtype=np.int64) * g_ext  # k * extend, for the Y scan

    S[_M, 0, 0] = 0
    _fill_y_row(S, P, 0, g_open, g_ext, k_ext)  # row 0 is one leading gap in seq_a

    for i in range(1, n + 1):
        up_m, up_x, up_y = S[_M, i - 1], S[_X, i - 1], S[_Y, i - 1]

        # M row: every candidate comes from row i-1, so it is a pure array op.
        dm, dx, dy = up_m[:-1], up_x[:-1], up_y[:-1]
        best = np.maximum(np.maximum(dm, dx), dy)
        P[_M, i, 1:] = _tie_mask(dm, dx, dy, best)
        S[_M, i, 1:] = best + sub[a_idx[i - 1], b_idx]
        _clamp_row(S, P, _M, i)

        # X row: also depends only on row i-1. Column 0 falls out of the same
        # formula because M[0,0] = 0 and every other boundary cell is -inf.
        vm, vx, vy = up_m + g_open, up_x + g_ext, up_y + g_open
        best = np.maximum(np.maximum(vm, vx), vy)
        P[_X, i] = _tie_mask(vm, vx, vy, best)
        S[_X, i] = best
        _clamp_row(S, P, _X, i)

        _fill_y_row(S, P, i, g_open, g_ext, k_ext)
    return S, P
