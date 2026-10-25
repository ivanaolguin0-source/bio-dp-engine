"""Hirschberg's linear-space global alignment (linear gap scores).

Full Needleman-Wunsch keeps every pointer so it can trace back, which costs
memory proportional to ``len(a) * len(b)``. Hirschberg (1975) observed that
the *score* of a prefix alignment needs only the previous row, and used that
to find, without any pointers, the column where an optimal alignment crosses
the middle row of the table. Splitting there yields two independent
subproblems of half the height; recursing on both and concatenating the
results gives an optimal alignment of the whole.

Cost: time stays O(nm) (about twice the full fill, since each level of the
recursion re-scores half the area of the level above), while working memory
falls to O(n + m).

Scope, stated plainly:

* **Linear gaps only.** With affine gaps an optimal alignment can cross the
  middle row inside a gap, and splitting there would charge that gap's
  opening score twice. Myers and Miller (1988) extend the method to affine
  gaps by tracking the gap state at the split; that is not implemented here,
  so affine input is rejected rather than answered approximately.
* **One optimal alignment, not all of them.** Enumerating or counting
  co-optimal alignments needs the pointer graph that this method exists to
  avoid storing. Use :func:`biodp.align.needleman_wunsch` for that.
"""

from __future__ import annotations

from typing import Final

import numpy as np
import numpy.typing as npt

from biodp.align.needleman_wunsch import needleman_wunsch
from biodp.align.traceback import Alignment
from biodp.scoring import GAP_SYMBOL, GapPenalty, SubstitutionMatrix
from biodp.scoring.substitution import IndexArray, IntMatrix

Row = npt.NDArray[np.int64]

# Subproblems at or below this many DP cells are solved with the full-table
# Needleman-Wunsch engine instead of being split further. The bound is a
# constant, so total working memory stays O(n + m); it only removes Python
# call overhead from the tiny leaves of the recursion. Pass ``base_cells=0``
# for textbook Hirschberg, which splits all the way down to single residues.
DEFAULT_BASE_CELLS: Final = 4096


def _require_linear(gap: GapPenalty) -> None:
    if not gap.is_linear:
        raise ValueError(
            f"Hirschberg's split is exact only for linear gaps (open == extend); got open={gap.open_score}, "
            f"extend={gap.extend_score}. Affine gaps need the Myers-Miller extension, which is not implemented. "
            "Use needleman_wunsch() for affine scoring."
        )


def last_row(a_idx: IndexArray, b_idx: IndexArray, sub: IntMatrix, g: int) -> Row:
    """Global alignment scores of all of ``a`` against every prefix of ``b``.

    Returns ``row`` where ``row[j]`` is the optimal score of ``a`` versus
    ``b[:j]``. Only one row of the table is ever held.

    Each new row is computed in two vectorised steps. First the candidates
    that come from the row above (diagonal and vertical moves). Then the
    horizontal dependency within the row, ``H[j] = max(C[j], H[j-1] + g)``,
    which unrolls to ``H[j] = j*g + max_{k <= j} (C[k] - k*g)``, a running
    maximum. Integer scores keep this exact.
    """
    m = len(b_idx)
    k_g = np.arange(m + 1, dtype=np.int64) * g
    row = k_g.copy()  # a = empty prefix: j leading gaps
    candidates = np.empty(m + 1, dtype=np.int64)
    for ai in a_idx:
        candidates[0] = row[0] + g
        candidates[1:] = np.maximum(row[:-1] + sub[ai, b_idx], row[1:] + g)
        row = k_g + np.maximum.accumulate(candidates - k_g)
    return np.asarray(row, dtype=np.int64)


def global_score(seq_a: str, seq_b: str, *, matrix: SubstitutionMatrix, gap: GapPenalty) -> int:
    """Optimal global alignment score in O(min(n, m)) memory, without the alignment.

    The row is laid along the shorter sequence. Swapping the sequences is
    exact provided the substitution matrix is transposed with them, because
    score(a, b, S) equals score(b, a, S transposed).
    """
    _require_linear(gap)
    a_idx = matrix.encode(seq_a, "seq_a")
    b_idx = matrix.encode(seq_b, "seq_b")
    sub = matrix.scores
    if len(b_idx) > len(a_idx):
        a_idx, b_idx, sub = b_idx, a_idx, sub.T
    return int(last_row(a_idx, b_idx, sub, gap.open_score)[-1])


def hirschberg(
    seq_a: str,
    seq_b: str,
    *,
    matrix: SubstitutionMatrix,
    gap: GapPenalty,
    base_cells: int = DEFAULT_BASE_CELLS,
) -> Alignment:
    """One optimal global alignment, computed in linear space.

    Args:
        seq_a: First sequence. Case-insensitive.
        seq_b: Second sequence. Case-insensitive.
        matrix: Substitution scores.
        gap: Must be linear (``open_score == extend_score``).
        base_cells: Subproblems this small are solved directly; see
            ``DEFAULT_BASE_CELLS``.

    Returns:
        An optimal alignment. Its score always equals the Needleman-Wunsch
        optimum; when several alignments tie, it may be a different one from
        Needleman-Wunsch's canonical choice.

    Raises:
        ValueError: affine gap scores were given.
        InvalidResidueError: a sequence contains a symbol outside the alphabet.
    """
    _require_linear(gap)
    if base_cells < 0:
        raise ValueError("base_cells cannot be negative")
    a_idx = matrix.encode(seq_a, "seq_a")
    b_idx = matrix.encode(seq_b, "seq_b")
    a, b = seq_a.upper(), seq_b.upper()
    g = gap.open_score
    sub = matrix.scores
    top: list[str] = []
    bottom: list[str] = []

    def solve(a_lo: int, a_hi: int, b_lo: int, b_hi: int) -> None:
        n, m = a_hi - a_lo, b_hi - b_lo
        if n == 0:
            top.append(GAP_SYMBOL * m)
            bottom.append(b[b_lo:b_hi])
            return
        if m == 0:
            top.append(a[a_lo:a_hi])
            bottom.append(GAP_SYMBOL * n)
            return
        if n == 1 or m == 1 or (n + 1) * (m + 1) <= base_cells:
            best = needleman_wunsch(a[a_lo:a_hi], b[b_lo:b_hi], matrix=matrix, gap=gap).best
            top.append(best.aligned_a)
            bottom.append(best.aligned_b)
            return
        # Divide: score the top half forwards and the bottom half backwards.
        mid = a_lo + n // 2
        forward = last_row(a_idx[a_lo:mid], b_idx[b_lo:b_hi], sub, g)
        backward = last_row(a_idx[mid:a_hi][::-1], b_idx[b_lo:b_hi][::-1], sub, g)[::-1]
        # forward[k] + backward[k] is the best score of any alignment that
        # crosses the middle row at column b_lo + k. The maximum is the
        # optimum, so splitting there loses nothing.
        k = int(np.argmax(forward + backward))
        solve(a_lo, mid, b_lo, b_lo + k)
        solve(mid, a_hi, b_lo + k, b_hi)

    # Recursion depth is about log2(len(a)), so even very long inputs stay
    # far below Python's recursion limit.
    solve(0, len(a), 0, len(b))
    aligned_a, aligned_b = "".join(top), "".join(bottom)
    return Alignment(aligned_a, aligned_b, _score_columns(aligned_a, aligned_b, matrix, g), _path(aligned_a, aligned_b))


def _score_columns(aligned_a: str, aligned_b: str, matrix: SubstitutionMatrix, g: int) -> int:
    """Linear-gap score of an alignment, computed column by column (vectorised)."""
    if not aligned_a:
        return 0
    top = np.frombuffer(aligned_a.encode("ascii"), dtype=np.uint8)
    bottom = np.frombuffer(aligned_b.encode("ascii"), dtype=np.uint8)
    gap_code = ord(GAP_SYMBOL)
    paired = (top != gap_code) & (bottom != gap_code)
    gaps = int(np.count_nonzero(~paired))
    pair_a = "".join(ch for ch, keep in zip(aligned_a, paired) if keep)
    pair_b = "".join(ch for ch, keep in zip(aligned_b, paired) if keep)
    pair_scores = matrix.scores[matrix.encode(pair_a), matrix.encode(pair_b)]
    return int(pair_scores.sum()) + gaps * g


def _path(aligned_a: str, aligned_b: str) -> tuple[tuple[int, int], ...]:
    i = j = 0
    cells = [(0, 0)]
    for x, y in zip(aligned_a, aligned_b):
        i += x != GAP_SYMBOL
        j += y != GAP_SYMBOL
        cells.append((i, j))
    return tuple(cells)
