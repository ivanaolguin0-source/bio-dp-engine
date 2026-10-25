"""The filled dynamic-programming tables shared by fill and traceback.

Alignment uses Gotoh's three-state formulation. Every cell ``(i, j)`` holds
three scores, one per *state*, where the state records which kind of column
ends the best partial alignment of ``seq_a[:i]`` against ``seq_b[:j]``:

    M  residue a_i opposite residue b_j     (diagonal move)
    X  residue a_i opposite a gap           (vertical move, consumes seq_a)
    Y  residue b_j opposite a gap           (horizontal move, consumes seq_b)

Each alignment corresponds to exactly one path through these states, so
counting state paths counts distinct alignments, with no double counting.
With linear gaps the three states collapse to the textbook single matrix
``H = max(M, X, Y)``; ``biodp.trace`` exports that projection for display.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Final

import numpy as np
import numpy.typing as npt

from biodp.scoring import GAP_SYMBOL, GapPenalty, SubstitutionMatrix

ScoreTensor = npt.NDArray[np.int64]
PointerTensor = npt.NDArray[np.uint8]


class State(IntEnum):
    M = 0
    X = 1
    Y = 2


# Pointer bitmask: which predecessor states achieved the cell's maximum.
# Several bits set means a tie, i.e. more than one optimal way to arrive.
FROM_M: Final = 1
FROM_X: Final = 2
FROM_Y: Final = 4
STATE_FLAGS: Final[tuple[tuple[State, int], ...]] = (
    (State.M, FROM_M),
    (State.X, FROM_X),
    (State.Y, FROM_Y),
)

# How far back the predecessor cell sits for each state's move.
CELL_STEP: Final[dict[State, tuple[int, int]]] = {
    State.M: (1, 1),
    State.X: (1, 0),
    State.Y: (0, 1),
}

# Minus infinity for int64 tables. A quarter of the int64 range leaves ample
# headroom so adding penalties to it never wraps around.
NEG_INF: Final = int(np.iinfo(np.int64).min // 4)
# Anything below this is "unreachable" no matter what was added to NEG_INF.
UNREACHABLE_BELOW: Final = NEG_INF // 2

# A node in the state graph: (state, i, j). Plain ints keep hashing cheap.
Node = tuple[int, int, int]
START: Final[Node] = (int(State.M), 0, 0)


@dataclass(frozen=True, eq=False)
class DPMatrices:
    """Completed score and pointer tensors for one global alignment problem.

    Attributes:
        seq_a: First sequence, uppercased. Indexes rows (``i``).
        seq_b: Second sequence, uppercased. Indexes columns (``j``).
        matrix: Substitution scores used for the fill.
        gap: Gap model used for the fill.
        scores: ``scores[state, i, j]``; unreachable cells hold ``NEG_INF``.
        pointers: ``pointers[state, i, j]`` bitmask of ``FROM_*`` flags.
            Unreachable cells and the start cell hold 0.
    """

    seq_a: str
    seq_b: str
    matrix: SubstitutionMatrix
    gap: GapPenalty
    scores: ScoreTensor
    pointers: PointerTensor

    @property
    def n(self) -> int:
        return len(self.seq_a)

    @property
    def m(self) -> int:
        return len(self.seq_b)

    @property
    def score(self) -> int:
        """Optimal global alignment score (best of the three end states)."""
        return int(np.max(self.scores[:, self.n, self.m]))

    @property
    def end_states(self) -> tuple[State, ...]:
        """Every state at ``(n, m)`` that attains the optimum, in priority order."""
        best = self.score
        return tuple(s for s in State if int(self.scores[s, self.n, self.m]) == best)

    def predecessors(self, node: Node) -> list[Node]:
        """Optimal predecessor nodes, ordered M, X, Y (the tie-break priority)."""
        state, i, j = node
        mask = int(self.pointers[state, i, j])
        di, dj = CELL_STEP[State(state)]
        return [(int(s), i - di, j - dj) for s, flag in STATE_FLAGS if mask & flag]

    def column(self, node: Node) -> tuple[str, str]:
        """The alignment column contributed by arriving at ``node``."""
        state, i, j = node
        if state == State.M:
            return self.seq_a[i - 1], self.seq_b[j - 1]
        if state == State.X:
            return self.seq_a[i - 1], GAP_SYMBOL
        return GAP_SYMBOL, self.seq_b[j - 1]
