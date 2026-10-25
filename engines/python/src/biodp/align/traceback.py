"""Reconstructing alignments from filled DP tables.

The pointer bitmasks define a directed acyclic graph whose nodes are
``(state, i, j)`` and whose edges point from a node to each predecessor that
tied for its maximum. Every path from an optimal end node back to the start
node ``(M, 0, 0)`` spells one optimal alignment, and different paths spell
different alignments. Three consequences follow, and this module provides one
operation for each:

* **Counting** the optimal alignments is counting paths in a DAG, which is
  itself a small dynamic program (``OptimalAlignments.count``). The count is
  exact and uses Python's arbitrary-precision integers, because it grows
  combinatorially: ``A*200`` against ``A*100`` has C(200, 100), about 9e58.
* **Enumerating** them is a depth-first walk that yields alignments lazily
  (``iter(OptimalAlignments)``), so asking for the first ten never touches
  the rest.
* **Sampling** one uniformly at random uses the path counts as weights
  (``OptimalAlignments.sample``). This is the unbiased way to pick a
  representative when the canonical tie-break would systematically push
  gaps to one side.

Canonical tie-break: walking backwards, prefer M, then X, then Y. Preferring
the diagonal near the end defers gaps toward the start of the alignment, so
the canonical alignment places gaps as far left (5') as the ties allow.
"""

from __future__ import annotations

import random
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np
import numpy.typing as npt

from biodp.align.model import START, DPMatrices, Node, State
from biodp.scoring import GAP_SYMBOL

# A persistent singly linked list of nodes, head = earliest node. Sharing tails
# lets the depth-first search branch without copying partial paths.
_Chain = tuple[Node, "_Chain"] | None

_MIDLINE_MATCH: Final = "|"
_MIDLINE_MISMATCH: Final = "."
_MIDLINE_GAP: Final = " "


@dataclass(frozen=True)
class Alignment:
    """One pairwise alignment.

    Attributes:
        aligned_a: ``seq_a`` with ``-`` inserted at gap positions.
        aligned_b: ``seq_b`` with ``-`` inserted at gap positions.
        score: Total alignment score under the scoring used to produce it.
        path: DP cells visited, from ``(0, 0)`` to ``(len(a), len(b))``.
            This is what the visualiser animates.
    """

    aligned_a: str
    aligned_b: str
    score: int
    path: tuple[tuple[int, int], ...]

    def __post_init__(self) -> None:
        if len(self.aligned_a) != len(self.aligned_b):
            raise ValueError("aligned strings must have equal length")
        if len(self.path) != len(self.aligned_a) + 1:
            raise ValueError("path must have exactly one more cell than the alignment has columns")

    @property
    def length(self) -> int:
        """Number of alignment columns."""
        return len(self.aligned_a)

    @property
    def matches(self) -> int:
        return sum(x == y != GAP_SYMBOL for x, y in zip(self.aligned_a, self.aligned_b))

    @property
    def mismatches(self) -> int:
        return sum(
            x != y and GAP_SYMBOL not in (x, y) for x, y in zip(self.aligned_a, self.aligned_b)
        )

    @property
    def gap_columns(self) -> int:
        return self.aligned_a.count(GAP_SYMBOL) + self.aligned_b.count(GAP_SYMBOL)

    @property
    def identity(self) -> float:
        """Identical columns divided by alignment length.

        Percent identity has several competing definitions (the denominator
        may be alignment length, the shorter sequence, or aligned columns
        only). This uses alignment length. Returns NaN for an empty alignment.
        """
        return self.matches / self.length if self.length else float("nan")

    def midline(self) -> str:
        """``|`` identical, ``.`` substitution, space for a gap column."""
        out = []
        for x, y in zip(self.aligned_a, self.aligned_b):
            if GAP_SYMBOL in (x, y):
                out.append(_MIDLINE_GAP)
            else:
                out.append(_MIDLINE_MATCH if x == y else _MIDLINE_MISMATCH)
        return "".join(out)

    def cigar(self, extended: bool = True) -> str:
        """CIGAR string treating ``seq_a`` as the reference and ``seq_b`` as the query.

        ``D`` marks a reference residue opposite a gap, ``I`` a query residue
        opposite a gap. With ``extended=True`` aligned pairs are split into
        ``=`` (identical) and ``X`` (substitution), as in SAM v1.4+; otherwise
        both are ``M``.
        """
        ops = []
        for x, y in zip(self.aligned_a, self.aligned_b):
            if y == GAP_SYMBOL:
                ops.append("D")
            elif x == GAP_SYMBOL:
                ops.append("I")
            elif not extended:
                ops.append("M")
            else:
                ops.append("=" if x == y else "X")
        runs: list[str] = []
        for op in ops:
            if runs and runs[-1][-1] == op:
                runs[-1] = f"{int(runs[-1][:-1]) + 1}{op}"
            else:
                runs.append(f"1{op}")
        return "".join(runs)

    def format(self, width: int = 60) -> str:
        """Three-line block layout, wrapped every ``width`` columns."""
        if width < 1:
            raise ValueError("width must be positive")
        mid = self.midline()
        blocks = [
            "\n".join((self.aligned_a[k : k + width], mid[k : k + width], self.aligned_b[k : k + width]))
            for k in range(0, max(self.length, 1), width)
        ]
        return "\n\n".join(blocks)

    def __str__(self) -> str:
        return self.format()


class OptimalAlignments:
    """A lazy, exact view over every co-optimal alignment in a filled DP.

    Construction does one pass over the *optimal subgraph* (nodes that lie on
    at least one optimal path). For typical inputs that subgraph is a thin
    band around the main path, so this is far cheaper than the fill. For
    inputs dominated by ties (homopolymers) it can approach the full table.
    """

    def __init__(self, dp: DPMatrices) -> None:
        self._dp = dp
        self._ends: tuple[Node, ...] = tuple((int(s), dp.n, dp.m) for s in dp.end_states)
        self._prefix_counts = _count_optimal_prefixes(dp, self._ends)

    @property
    def dp(self) -> DPMatrices:
        return self._dp

    @property
    def score(self) -> int:
        return self._dp.score

    @property
    def count(self) -> int:
        """Exact number of distinct optimal alignments.

        Exposed as a property rather than ``__len__`` because ``len()`` must
        fit in a machine integer and this number frequently does not.
        """
        return sum(self._prefix_counts[end] for end in self._ends)

    def __iter__(self) -> Iterator[Alignment]:
        """Yield every optimal alignment, canonical one first, each exactly once."""
        dp = self._dp
        stack: list[tuple[Node, _Chain]] = [(end, None) for end in reversed(self._ends)]
        while stack:
            node, chain = stack.pop()
            if node == START:
                yield self._build(chain)
                continue
            extended: _Chain = (node, chain)
            for pred in reversed(dp.predecessors(node)):
                stack.append((pred, extended))

    def best(self) -> Alignment:
        """The canonical optimal alignment (see module docstring for the tie-break)."""
        return next(iter(self))

    def take(self, limit: int) -> list[Alignment]:
        """Up to ``limit`` alignments in canonical enumeration order."""
        if limit < 0:
            raise ValueError("limit cannot be negative")
        out: list[Alignment] = []
        for alignment in self:
            if len(out) == limit:
                break
            out.append(alignment)
        return out

    def sample(self, rng: random.Random | None = None) -> Alignment:
        """Draw one optimal alignment uniformly at random.

        At each step the walk picks a predecessor with probability
        proportional to the number of optimal prefixes ending there, so the
        product of the step probabilities telescopes to ``1 / count`` for
        every alignment.
        """
        rng = rng or random.Random()
        node = self._weighted_choice(rng, self._ends)
        chain: _Chain = None
        while node != START:
            chain = (node, chain)
            node = self._weighted_choice(rng, self._dp.predecessors(node))
        return self._build(chain)

    def on_path_cells(self) -> npt.NDArray[np.bool_]:
        """Boolean ``(n+1, m+1)`` mask of cells used by at least one optimal alignment."""
        mask = np.zeros((self._dp.n + 1, self._dp.m + 1), dtype=np.bool_)
        for _, i, j in self._prefix_counts:
            mask[i, j] = True
        return mask

    # ------------------------------------------------------------------ internals

    def _weighted_choice(self, rng: random.Random, nodes: Sequence[Node]) -> Node:
        weights = [self._prefix_counts[node] for node in nodes]
        ticket = rng.randrange(sum(weights))
        for node, weight in zip(nodes, weights):
            if ticket < weight:
                return node
            ticket -= weight
        raise AssertionError("unreachable: ticket exceeded total weight")

    def _build(self, chain: _Chain) -> Alignment:
        dp = self._dp
        top: list[str] = []
        bottom: list[str] = []
        path: list[tuple[int, int]] = [(0, 0)]
        while chain is not None:
            node, chain = chain
            x, y = dp.column(node)
            top.append(x)
            bottom.append(y)
            path.append((node[1], node[2]))
        return Alignment("".join(top), "".join(bottom), dp.score, tuple(path))


def _count_optimal_prefixes(dp: DPMatrices, ends: Sequence[Node]) -> dict[Node, int]:
    """Count optimal prefix paths for every node on the optimal subgraph.

    Pass 1 walks backwards from the end nodes, bucketing reachable nodes by
    anti-diagonal ``d = i + j``. Every predecessor lies on a strictly smaller
    anti-diagonal, so processing buckets in decreasing ``d`` visits each node
    after all of its successors. Pass 2 runs forward over the same buckets:
    ``count(START) = 1`` and ``count(v) = sum(count(p) for p in preds(v))``.

    Iterative by construction; recursion would overflow Python's stack on
    sequences only a few thousand residues long.
    """
    diagonals: list[set[Node]] = [set() for _ in range(dp.n + dp.m + 1)]
    diagonals[dp.n + dp.m].update(ends)
    for d in range(dp.n + dp.m, 0, -1):
        for node in diagonals[d]:
            for pred in dp.predecessors(node):
                diagonals[pred[1] + pred[2]].add(pred)

    counts: dict[Node, int] = {}
    for bucket in diagonals:
        for node in bucket:
            if node == START:
                counts[node] = 1
            else:
                counts[node] = sum(counts[pred] for pred in dp.predecessors(node))
    if START not in counts:
        raise AssertionError("optimal subgraph does not reach the start node; pointer tables are inconsistent")
    return counts


def rescore(aligned_a: str, aligned_b: str, dp: DPMatrices) -> int:
    """Independently score an alignment column by column under ``dp``'s scoring.

    Used by the test suite to confirm that every enumerated alignment really
    achieves the reported optimum. A change of gap *kind* (a gap in ``a``
    directly after a gap in ``b``) opens a new gap, matching the recurrences.
    """
    total = 0
    previous: State | None = None
    for x, y in zip(aligned_a, aligned_b):
        if x == GAP_SYMBOL and y == GAP_SYMBOL:
            raise ValueError("a column cannot be gap-against-gap")
        if y == GAP_SYMBOL:
            state = State.X
        elif x == GAP_SYMBOL:
            state = State.Y
        else:
            state = State.M
        if state == State.M:
            total += dp.matrix.score(x, y)
        else:
            total += dp.gap.extend_score if previous == state else dp.gap.open_score
        previous = state
    return total
