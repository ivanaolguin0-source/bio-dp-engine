"""Independent test oracles. Nothing here imports biodp.

The exhaustive oracle enumerates every possible alignment of two short
sequences, scores each column by column, and keeps the best. It shares no
code and no recurrence with the engine, so agreement between the two is
real evidence rather than a program agreeing with itself.

Cost grows with the Delannoy numbers (48,639 alignments for 7 x 7), so keep
inputs to roughly 10 residues or fewer.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass

Column = tuple[str, str]
ScoreFn = Callable[[str, str], int]


def all_alignments(a: str, b: str) -> Iterator[tuple[Column, ...]]:
    """Every column sequence that spells ``a`` on top and ``b`` underneath."""

    def extend(i: int, j: int, prefix: tuple[Column, ...]) -> Iterator[tuple[Column, ...]]:
        if i == len(a) and j == len(b):
            yield prefix
            return
        if i < len(a) and j < len(b):
            yield from extend(i + 1, j + 1, prefix + ((a[i], b[j]),))
        if i < len(a):
            yield from extend(i + 1, j, prefix + ((a[i], "-"),))
        if j < len(b):
            yield from extend(i, j + 1, prefix + (("-", b[j]),))

    yield from extend(0, 0, ())


def score_columns(columns: tuple[Column, ...], sub: ScoreFn, gap_open: int, gap_extend: int) -> int:
    """Affine scoring: each maximal run of same-kind gap columns is one gap."""
    total = 0
    previous_kind = "pair"
    for x, y in columns:
        kind = "gap_in_b" if y == "-" else "gap_in_a" if x == "-" else "pair"
        if kind == "pair":
            total += sub(x, y)
        else:
            total += gap_extend if kind == previous_kind else gap_open
        previous_kind = kind
    return total


@dataclass(frozen=True)
class OracleResult:
    score: int
    alignments: frozenset[tuple[str, str]]


def exhaustive_global(a: str, b: str, sub: ScoreFn, gap_open: int, gap_extend: int) -> OracleResult:
    best: int | None = None
    winners: set[tuple[str, str]] = set()
    for columns in all_alignments(a, b):
        s = score_columns(columns, sub, gap_open, gap_extend)
        if best is None or s > best:
            best, winners = s, set()
        if s == best:
            winners.add(("".join(c[0] for c in columns), "".join(c[1] for c in columns)))
    assert best is not None
    return OracleResult(best, frozenset(winners))


def parse_ncbi_matrix(text: str) -> dict[tuple[str, str], int]:
    """Minimal standalone reader for NCBI-format matrices."""
    rows = [ln.split() for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    header = rows[0]
    return {(row[0], col): int(v) for row in rows[1:] for col, v in zip(header, row[1:])}
