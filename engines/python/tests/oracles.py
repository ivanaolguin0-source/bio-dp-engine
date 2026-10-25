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


# ---------------------------------------------------------------------- RNA folding oracles

RnaPairs = tuple[tuple[int, int], ...]
WC = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
GU = {("G", "U"), ("U", "G")}


def rna_can_pair(x: str, y: str, wobble: bool) -> bool:
    return (x, y) in WC or (wobble and (x, y) in GU)


def all_rna_structures(seq: str, min_loop: int, wobble: bool) -> list[RnaPairs]:
    """Every pseudoknot-free secondary structure, built from the FIRST base.

    The engine decomposes on the last base of each interval; this oracle
    decides the first base's partner instead, so the two share no recurrence.
    """
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def build(i: int, j: int) -> tuple[RnaPairs, ...]:
        if j < i:
            return ((),)
        out: list[RnaPairs] = list(build(i + 1, j))  # i unpaired
        for k in range(i + min_loop + 1, j + 1):
            if rna_can_pair(seq[i], seq[k], wobble):
                for inner in build(i + 1, k - 1):
                    for right in build(k + 1, j):
                        out.append(((i, k),) + inner + right)
        return tuple(out)

    return [tuple(sorted(s)) for s in build(0, len(seq) - 1)]


def max_pair_structures(seq: str, min_loop: int, wobble: bool) -> tuple[int, frozenset[RnaPairs]]:
    structures = all_rna_structures(seq, min_loop, wobble)
    best = max(len(s) for s in structures)
    return best, frozenset(s for s in structures if len(s) == best)


def classic_nussinov_score(seq: str, min_loop: int, wobble: bool) -> int:
    """The four-case textbook recurrence (Durbin et al. 1998), score only, plain loops."""
    n = len(seq)
    if n == 0:
        return 0
    N = [[0] * n for _ in range(n)]
    for d in range(min_loop + 1, n):
        for i in range(n - d):
            j = i + d
            best = max(N[i + 1][j], N[i][j - 1])
            if rna_can_pair(seq[i], seq[j], wobble):
                best = max(best, N[i + 1][j - 1] + 1)
            for k in range(i + 1, j):
                best = max(best, N[i][k] + N[k + 1][j])
            N[i][j] = best
    return N[0][n - 1]


def is_valid_structure(seq: str, pairs: RnaPairs, min_loop: int, wobble: bool) -> bool:
    used: set[int] = set()
    for i, j in pairs:
        if not (0 <= i < j < len(seq)) or i in used or j in used:
            return False
        if j - i - 1 < min_loop or not rna_can_pair(seq[i], seq[j], wobble):
            return False
        used.update((i, j))
    for a, b in pairs:
        for c, d in pairs:
            if a < c < b < d:
                return False  # crossing pairs form a pseudoknot
    return True
