"""Nussinov RNA secondary structure prediction: maximise the number of base pairs.

A secondary structure is a set of base pairs (i, j) in which every base pairs
at most once and no two pairs cross (i < k < j < l is forbidden for pairs
(i, j) and (k, l)). Forbidding crossings, which excludes pseudoknots, is what
makes the problem decompose into independent intervals and therefore solvable
by dynamic programming.

What this is and is not
-----------------------
Nussinov maximises a pair *count*. Real RNA folds toward minimum free energy,
where stacking between adjacent pairs and the cost of loops dominate. Energy
models (Zuker's algorithm, as in the ViennaRNA package) are what practical
predictors use. Nussinov is the cleanest example of interval DP and a
standard teaching step on the way there; treat its output as a toy model.

Recurrence
----------
Durbin et al. (1998) present Nussinov with four cases: i unpaired, j unpaired,
i pairs with j, or a bifurcation into two sub-intervals. That form finds the
right optimum, but the same structure can usually be built through several of
those cases, so counting derivations overcounts structures. This engine uses
the standard *unambiguous* equivalent, which decomposes on the last base:

    N(i, j) = max( N(i, j-1),                                 # j unpaired
                   max over k of N(i, k-1) + 1 + N(k+1, j-1) )  # j pairs with k

where k ranges over i <= k < j with a valid pair (k, j) enclosing a hairpin
loop of at least ``min_loop`` bases. Every structure has exactly one
derivation, so the number of optimal derivations is exactly the number of
optimal structures. Tests confirm that both forms give the same optimum.
"""

from __future__ import annotations

import random
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import Final

import numpy as np
import numpy.typing as npt

RNA_ALPHABET: Final = "ACGU"
WATSON_CRICK: Final = frozenset({("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")})
WOBBLE: Final = frozenset({("G", "U"), ("U", "G")})
DEFAULT_MIN_LOOP: Final = 3
# n^2 cells, each holding an int64 score and an arbitrary-precision count.
DEFAULT_MAX_LENGTH: Final = 2000

IntGrid = npt.NDArray[np.int64]
Pair = tuple[int, int]


def pair_label(x: str, y: str) -> str:
    """'GC', 'AU' or 'GU' for a valid pair in either orientation."""
    return {frozenset("GC"): "GC", frozenset("AU"): "AU", frozenset("GU"): "GU"}[frozenset((x, y))]


class InvalidRnaError(ValueError):
    """The sequence contains a symbol other than A, C, G, U (case-insensitive)."""


@dataclass(frozen=True)
class FoldingRules:
    """Which pairs are allowed.

    Attributes:
        min_loop: Minimum number of unpaired bases enclosed by a hairpin, so a
            pair (i, j) requires j - i - 1 >= min_loop. Three is the usual
            choice; zero allows sterically impossible hairpins.
        wobble: Whether G-U pairs are allowed in addition to Watson-Crick pairs.
    """

    min_loop: int = DEFAULT_MIN_LOOP
    wobble: bool = True

    def __post_init__(self) -> None:
        if isinstance(self.min_loop, bool) or not isinstance(self.min_loop, int) or self.min_loop < 0:
            raise ValueError("min_loop must be a non-negative integer")

    def can_pair(self, x: str, y: str) -> bool:
        return (x, y) in WATSON_CRICK or (self.wobble and (x, y) in WOBBLE)


def validate_rna(sequence: str) -> str:
    """Uppercase and check a sequence; T is rejected rather than silently converted."""
    seq = sequence.upper()
    for position, base in enumerate(seq):
        if base not in RNA_ALPHABET:
            hint = " (DNA 'T' is not converted automatically; replace it with 'U')" if base == "T" else ""
            raise InvalidRnaError(f"symbol {base!r} at position {position} (0-based) is not an RNA base{hint}")
    return seq


@dataclass(frozen=True)
class Structure:
    """One secondary structure.

    Attributes:
        sequence: The folded sequence, uppercase.
        pairs: Base pairs (i, j), 0-based, i < j, sorted by i.
        visited: Table cells (i, j) visited by the traceback that produced
            this structure, in order. The visualiser animates these.
    """

    sequence: str
    pairs: tuple[Pair, ...]
    visited: tuple[Pair, ...] = ()

    @property
    def dot_bracket(self) -> str:
        """Standard notation: '(' and ')' for paired bases, '.' for unpaired."""
        chars = ["."] * len(self.sequence)
        for i, j in self.pairs:
            chars[i], chars[j] = "(", ")"
        return "".join(chars)

    @property
    def pair_count(self) -> int:
        return len(self.pairs)

    def pair_types(self) -> dict[str, int]:
        """Counts of GC, AU and GU pairs, regardless of which base comes first."""
        counts = {"GC": 0, "AU": 0, "GU": 0}
        for i, j in self.pairs:
            counts[pair_label(self.sequence[i], self.sequence[j])] += 1
        return counts


@dataclass(frozen=True, eq=False)
class NussinovResult:
    """Completed tables plus lazy access to every optimal structure.

    ``scores[i, j]`` is the maximum number of pairs in the interval i..j
    (inclusive) for i <= j, and -1 below the diagonal, where it has no meaning.
    """

    sequence: str
    rules: FoldingRules
    scores: IntGrid
    _counts: list[list[int]]

    @property
    def n(self) -> int:
        return len(self.sequence)

    @property
    def score(self) -> int:
        """Maximum number of base pairs for the whole sequence."""
        return self.N(0, self.n - 1)

    @property
    def optimal_count(self) -> int:
        """Exact number of distinct optimal structures."""
        return self.C(0, self.n - 1)

    # Interval accessors that treat empty intervals (j < i) as 0 pairs, 1 structure.
    def N(self, i: int, j: int) -> int:
        return 0 if j < i else int(self.scores[i, j])

    def C(self, i: int, j: int) -> int:
        return 1 if j < i else self._counts[i][j]

    def choices(self, i: int, j: int) -> tuple[bool, list[int]]:
        """Optimal decisions at (i, j): whether 'j unpaired' is optimal, and every optimal partner k of j."""
        best = self.N(i, j)
        unpaired = j > i and self.N(i, j - 1) == best
        partners = [
            k
            for k in range(i, j - self.rules.min_loop)
            if self.rules.can_pair(self.sequence[k], self.sequence[j])
            and self.N(i, k - 1) + 1 + self.N(k + 1, j - 1) == best
        ]
        return unpaired, partners

    def visited_for(self, pairs: tuple[Pair, ...]) -> tuple[Pair, ...]:
        """Cells the traceback visits to produce a given optimal structure, in order."""
        partner = {j: i for i, j in pairs}
        return self._traceback(lambda i, j: partner[j] if partner.get(j, -1) >= i else None).visited

    def optimal_cells(self) -> frozenset[Pair]:
        """Every table cell that lies on the traceback of at least one optimal structure."""
        seen: set[Pair] = set()
        stack = [(0, self.n - 1)] if self.n else []
        while stack:
            i, j = stack.pop()
            if j < i or (i, j) in seen:
                continue
            seen.add((i, j))
            if self.N(i, j) == 0:
                continue
            unpaired, partners = self.choices(i, j)
            if unpaired:
                stack.append((i, j - 1))
            for k in partners:
                stack.extend([(i, k - 1), (k + 1, j - 1)])
        return frozenset(seen)

    def best(self) -> Structure:
        """Canonical optimal structure: at each step prefer 'j unpaired', then the smallest partner k."""
        def canonical(i: int, j: int) -> int | None:
            unpaired, partners = self.choices(i, j)
            return None if unpaired else partners[0]

        return self._traceback(canonical)

    def sample(self, rng: random.Random | None = None) -> Structure:
        """One optimal structure drawn uniformly at random, using the structure counts as weights."""
        rng = rng or random.Random()

        def pick_at(i: int, j: int) -> int | None:
            unpaired, partners = self.choices(i, j)
            options: list[int | None] = ([None] if unpaired else []) + list(partners)
            weights = [self.C(i, j - 1) if k is None else self.C(i, k - 1) * self.C(k + 1, j - 1) for k in options]
            ticket = rng.randrange(sum(weights))
            for k, w in zip(options, weights):
                if ticket < w:
                    return k
                ticket -= w
            raise AssertionError("unreachable")

        return self._traceback(pick_at)

    def structures(self) -> Iterator[Structure]:
        """Every optimal structure exactly once, canonical first, generated lazily.

        Nested generators recurse once per nested interval, so this is meant
        for the sequence lengths the visualiser shows (a few hundred bases).
        """

        def gen(i: int, j: int) -> Iterator[tuple[Pair, ...]]:
            if j < i or self.N(i, j) == 0:
                yield ()
                return
            unpaired, partners = self.choices(i, j)
            if unpaired:
                yield from gen(i, j - 1)
            for k in partners:
                for left in gen(i, k - 1):
                    for inner in gen(k + 1, j - 1):
                        yield left + ((k, j),) + inner

        for pairs in gen(0, self.n - 1):
            yield Structure(self.sequence, tuple(sorted(pairs)))

    def take(self, limit: int) -> list[Structure]:
        out: list[Structure] = []
        for s in self.structures():
            if len(out) == limit:
                break
            out.append(s)
        return out

    def _traceback(self, pick: Callable[[int, int], int | None]) -> Structure:
        """Iterative traceback. ``pick(i, j)`` returns j's partner k, or None for 'j unpaired'."""
        pairs: list[Pair] = []
        visited: list[Pair] = []
        stack = [(0, self.n - 1)] if self.n else []
        while stack:
            i, j = stack.pop()
            if j < i:
                continue
            visited.append((i, j))
            if self.N(i, j) == 0:
                continue
            k = pick(i, j)
            if k is None:
                stack.append((i, j - 1))
            else:
                pairs.append((k, j))
                # Inner interval first, then the left one; the stack pops left first.
                stack.append((k + 1, j - 1))
                stack.append((i, k - 1))
        return Structure(self.sequence, tuple(sorted(pairs)), tuple(visited))


def nussinov(
    sequence: str,
    *,
    rules: FoldingRules | None = None,
    max_length: int = DEFAULT_MAX_LENGTH,
) -> NussinovResult:
    """Fold an RNA sequence by maximising base pairs.

    Args:
        sequence: RNA bases (A, C, G, U), case-insensitive.
        rules: Pairing rules; defaults to Watson-Crick plus G-U wobble with a
            minimum hairpin loop of 3.
        max_length: Refuse longer inputs (time is O(n^3), memory O(n^2)).

    Returns:
        The filled tables, exposing the optimal pair count, the exact number of
        optimal structures, the canonical structure, enumeration and sampling.
    """
    rules = rules or FoldingRules()
    seq = validate_rna(sequence)
    n = len(seq)
    if n > max_length:
        raise ValueError(f"sequence length {n} exceeds max_length {max_length}; Nussinov is O(n^3) time and O(n^2) memory")

    codes = np.array([RNA_ALPHABET.index(b) for b in seq], dtype=np.intp)
    allowed = np.zeros((4, 4), dtype=bool)
    for x in RNA_ALPHABET:
        for y in RNA_ALPHABET:
            allowed[RNA_ALPHABET.index(x), RNA_ALPHABET.index(y)] = rules.can_pair(x, y)
    pairable = allowed[codes[:, None], codes[None, :]] if n else np.zeros((0, 0), dtype=bool)

    # S[i, j] for i <= j; padded so that S[i, i-1] (empty interval) reads as 0.
    # Row/column offset: P[i, j + 1] stores N(i, j); P[i, i] stores N(i, i-1) = 0.
    P = np.zeros((n + 1, n + 1), dtype=np.int64)
    counts: list[list[int]] = [[1] * (n + 1) for _ in range(n + 1)]  # same offset; empty = 1
    h = rules.min_loop

    # Fill by interval length so every sub-interval is ready before it is read.
    for d in range(0, n):
        for i in range(0, n - d):
            j = i + d
            best = P[i, j]  # N(i, j-1): j unpaired
            count_unpaired = counts[i][j]
            ks = np.arange(i, max(i, j - h))
            ks = ks[pairable[ks, j]] if ks.size else ks
            if ks.size:
                # N(i, k-1) is P[i, k]; N(k+1, j-1) is P[k+1, j].
                cand = P[i, ks] + 1 + P[ks + 1, j]
                top = int(cand.max())
            else:
                top = -1
            value = max(int(best), top)
            total = count_unpaired if int(best) == value else 0
            if top == value:
                for k in ks[cand == value].tolist():
                    total += counts[i][k] * counts[k + 1][j]
            P[i, j + 1] = value
            counts[i][j + 1] = total

    scores = np.full((n, n), -1, dtype=np.int64)
    out_counts = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            scores[i, j] = P[i, j + 1]
            out_counts[i][j] = counts[i][j + 1]
    scores.setflags(write=False)
    return NussinovResult(seq, rules, scores, out_counts)


def base_pair_metrics(predicted: Structure | tuple[Pair, ...], reference: tuple[Pair, ...]) -> dict[str, float]:
    """Sensitivity (recall) and positive predictive value (precision) of predicted pairs against a reference."""
    pred = set(predicted.pairs if isinstance(predicted, Structure) else predicted)
    ref = set(reference)
    hit = len(pred & ref)
    return {
        "true_positives": float(hit),
        "sensitivity": hit / len(ref) if ref else float("nan"),
        "ppv": hit / len(pred) if pred else float("nan"),
    }


def parse_dot_bracket(structure: str) -> tuple[Pair, ...]:
    """Pairs from a dot-bracket string; raises on unbalanced brackets."""
    stack: list[int] = []
    pairs: list[Pair] = []
    for j, ch in enumerate(structure):
        if ch == "(":
            stack.append(j)
        elif ch == ")":
            if not stack:
                raise ValueError(f"unbalanced ')' at position {j}")
            pairs.append((stack.pop(), j))
        elif ch != ".":
            raise ValueError(f"unexpected character {ch!r} at position {j}")
    if stack:
        raise ValueError(f"unbalanced '(' at position {stack[-1]}")
    return tuple(sorted(pairs))
