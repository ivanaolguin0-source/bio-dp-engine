"""Substitution matrices: the score for placing residue x opposite residue y.

Scores are integers by design. The standard matrices (BLOSUM, PAM) are
distributed as integer log-odds values, and integer arithmetic makes tie
detection exact. Exact ties are what let the engine report *every*
co-optimal alignment instead of whichever one floating point rounding
happened to favour.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

import numpy as np
import numpy.typing as npt

IntMatrix = npt.NDArray[np.int64]
IndexArray = npt.NDArray[np.intp]

_ASCII_SIZE: Final = 128
GAP_SYMBOL: Final = "-"


class InvalidResidueError(ValueError):
    """A sequence contains a symbol that the scoring alphabet does not define."""

    def __init__(self, label: str, position: int, symbol: str, alphabet: str) -> None:
        self.label = label
        self.position = position
        self.symbol = symbol
        self.alphabet = alphabet
        super().__init__(
            f"{label}: symbol {symbol!r} at position {position} (0-based) "
            f"is not in the scoring alphabet {alphabet!r}"
        )


@dataclass(frozen=True, eq=False)
class SubstitutionMatrix:
    """A square integer scoring table over a residue alphabet.

    Attributes:
        name: Human-readable identifier, carried into DPTrace output.
        alphabet: Residue symbols, one per row/column, stored uppercase.
        scores: ``scores[i, j]`` is the score for ``alphabet[i]`` aligned
            to ``alphabet[j]``. Stored as a read-only ``int64`` copy.

    Input sequences are case-folded before lookup, so soft-masked
    (lowercase) genomic sequence aligns like its uppercase form.
    """

    name: str
    alphabet: str
    scores: IntMatrix
    _lookup: IndexArray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        alphabet = self.alphabet.upper()
        if not alphabet:
            raise ValueError("alphabet must contain at least one symbol")
        if not alphabet.isascii() or not alphabet.isprintable():
            raise ValueError("alphabet must consist of printable ASCII symbols")
        if len(set(alphabet)) != len(alphabet):
            raise ValueError(f"alphabet {alphabet!r} contains duplicate symbols")
        if GAP_SYMBOL in alphabet:
            raise ValueError(f"{GAP_SYMBOL!r} is reserved for gaps and cannot be a residue")

        scores = np.array(self.scores)  # defensive copy
        k = len(alphabet)
        if scores.shape != (k, k):
            raise ValueError(f"scores must have shape ({k}, {k}) to match the alphabet, got {scores.shape}")
        if not np.issubdtype(scores.dtype, np.integer):
            raise TypeError(
                "substitution scores must be integers; integer scoring keeps tie detection exact "
                f"(got dtype {scores.dtype})"
            )
        scores = scores.astype(np.int64)
        scores.setflags(write=False)

        lookup = np.full(_ASCII_SIZE, -1, dtype=np.intp)
        for index, symbol in enumerate(alphabet):
            lookup[ord(symbol)] = index
        lookup.setflags(write=False)

        object.__setattr__(self, "alphabet", alphabet)
        object.__setattr__(self, "scores", scores)
        object.__setattr__(self, "_lookup", lookup)

    # ------------------------------------------------------------------ constructors

    @classmethod
    def match_mismatch(
        cls, match: int, mismatch: int, alphabet: str = "ACGT", name: str | None = None
    ) -> SubstitutionMatrix:
        """Build the simple nucleotide-style table: one score on the diagonal, one off it."""
        k = len(alphabet)
        scores = np.full((k, k), mismatch, dtype=np.int64)
        np.fill_diagonal(scores, match)
        return cls(name or f"match{match:+d}/mismatch{mismatch:+d}", alphabet, scores)

    @classmethod
    def from_ncbi_text(cls, text: str, name: str) -> SubstitutionMatrix:
        """Parse the whitespace-delimited matrix format used by NCBI and EMBOSS.

        Lines starting with ``#`` are comments. The first remaining line lists
        column symbols; each later line is a row symbol followed by one integer
        per column. Rows may appear in any order.
        """
        rows = [
            line.split()
            for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not rows:
            raise ValueError(f"{name}: no matrix data found")
        header, body = rows[0], rows[1:]
        k = len(header)
        by_symbol: dict[str, list[int]] = {}
        for line_no, row in enumerate(body, start=2):
            if len(row) != k + 1:
                raise ValueError(f"{name}: data row {line_no} has {len(row) - 1} values, expected {k}")
            symbol = row[0].upper()
            if symbol in by_symbol:
                raise ValueError(f"{name}: row {symbol!r} appears twice")
            by_symbol[symbol] = [int(value) for value in row[1:]]
        columns = [symbol.upper() for symbol in header]
        if set(by_symbol) != set(columns):
            raise ValueError(f"{name}: row symbols do not match column symbols")
        scores = np.array([by_symbol[symbol] for symbol in columns], dtype=np.int64)
        return cls(name, "".join(columns), scores)

    @classmethod
    def from_ncbi_file(cls, path: str | Path, name: str | None = None) -> SubstitutionMatrix:
        """Load an NCBI-format matrix file; the name defaults to the file stem."""
        path = Path(path)
        return cls.from_ncbi_text(path.read_text(encoding="ascii"), name or path.stem)

    # ------------------------------------------------------------------ lookups

    def encode(self, sequence: str, label: str = "sequence") -> IndexArray:
        """Map a sequence to row/column indices, rejecting unknown symbols.

        Raises:
            InvalidResidueError: naming the first offending symbol and its
                0-based position, so the caller can report it precisely.
        """
        if not sequence.isascii():
            position = next(i for i, ch in enumerate(sequence) if not ch.isascii())
            raise InvalidResidueError(label, position, sequence[position], self.alphabet)
        normalized = sequence.upper()  # length-preserving for ASCII
        codes = np.frombuffer(normalized.encode("ascii"), dtype=np.uint8)
        indices = self._lookup[codes]
        bad = np.flatnonzero(indices < 0)
        if bad.size:
            position = int(bad[0])
            raise InvalidResidueError(label, position, normalized[position], self.alphabet)
        return indices

    def score(self, x: str, y: str) -> int:
        """Score a single residue pair."""
        i, j = self.encode(x + y, "pair")
        return int(self.scores[i, j])
