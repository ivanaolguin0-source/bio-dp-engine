"""Gap scoring models."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GapPenalty:
    """Affine gap scoring.

    A gap of length ``L >= 1`` scores ``open_score + (L - 1) * extend_score``.
    Both values are scores, not costs, so they are zero or negative.

    This is the same convention Biopython's ``PairwiseAligner`` uses for
    ``open_gap_score`` / ``extend_gap_score``, which keeps cross-validation
    against it straightforward. Note that some tools (EMBOSS ``needle``, for
    one) phrase the parameters as positive penalties; convert before comparing.

    Linear gap scoring is the special case ``open_score == extend_score``.
    """

    open_score: int
    extend_score: int

    def __post_init__(self) -> None:
        for label, value in (("open_score", self.open_score), ("extend_score", self.extend_score)):
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"{label} must be an int, got {type(value).__name__}")
            if value > 0:
                raise ValueError(f"{label} must be <= 0 (it is a score, not a penalty); got {value}")

    @classmethod
    def linear(cls, score: int) -> GapPenalty:
        """Every gap position scores the same amount."""
        return cls(score, score)

    @property
    def is_linear(self) -> bool:
        return self.open_score == self.extend_score

    def score(self, length: int) -> int:
        """Total score of one contiguous gap of ``length`` positions."""
        if length < 0:
            raise ValueError("gap length cannot be negative")
        if length == 0:
            return 0
        return self.open_score + (length - 1) * self.extend_score
