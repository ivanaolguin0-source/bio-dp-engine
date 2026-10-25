from __future__ import annotations

import numpy as np
import pytest
from helpers import named_matrix

from biodp.scoring import GapPenalty, InvalidResidueError, SubstitutionMatrix


class TestSubstitutionMatrix:
    def test_encode_is_case_insensitive(self) -> None:
        m = SubstitutionMatrix.match_mismatch(1, -1)
        assert m.encode("acgt").tolist() == m.encode("ACGT").tolist() == [0, 1, 2, 3]

    def test_encode_reports_first_bad_symbol_and_position(self) -> None:
        m = SubstitutionMatrix.match_mismatch(1, -1)
        with pytest.raises(InvalidResidueError) as info:
            m.encode("ACGNT", "seq_a")
        assert (info.value.label, info.value.position, info.value.symbol) == ("seq_a", 3, "N")

    def test_encode_rejects_non_ascii_with_correct_position(self) -> None:
        m = SubstitutionMatrix.match_mismatch(1, -1)
        with pytest.raises(InvalidResidueError) as info:
            m.encode("ACß")
        assert info.value.position == 2

    def test_encode_empty_sequence(self) -> None:
        assert SubstitutionMatrix.match_mismatch(1, -1).encode("").size == 0

    def test_scores_are_read_only_copies(self) -> None:
        source = np.array([[1, -1], [-1, 1]])
        m = SubstitutionMatrix("tiny", "AC", source)
        source[0, 0] = 99
        assert m.score("A", "A") == 1
        with pytest.raises(ValueError):
            m.scores[0, 0] = 5

    @pytest.mark.parametrize(
        ("alphabet", "scores", "error"),
        [
            ("AA", [[1, 0], [0, 1]], ValueError),
            ("A-", [[1, 0], [0, 1]], ValueError),
            ("AC", [[1, 0, 0], [0, 1, 0]], ValueError),
            ("AC", [[1.5, 0.0], [0.0, 1.0]], TypeError),
            ("", [], ValueError),
        ],
    )
    def test_rejects_malformed_matrices(self, alphabet: str, scores: list[list[float]], error: type[Exception]) -> None:
        with pytest.raises(error):
            SubstitutionMatrix("bad", alphabet, np.array(scores))

    def test_blosum62_loader_matches_biopython(self) -> None:
        substitution_matrices = pytest.importorskip("Bio.Align.substitution_matrices")
        ours = named_matrix("BLOSUM62")
        theirs = substitution_matrices.load("BLOSUM62")
        assert ours.alphabet == theirs.alphabet
        for x in ours.alphabet:
            for y in ours.alphabet:
                assert ours.score(x, y) == int(theirs[x][y])

    def test_ncbi_parser_rejects_ragged_rows(self) -> None:
        with pytest.raises(ValueError, match="row 3"):
            SubstitutionMatrix.from_ncbi_text("  A C\nA 1 0\nC 0\n", "ragged")


class TestGapPenalty:
    def test_affine_length_scoring(self) -> None:
        gap = GapPenalty(-10, -1)
        assert [gap.score(k) for k in range(4)] == [0, -10, -11, -12]

    def test_linear_is_special_case(self) -> None:
        assert GapPenalty.linear(-2).is_linear
        assert not GapPenalty(-3, -1).is_linear

    @pytest.mark.parametrize(("o", "e"), [(1, -1), (-1, 2)])
    def test_rejects_positive_scores(self, o: int, e: int) -> None:
        with pytest.raises(ValueError):
            GapPenalty(o, e)

    def test_rejects_non_integers(self) -> None:
        with pytest.raises(TypeError):
            GapPenalty(-1.5, -1)  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            GapPenalty(True, -1)
