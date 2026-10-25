"""Hirschberg's linear-space alignment, checked against the full-table engine.

The full-table Needleman-Wunsch engine is itself verified against brute force
and Biopython, so it serves as the reference here. Most tests pass
``base_cells=0`` so the divide step runs all the way down to single
residues; with the default threshold, small inputs would never split at all
and the tests would only be exercising Needleman-Wunsch again.
"""

from __future__ import annotations

import random
import tracemalloc
from typing import Any

import numpy as np
import pytest
from helpers import NW_GOLDEN, load_json, named_matrix, scoring_from_fixture
from hypothesis import given, settings
from hypothesis import strategies as st

from biodp.align import global_score, hirschberg, needleman_wunsch
from biodp.align.hirschberg import last_row
from biodp.scoring import GapPenalty, InvalidResidueError, SubstitutionMatrix

DNA = SubstitutionMatrix.match_mismatch(1, -1)


def is_valid_alignment(aligned_a: str, aligned_b: str, a: str, b: str) -> bool:
    return (
        len(aligned_a) == len(aligned_b)
        and aligned_a.replace("-", "") == a.upper()
        and aligned_b.replace("-", "") == b.upper()
        and all(not (x == "-" and y == "-") for x, y in zip(aligned_a, aligned_b))
    )


# ---------------------------------------------------------------------- golden fixtures

LINEAR = [p for p in NW_GOLDEN if load_json(p)["input"]["gap"]["open"] == load_json(p)["input"]["gap"]["extend"]]
AFFINE = [p for p in NW_GOLDEN if p not in LINEAR]


@pytest.mark.parametrize("path", LINEAR, ids=lambda p: p.stem)
def test_linear_fixture_matches_full_table(path: Any) -> None:
    """The build-order checkpoint: same optimal score as Needleman-Wunsch on every fixture."""
    fixture = load_json(path)
    matrix, gap = scoring_from_fixture(fixture)
    a, b = fixture["input"]["seq_a"], fixture["input"]["seq_b"]
    result = hirschberg(a, b, matrix=matrix, gap=gap, base_cells=0)
    assert result.score == fixture["expected"]["score"]
    assert global_score(a, b, matrix=matrix, gap=gap) == fixture["expected"]["score"]
    if fixture["expected"]["alignments"] is not None:
        assert [result.aligned_a, result.aligned_b] in fixture["expected"]["alignments"]


@pytest.mark.parametrize("path", AFFINE, ids=lambda p: p.stem)
def test_affine_fixture_is_refused(path: Any) -> None:
    fixture = load_json(path)
    matrix, gap = scoring_from_fixture(fixture)
    with pytest.raises(ValueError, match="Myers-Miller"):
        hirschberg(fixture["input"]["seq_a"], fixture["input"]["seq_b"], matrix=matrix, gap=gap)


def test_fixture_split_is_real() -> None:
    assert len(LINEAR) >= 4 and len(AFFINE) >= 1


# ---------------------------------------------------------------------- properties


@settings(max_examples=300, deadline=None)
@given(
    a=st.text(alphabet="ACG", max_size=12),
    b=st.text(alphabet="ACG", max_size=12),
    match=st.integers(min_value=0, max_value=5),
    mismatch=st.integers(min_value=-6, max_value=0),
    g=st.integers(min_value=-6, max_value=0),
)
def test_result_is_one_of_the_cooptimal_alignments(a: str, b: str, match: int, mismatch: int, g: int) -> None:
    matrix = SubstitutionMatrix.match_mismatch(match, mismatch)
    gap = GapPenalty.linear(g)
    reference = needleman_wunsch(a, b, matrix=matrix, gap=gap)
    result = hirschberg(a, b, matrix=matrix, gap=gap, base_cells=0)
    assert result.score == reference.score
    assert is_valid_alignment(result.aligned_a, result.aligned_b, a, b)
    if reference.optimal_count <= 2000:
        optimal = {(x.aligned_a, x.aligned_b) for x in reference.alignments}
        assert (result.aligned_a, result.aligned_b) in optimal


@settings(max_examples=100, deadline=None)
@given(a=st.text(alphabet="ACGT", max_size=60), b=st.text(alphabet="ACGT", max_size=60), cells=st.sampled_from([0, 16, 4096]))
def test_base_threshold_never_changes_the_score(a: str, b: str, cells: int) -> None:
    gap = GapPenalty.linear(-2)
    expected = needleman_wunsch(a, b, matrix=DNA, gap=gap).score
    assert hirschberg(a, b, matrix=DNA, gap=gap, base_cells=cells).score == expected


@settings(max_examples=100, deadline=None)
@given(a=st.text(alphabet="ACGT", max_size=40), b=st.text(alphabet="ACGT", max_size=40))
def test_last_row_equals_full_table_bottom_row(a: str, b: str) -> None:
    gap = GapPenalty.linear(-1)
    dp = needleman_wunsch(a, b, matrix=DNA, gap=gap).dp
    full_bottom = dp.scores[:, len(a), :].max(axis=0)
    row = last_row(DNA.encode(a), DNA.encode(b), DNA.scores, -1)
    np.testing.assert_array_equal(row, full_bottom)


def test_path_runs_corner_to_corner_in_unit_steps() -> None:
    result = hirschberg("GATTACAGATTACA", "GCATGCUGCA".replace("U", "T"), matrix=DNA, gap=GapPenalty.linear(-1), base_cells=0)
    path = result.path
    assert path[0] == (0, 0) and path[-1] == (14, 10)
    for (i0, j0), (i1, j1) in zip(path, path[1:]):
        assert (i1 - i0, j1 - j0) in {(1, 1), (1, 0), (0, 1)}


# ---------------------------------------------------------------------- larger inputs


@pytest.mark.parametrize("seed", range(10))
def test_protein_blosum62_agrees_with_full_table(seed: int) -> None:
    rng = random.Random(seed)
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    a = "".join(rng.choice(alphabet) for _ in range(rng.randint(50, 200)))
    b = "".join(rng.choice(alphabet) for _ in range(rng.randint(50, 200)))
    blosum = named_matrix("BLOSUM62")
    gap = GapPenalty.linear(-4)
    expected = needleman_wunsch(a, b, matrix=blosum, gap=gap).score
    assert hirschberg(a, b, matrix=blosum, gap=gap).score == expected
    assert global_score(a, b, matrix=blosum, gap=gap) == expected


def test_asymmetric_matrix_survives_the_swap_in_global_score() -> None:
    # A deliberately asymmetric matrix: aligning A over C is not the same as C over A.
    matrix = SubstitutionMatrix("asym", "AC", np.array([[3, -5], [1, 2]]))
    gap = GapPenalty.linear(-2)
    for a, b in [("AACCA", "CA"), ("C", "AAAAC"), ("ACAC", "CACACA")]:
        expected = needleman_wunsch(a, b, matrix=matrix, gap=gap).score
        assert global_score(a, b, matrix=matrix, gap=gap) == expected
        assert hirschberg(a, b, matrix=matrix, gap=gap, base_cells=0).score == expected


def test_memory_stays_linear_where_the_full_table_would_not_fit() -> None:
    # 3,000 x 3,000 is 9 million cells: about 243 MB of tables for the full
    # engine. Hirschberg should need a tiny fraction of that.
    rng = random.Random(7)
    a = "".join(rng.choice("ACGT") for _ in range(3000))
    b = list(a)
    for _ in range(300):
        b[rng.randrange(len(b))] = rng.choice("ACGT")
    gap = GapPenalty.linear(-2)
    tracemalloc.start()
    result = hirschberg(a, "".join(b), matrix=DNA, gap=gap)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert result.score == global_score(a, "".join(b), matrix=DNA, gap=gap)
    full_table_bytes = 27 * 3001 * 3001
    assert peak < full_table_bytes / 50, f"peak {peak:,} bytes is not linear-space"


# ---------------------------------------------------------------------- edge cases


@pytest.mark.parametrize(("a", "b"), [("", ""), ("", "ACG"), ("ACG", ""), ("A", "A"), ("A", "C")])
def test_edge_cases(a: str, b: str) -> None:
    gap = GapPenalty.linear(-1)
    result = hirschberg(a, b, matrix=DNA, gap=gap, base_cells=0)
    assert result.score == needleman_wunsch(a, b, matrix=DNA, gap=gap).score
    assert is_valid_alignment(result.aligned_a, result.aligned_b, a, b)


def test_affine_is_refused_with_guidance() -> None:
    with pytest.raises(ValueError, match="Myers-Miller"):
        hirschberg("ACGT", "AGT", matrix=DNA, gap=GapPenalty(-5, -1))
    with pytest.raises(ValueError, match="Myers-Miller"):
        global_score("ACGT", "AGT", matrix=DNA, gap=GapPenalty(-5, -1))


def test_invalid_residue_is_reported() -> None:
    with pytest.raises(InvalidResidueError, match="seq_b"):
        hirschberg("ACGT", "ACXT", matrix=DNA, gap=GapPenalty.linear(-1))


def test_negative_base_cells_rejected() -> None:
    with pytest.raises(ValueError):
        hirschberg("A", "A", matrix=DNA, gap=GapPenalty.linear(-1), base_cells=-1)
