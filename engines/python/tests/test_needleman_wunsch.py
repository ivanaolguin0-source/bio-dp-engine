"""Correctness of the Needleman-Wunsch engine.

Evidence comes from three directions that share no code with the engine:
hand-chosen golden fixtures (expected values from exhaustive enumeration),
property-based comparison against that same exhaustive oracle on random
inputs, and comparison against Biopython on inputs too long to enumerate.
"""

from __future__ import annotations

import math
import random
from typing import Any

import numpy as np
import pytest
from helpers import NW_GOLDEN, load_json, named_matrix, scoring_from_fixture
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from oracles import exhaustive_global

from biodp.align import fill_matrices, needleman_wunsch, rescore
from biodp.scoring import GapPenalty, InvalidResidueError, SubstitutionMatrix

DNA = SubstitutionMatrix.match_mismatch(1, -1)

# ---------------------------------------------------------------------- golden fixtures


@pytest.mark.parametrize("path", NW_GOLDEN, ids=lambda p: p.stem)
@pytest.mark.parametrize("kernel", ["vectorized", "reference"])
def test_golden_fixture(path: Any, kernel: Any) -> None:
    fixture = load_json(path)
    matrix, gap = scoring_from_fixture(fixture)
    result = needleman_wunsch(
        fixture["input"]["seq_a"], fixture["input"]["seq_b"], matrix=matrix, gap=gap, kernel=kernel
    )
    expected = fixture["expected"]
    assert result.score == expected["score"]
    assert result.optimal_count == int(expected["optimal_count"])
    if expected["alignments"] is not None:
        got = sorted([a.aligned_a, a.aligned_b] for a in result.alignments)
        assert got == expected["alignments"]


def test_golden_fixtures_exist() -> None:
    assert len(NW_GOLDEN) >= 3, "golden fixtures missing; run scripts/generate_nw_golden.py"


# ---------------------------------------------------------------------- edge cases


def test_both_empty() -> None:
    result = needleman_wunsch("", "", matrix=DNA, gap=GapPenalty.linear(-1))
    assert (result.score, result.optimal_count) == (0, 1)
    assert result.best.aligned_a == result.best.aligned_b == ""
    assert result.best.path == ((0, 0),)


@pytest.mark.parametrize(("a", "b"), [("", "ACGT"), ("ACGT", "")])
def test_one_empty_is_a_single_affine_gap(a: str, b: str) -> None:
    gap = GapPenalty(-5, -2)
    result = needleman_wunsch(a, b, matrix=DNA, gap=gap)
    assert result.score == gap.score(4)
    assert result.optimal_count == 1


def test_identical_sequences_align_without_gaps() -> None:
    result = needleman_wunsch("ACGTACGT", "acgtacgt", matrix=DNA, gap=GapPenalty(-3, -1))
    assert result.best.cigar() == "8="
    assert result.optimal_count == 1


def test_invalid_residue_is_reported_with_sequence_label() -> None:
    with pytest.raises(InvalidResidueError, match="seq_b"):
        needleman_wunsch("ACGT", "ACXT", matrix=DNA, gap=GapPenalty.linear(-1))


def test_size_guard() -> None:
    with pytest.raises(ValueError, match="Hirschberg"):
        needleman_wunsch("A" * 100, "A" * 100, matrix=DNA, gap=GapPenalty.linear(-1), max_cells=10_000)


def test_unknown_kernel() -> None:
    with pytest.raises(ValueError, match="kernel"):
        fill_matrices("A", "A", matrix=DNA, gap=GapPenalty.linear(-1), kernel="gpu")  # type: ignore[arg-type]


# ---------------------------------------------------------------------- multiple optimal paths


def test_count_is_exact_beyond_machine_integers() -> None:
    # Every optimal alignment pairs all 100 short-sequence residues and gaps the
    # other 100, so the count is the number of ways to choose the gapped ones.
    result = needleman_wunsch("A" * 200, "A" * 100, matrix=DNA, gap=GapPenalty.linear(-1))
    assert result.optimal_count == math.comb(200, 100)
    assert result.optimal_count > 2**63


def test_canonical_tie_break_left_aligns_gaps() -> None:
    result = needleman_wunsch("AAA", "AA", matrix=DNA, gap=GapPenalty.linear(-1))
    assert (result.best.aligned_a, result.best.aligned_b) == ("AAA", "-AA")


def test_enumeration_is_lazy() -> None:
    result = needleman_wunsch("A" * 200, "A" * 100, matrix=DNA, gap=GapPenalty.linear(-1))
    first = result.alignments.take(5)
    assert len({(a.aligned_a, a.aligned_b) for a in first}) == 5


def test_sampling_is_uniform() -> None:
    # 6 optimal alignments; with 6000 draws each should appear about 1000 times.
    result = needleman_wunsch("AAAA", "AA", matrix=DNA, gap=GapPenalty.linear(-1))
    assert result.optimal_count == math.comb(4, 2)
    rng = random.Random(20261025)
    tallies: dict[str, int] = {}
    for _ in range(6000):
        key = result.alignments.sample(rng).aligned_b
        tallies[key] = tallies.get(key, 0) + 1
    assert len(tallies) == 6
    # Chi-square with 5 degrees of freedom; 20.5 is roughly the 0.999 quantile.
    chi2 = sum((observed - 1000) ** 2 / 1000 for observed in tallies.values())
    assert chi2 < 20.5


# ---------------------------------------------------------------------- property tests

small_seq = st.text(alphabet="ACG", min_size=0, max_size=5)
scores = st.integers(min_value=-6, max_value=0)


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    a=small_seq,
    b=small_seq,
    match=st.integers(min_value=0, max_value=5),
    mismatch=scores,
    g_open=scores,
    g_ext=scores,
)
def test_matches_exhaustive_oracle(a: str, b: str, match: int, mismatch: int, g_open: int, g_ext: int) -> None:
    matrix = SubstitutionMatrix.match_mismatch(match, mismatch)
    gap = GapPenalty(g_open, g_ext)
    oracle = exhaustive_global(a, b, matrix.score, g_open, g_ext)
    result = needleman_wunsch(a, b, matrix=matrix, gap=gap)

    assert result.score == oracle.score
    assert result.optimal_count == len(oracle.alignments)
    enumerated = [(x.aligned_a, x.aligned_b) for x in result.alignments]
    assert len(enumerated) == len(set(enumerated)), "enumeration yielded a duplicate"
    assert set(enumerated) == oracle.alignments


@settings(max_examples=150, deadline=None)
@given(
    a=st.text(alphabet="ACGT", max_size=30),
    b=st.text(alphabet="ACGT", max_size=30),
    g_open=scores,
    g_ext=scores,
)
def test_kernels_are_bit_identical(a: str, b: str, g_open: int, g_ext: int) -> None:
    gap = GapPenalty(g_open, g_ext)
    fast = fill_matrices(a, b, matrix=DNA, gap=gap, kernel="vectorized")
    slow = fill_matrices(a, b, matrix=DNA, gap=gap, kernel="reference")
    np.testing.assert_array_equal(fast.scores, slow.scores)
    np.testing.assert_array_equal(fast.pointers, slow.pointers)


@settings(max_examples=100, deadline=None)
@given(a=st.text(alphabet="ACGT", max_size=25), b=st.text(alphabet="ACGT", max_size=25))
def test_every_reported_alignment_rescores_to_optimum(a: str, b: str) -> None:
    result = needleman_wunsch(a, b, matrix=DNA, gap=GapPenalty(-3, -1))
    for alignment in result.alignments.take(50):
        assert alignment.aligned_a.replace("-", "") == result.dp.seq_a
        assert alignment.aligned_b.replace("-", "") == result.dp.seq_b
        assert rescore(alignment.aligned_a, alignment.aligned_b, result.dp) == result.score
        assert alignment.path[0] == (0, 0) and alignment.path[-1] == (len(a), len(b))


# ---------------------------------------------------------------------- Biopython cross-check

PROTEIN = "ARNDCQEGHILKMFPSTWYV"


def _biopython_aligner(matrix: Any, g_open: int, g_ext: int) -> Any:
    align = pytest.importorskip("Bio.Align")
    aligner = align.PairwiseAligner(mode="global")
    if isinstance(matrix, tuple):
        aligner.match_score, aligner.mismatch_score = matrix
    else:
        aligner.substitution_matrix = matrix
    aligner.open_gap_score, aligner.extend_gap_score = g_open, g_ext
    return aligner


@pytest.mark.parametrize("seed", range(40))
def test_dna_agrees_with_biopython(seed: int) -> None:
    rng = random.Random(seed)
    a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 120)))
    b = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 120)))
    match, mismatch = rng.randint(1, 5), rng.randint(-6, -1)
    g_open, g_ext = rng.randint(-12, -1), rng.randint(-4, 0)
    aligner = _biopython_aligner((match, mismatch), g_open, g_ext)
    theirs = aligner.align(a, b)
    ours = needleman_wunsch(
        a, b, matrix=SubstitutionMatrix.match_mismatch(match, mismatch), gap=GapPenalty(g_open, g_ext)
    )
    assert ours.score == int(theirs.score)
    assert ours.optimal_count == len(theirs)


@pytest.mark.parametrize("seed", range(20))
def test_protein_blosum62_agrees_with_biopython(seed: int) -> None:
    substitution_matrices = pytest.importorskip("Bio.Align.substitution_matrices")
    rng = random.Random(1000 + seed)
    a = "".join(rng.choice(PROTEIN) for _ in range(rng.randint(1, 150)))
    b = "".join(rng.choice(PROTEIN) for _ in range(rng.randint(1, 150)))
    aligner = _biopython_aligner(substitution_matrices.load("BLOSUM62"), -11, -1)
    theirs = aligner.align(a, b)
    ours = needleman_wunsch(a, b, matrix=named_matrix("BLOSUM62"), gap=GapPenalty(-11, -1))
    assert ours.score == int(theirs.score)
    assert ours.optimal_count == len(theirs)
