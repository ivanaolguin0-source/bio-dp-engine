"""Nussinov folding, checked against two references that share no code with the engine.

1. Exhaustive enumeration of every secondary structure, built from the first
   base of each interval (the engine decomposes on the last base).
2. The classic four-case recurrence from Durbin et al. (1998), for scores on
   inputs too long to enumerate.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

import pytest
from helpers import FIXTURES, load_json
from hypothesis import given, settings
from hypothesis import strategies as st
from oracles import classic_nussinov_score, is_valid_structure, max_pair_structures

from biodp.rna import FoldingRules, InvalidRnaError, base_pair_metrics, nussinov, parse_dot_bracket

GOLDEN = sorted((FIXTURES / "golden" / "nussinov").glob("*.json"))

# Yeast tRNA-Phe, unmodified 76-nt sequence and reference cloverleaf, as listed
# in the tRNA database at tpsic.igcz.poznan.pl (entry 520).
TRNA_PHE = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
TRNA_PHE_REFERENCE = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."


def rules_of(fixture: dict[str, Any]) -> FoldingRules:
    return FoldingRules(min_loop=fixture["input"]["min_loop"], wobble=fixture["input"]["wobble"])


# ---------------------------------------------------------------------- golden fixtures


def test_golden_fixtures_exist() -> None:
    assert len(GOLDEN) >= 8, "run scripts/generate_nussinov_golden.py"


@pytest.mark.parametrize("path", GOLDEN, ids=lambda p: Path(p).stem)
def test_golden_fixture(path: Path) -> None:
    fixture = load_json(path)
    result = nussinov(fixture["input"]["sequence"], rules=rules_of(fixture))
    expected = fixture["expected"]
    assert result.score == expected["score"]
    assert result.optimal_count == int(expected["optimal_count"])
    if expected["structures"] is not None:
        got = sorted([list(p) for p in s.pairs] for s in result.structures())
        assert got == expected["structures"]


# ---------------------------------------------------------------------- property tests

rna = st.text(alphabet="ACGU", max_size=12)


@settings(max_examples=300, deadline=None)
@given(seq=rna, min_loop=st.integers(min_value=0, max_value=4), wobble=st.booleans())
def test_matches_exhaustive_enumeration(seq: str, min_loop: int, wobble: bool) -> None:
    result = nussinov(seq, rules=FoldingRules(min_loop=min_loop, wobble=wobble))
    best, optimal = max_pair_structures(seq, min_loop, wobble)
    assert result.score == best
    assert result.optimal_count == len(optimal)
    produced = [s.pairs for s in result.structures()]
    assert len(produced) == len(set(produced)), "a structure was generated twice"
    assert set(produced) == optimal


@settings(max_examples=60, deadline=None)
@given(seq=st.text(alphabet="ACGU", max_size=60), min_loop=st.integers(0, 4), wobble=st.booleans())
def test_score_matches_classic_four_case_recurrence(seq: str, min_loop: int, wobble: bool) -> None:
    assert nussinov(seq, rules=FoldingRules(min_loop, wobble)).score == classic_nussinov_score(seq, min_loop, wobble)


@settings(max_examples=100, deadline=None)
@given(seq=st.text(alphabet="ACGU", max_size=40), min_loop=st.integers(0, 4), wobble=st.booleans())
def test_traceback_structures_are_valid_and_optimal(seq: str, min_loop: int, wobble: bool) -> None:
    result = nussinov(seq, rules=FoldingRules(min_loop, wobble))
    for structure in [result.best(), result.sample(random.Random(1)), *result.take(20)]:
        assert is_valid_structure(result.sequence, structure.pairs, min_loop, wobble)
        assert structure.pair_count == result.score
        assert parse_dot_bracket(structure.dot_bracket) == structure.pairs


# ---------------------------------------------------------------------- ties


def test_canonical_is_first_enumerated() -> None:
    result = nussinov("GGGGAAAACCCCAAAAGGGG")
    assert result.best().pairs == next(result.structures()).pairs


def test_sampling_is_uniform() -> None:
    result = nussinov("AAAAAAUUUUUU")  # 9 optimal structures
    assert result.optimal_count == 9
    rng = random.Random(20261025)
    tallies: dict[str, int] = {}
    draws = 9000
    for _ in range(draws):
        key = result.sample(rng).dot_bracket
        tallies[key] = tallies.get(key, 0) + 1
    assert len(tallies) == 9
    expected = draws / 9
    chi2 = sum((v - expected) ** 2 / expected for v in tallies.values())
    assert chi2 < 26.1  # roughly the 0.999 quantile of chi-square with 8 degrees of freedom


def test_traceback_visits_cells_inside_the_upper_triangle() -> None:
    result = nussinov(TRNA_PHE)
    visited = result.best().visited
    assert visited[0] == (0, len(TRNA_PHE) - 1)
    assert all(0 <= i <= j < len(TRNA_PHE) for i, j in visited)


# ---------------------------------------------------------------------- a real RNA


def test_trna_phe_against_reference_cloverleaf() -> None:
    """Nussinov over-pairs a real tRNA and recovers only part of its true structure.

    This is the honest headline for the method: maximising pair count is not
    how RNA folds. Score and the reference are checked independently; the
    sensitivity and PPV values are regression checks on this engine.
    """
    reference = parse_dot_bracket(TRNA_PHE_REFERENCE)
    assert len(TRNA_PHE) == len(TRNA_PHE_REFERENCE) == 76
    assert len(reference) == 21
    result = nussinov(TRNA_PHE)
    assert result.score == classic_nussinov_score(TRNA_PHE, 3, True)
    assert result.score > len(reference)
    metrics = base_pair_metrics(result.best(), reference)
    assert metrics["true_positives"] == 10
    assert metrics["sensitivity"] == pytest.approx(10 / 21)


# ---------------------------------------------------------------------- edge cases and inputs


@pytest.mark.parametrize("seq", ["", "A", "GC", "GAUC"])
def test_tiny_inputs(seq: str) -> None:
    result = nussinov(seq)
    assert (result.score, result.optimal_count) == (0, 1)
    assert result.best().dot_bracket == "." * len(seq)


def test_lowercase_is_accepted() -> None:
    assert nussinov("gggaaaucc").best().dot_bracket == "(((...)))"


def test_dna_t_is_rejected_with_a_hint() -> None:
    with pytest.raises(InvalidRnaError, match="replace it with 'U'"):
        nussinov("GGGAAATCC")


def test_other_symbols_rejected_with_position() -> None:
    with pytest.raises(InvalidRnaError, match="position 3"):
        nussinov("GGGNAAUCC")


def test_rules_validation() -> None:
    with pytest.raises(ValueError):
        FoldingRules(min_loop=-1)


def test_length_guard() -> None:
    with pytest.raises(ValueError, match="max_length"):
        nussinov("A" * 50, max_length=10)


@pytest.mark.parametrize("bad", ["(()", "())", "(.x)"])
def test_parse_dot_bracket_rejects_malformed(bad: str) -> None:
    with pytest.raises(ValueError):
        parse_dot_bracket(bad)


def test_pair_types() -> None:
    structure = nussinov("GGGAAAUUU").best()  # wobble on: G-U pairs
    assert structure.pair_types() == {"GC": 0, "AU": 0, "GU": 3}
