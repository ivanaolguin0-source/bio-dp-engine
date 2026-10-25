from __future__ import annotations

import json

import jsonschema
import numpy as np
import pytest
from helpers import SCHEMA_PATH, load_json

from biodp.align import needleman_wunsch
from biodp.scoring import GapPenalty, SubstitutionMatrix
from biodp.trace import dumps, to_dptrace

DNA = SubstitutionMatrix.match_mismatch(1, -1, "ACGTU")
VALIDATOR = jsonschema.Draft202012Validator(load_json(SCHEMA_PATH))


@pytest.mark.parametrize(
    ("a", "b", "gap"),
    [("GATTACA", "GCATGCU", GapPenalty.linear(-1)), ("", "", GapPenalty.linear(-1)), ("ACGT", "", GapPenalty(-3, -1))],
)
def test_trace_validates_against_contract(a: str, b: str, gap: GapPenalty) -> None:
    document = to_dptrace(needleman_wunsch(a, b, matrix=DNA, gap=gap))
    VALIDATOR.validate(json.loads(dumps(document)))


def textbook_single_matrix(a: str, b: str, match: int, mismatch: int, g: int) -> tuple[list[list[int]], list[list[int]]]:
    """The classic one-matrix linear-gap fill, written independently for comparison."""
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        H[i][0], D[i][0] = i * g, 2
    for j in range(1, m + 1):
        H[0][j], D[0][j] = j * g, 4
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up, left = H[i - 1][j] + g, H[i][j - 1] + g
            H[i][j] = best = max(diag, up, left)
            D[i][j] = (diag == best) * 1 | (up == best) * 2 | (left == best) * 4
    return H, D


def test_projection_equals_textbook_table_for_linear_gaps() -> None:
    a, b = "GATTACA", "GCATGCU"
    document = to_dptrace(needleman_wunsch(a, b, matrix=DNA, gap=GapPenalty.linear(-1)))
    H, D = textbook_single_matrix(a, b, 1, -1, -1)
    assert document["projection"]["scores"] == H
    assert document["projection"]["moves"] == D


def test_optimal_count_is_a_string_and_list_is_truncated() -> None:
    result = needleman_wunsch("A" * 40, "A" * 20, matrix=DNA, gap=GapPenalty.linear(-1))
    document = to_dptrace(result, max_alignments=4)
    assert document["result"]["optimal_count"] == str(result.optimal_count)
    assert len(document["result"]["alignments"]) == 4
    assert document["result"]["alignments_truncated"] is True


def test_on_path_mask_covers_every_listed_path() -> None:
    result = needleman_wunsch("GATTACA", "GCATGCU", matrix=DNA, gap=GapPenalty.linear(-1))
    document = to_dptrace(result)
    mask = np.array(document["projection"]["on_optimal_path"])
    union = np.zeros_like(mask)
    for alignment in document["result"]["alignments"]:
        for i, j in alignment["path"]:
            union[i, j] = True
    assert not document["result"]["alignments_truncated"]
    np.testing.assert_array_equal(mask, union)


def test_unreachable_cells_serialise_as_null() -> None:
    document = to_dptrace(needleman_wunsch("AC", "AG", matrix=DNA, gap=GapPenalty(-3, -1)))
    m_layer = next(layer for layer in document["layers"] if layer["name"] == "M")
    assert m_layer["scores"][0][1] is None and m_layer["scores"][1][0] is None


def test_visualisation_size_cap() -> None:
    result = needleman_wunsch("A" * 400, "A" * 400, matrix=DNA, gap=GapPenalty(-5, -1))
    with pytest.raises(ValueError, match="capped"):
        to_dptrace(result)


# ---------------------------------------------------------------------- Nussinov traces

from biodp.rna import nussinov  # noqa: E402

TRNA = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"


@pytest.mark.parametrize("seq", ["GGGAAAUCC", "A", "GGGGAAAACCCCAAAAGGGG", TRNA])
def test_nussinov_trace_validates(seq: str) -> None:
    document = to_dptrace(nussinov(seq))
    VALIDATOR.validate(json.loads(dumps(document)))
    assert document["fill_order"]["kind"] == "interval-length"
    assert document["shape"] == [len(seq), len(seq)]


def test_nussinov_trace_structure_of_table() -> None:
    seq = "GGGGAAAACCCCAAAAGGGG"
    result = nussinov(seq)
    document = to_dptrace(result)
    scores, moves = document["projection"]["scores"], document["projection"]["moves"]
    n = len(seq)
    for i in range(n):
        for j in range(n):
            if j < i:
                assert scores[i][j] is None and moves[i][j] == 0
            else:
                assert scores[i][j] == result.N(i, j)
                assert (moves[i][j] == 0) == (result.N(i, j) == 0)
    assert scores[0][n - 1] == document["result"]["score"] == 4
    assert document["result"]["optimal_count"] == "70"
    assert document["result"]["structures_truncated"] is True


def test_nussinov_listed_structures_replay_consistently() -> None:
    result = nussinov("AAAAAAUUUUUU")
    document = to_dptrace(result)
    on_path = document["projection"]["on_optimal_path"]
    assert len(document["result"]["structures"]) == 9
    for entry in document["result"]["structures"]:
        pairs = [tuple(p) for p in entry["pairs"]]
        assert entry["dot_bracket"].count("(") == len(pairs) == result.score
        assert entry["visited"][0] == [0, 11]
        assert all(on_path[i][j] for i, j in entry["visited"])
        # Every pair (k, j) must be decided at a visited cell whose right end is j.
        visited_right_ends = {j for _, j in entry["visited"]}
        assert all(j in visited_right_ends for _, j in pairs)


def test_nussinov_empty_sequence_is_refused() -> None:
    with pytest.raises(ValueError, match="empty"):
        to_dptrace(nussinov(""))
