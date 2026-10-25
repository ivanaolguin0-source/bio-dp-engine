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
