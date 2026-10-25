# Hirschberg's Linear-Space Alignment

Finds an optimal global alignment of two sequences while holding only a couple of table rows in memory, instead of the whole table.

| | |
|---|---|
| **Problem** | Optimal global pairwise alignment with linear gap scores, in linear space |
| **Input** | Two sequences of lengths $n$ and $m$; a substitution matrix; one gap score $g \le 0$ |
| **Output** | One optimal alignment and its score; or, from `global_score`, the score alone |
| **Time / space** | $O(nm)$ / $O(n + m)$ |
| **Engine** | Python, [`engines/python/src/biodp/align/hirschberg.py`](../../engines/python/src/biodp/align/hirschberg.py) |
| **Status** | Implemented and tested; affine gaps (Myers-Miller) not implemented |

## 1. Biological significance

The biology is the same as for [Needleman-Wunsch](needleman-wunsch.md): a global alignment is a hypothesis about which positions of two related sequences descend from a common ancestor. What changes is the size of problem that can be answered.

Full-table Needleman-Wunsch stores a pointer for every cell so that it can trace back. In this repository that is 27 bytes per cell, so two 10,000-residue sequences (a long gene, or a small viral genome against a relative) need about 2.7 GB before any work is done. Hirschberg's method produces an optimal alignment of the same pair in under 2 MB.

The trade is real and should be stated: it takes somewhat longer, and it returns *one* optimal alignment. Counting or listing every co-optimal alignment needs the pointer graph this method exists to avoid storing.

## 2. Problem statement

Sequences $a = a_1 \dots a_n$ and $b = b_1 \dots b_m$, substitution score $s(x, y)$, and a linear gap score $g$: every gap column scores $g$, so a gap of length $L$ scores $L\,g$. Find an alignment of maximum total score. This is the linear-gap special case of the Needleman-Wunsch problem.

## 3. DP formulation

### 3.1 Subproblems

Hirschberg's method uses the ordinary Needleman-Wunsch recurrence, but only ever asks it for one row. Define

```math
F_{p}(k) = \text{best score of } a_1 \dots a_p \text{ against } b_1 \dots b_k ,
\qquad
B_{p}(k) = \text{best score of } a_{p+1} \dots a_n \text{ against } b_{k+1} \dots b_m .
```

$F$ is the forward table's row $p$. $B$ is the same thing computed on both sequences reversed, then read backwards.

### 3.2 Recurrence relation

Each row comes from the previous row only:

```math
H(i, j) = \max \begin{cases}
  H(i{-}1, j{-}1) + s(a_i, b_j) \\
  H(i{-}1, j) + g \\
  H(i, j{-}1) + g
\end{cases}
```

so a row of length $m + 1$ can be updated in place, keeping $O(m)$ memory.

The key step is the split. Take the middle row $p = \lfloor n/2 \rfloor$. Every alignment path crosses row $p$, at some column $k$, and the part above the crossing and the part below are independent. Therefore

```math
S^{*} = \max_{0 \le k \le m} \bigl( F_{p}(k) + B_{p}(k) \bigr),
\qquad
k^{*} = \arg\max_{0 \le k \le m} \bigl( F_{p}(k) + B_{p}(k) \bigr).
```

An optimal alignment of the whole is an optimal alignment of $(a_1 \dots a_p,\ b_1 \dots b_{k^*})$ followed by an optimal alignment of $(a_{p+1} \dots a_n,\ b_{k^*+1} \dots b_m)$. Recurse on both halves.

### 3.3 Base cases

| Subproblem | Answer |
|---|---|
| $a$ empty | all of $b$ against gaps |
| $b$ empty | all of $a$ against gaps |
| $a$ or $b$ has one residue | solved directly by the full-table engine (a $2 \times (m+1)$ table) |
| fewer than `base_cells` cells (default 4,096) | solved directly by the full-table engine |

The last rule is a practical one. Splitting all the way down to single residues is correct but spends most of its time in Python function calls on tiny problems. Because the threshold is a constant, it does not change the $O(n + m)$ memory bound. The tests set `base_cells=0` so that they exercise the textbook version.

### 3.4 Fill order and where the answer lives

Each row update is vectorised. The diagonal and vertical candidates come from the row above, and the horizontal dependency within the row unrolls to a running maximum:

```math
H(i, j) = j\,g + \max_{0 \le k \le j} \bigl( C(k) - k\,g \bigr),
\qquad
C(k) = \max\{ H(i{-}1, k{-}1) + s(a_i, b_k),\; H(i{-}1, k) + g \}.
```

`numpy.maximum.accumulate` computes that in one pass. The score is $F_n(m)$, the last entry of the last row; the alignment is the concatenation of the base-case alignments in left-to-right order.

## 4. Matrix backtracking

### 4.1 From table to answer

There is no traceback in the usual sense, because no pointers are kept. The recursion *is* the traceback: each split fixes one point that an optimal path passes through, and the base cases fill in the short segments between those points.

### 4.2 Ties and multiple optimal solutions

When $F_p(k) + B_p(k)$ attains its maximum at several columns, each is the crossing point of at least one optimal alignment. The engine takes the first (smallest $k$), so its result is always optimal but is not necessarily the same alignment as the Needleman-Wunsch engine's canonical one. The tests check exactly this: on random inputs, the Hirschberg result must be a member of the full set of co-optimal alignments that the Needleman-Wunsch engine enumerates.

### 4.3 Worked example

GATTACA against GCATGCU, match $+1$, mismatch $-1$, gap $-1$, computed by the engine. The middle row is $p = 3$, splitting GATTACA into GAT and TACA:

| $k$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| $F_3(k)$: GAT vs GCATGCU prefixes | −3 | −1 | −1 | 0 | 2 | 1 | 0 | −1 |
| $B_3(k)$: TACA vs GCATGCU suffixes | −3 | −2 | −1 | 0 | −2 | −2 | −4 | −4 |
| Sum | −6 | −3 | −2 | **0** | **0** | −1 | −4 | −5 |

$F_3$ is exactly row T (row 3) of the full table shown in the [Needleman-Wunsch write-up](needleman-wunsch.md#43-worked-example). The maximum, 0, is the optimal score, and it occurs at two columns. That tie is the three co-optimal alignments showing through: the canonical alignment crosses row 3 at column 3, and the second one crosses at column 4. The engine splits at $k^* = 3$ and, after recursing, returns

```
G-ATTACA
| | |.|.
GCA-TGCU
```

## 5. Time and space

### 5.1 Complexity

Space: the two rows used at any split hold $O(m)$ integers, the recursion is about $\log_2 n$ levels deep, and the output alignment is $O(n + m)$. Total $O(n + m)$.

Time: the top level scores the whole $n \times m$ area once (half forwards, half backwards). Its two subproblems have heights $n/2$ and widths summing to $m$, so together they cover half the area; their children a quarter; and so on.

```math
T(n, m) \le c\,nm \left(1 + \tfrac{1}{2} + \tfrac{1}{4} + \cdots\right) = 2\,c\,nm .
```

### 5.2 Optimisations

| Technique | Status |
|---|---|
| Vectorised row update with a running maximum | Implemented |
| Score-only mode (`global_score`), with the row laid along the shorter sequence | Implemented; swapping the sequences transposes the substitution matrix, so asymmetric matrices stay correct (tested) |
| Direct solve below a constant cell threshold | Implemented (`base_cells`) |
| Myers-Miller extension to affine gaps | Not implemented; affine input is rejected with an explanation |

### 5.3 Measured performance

Single run, random DNA with 10 % substitutions, linear gap $-2$, Python 3.13, NumPy 2.5.3, cloud container. Peak memory measured with `tracemalloc`.

| $n = m$ | Hirschberg time | Hirschberg peak memory | Full-table peak memory | Full-table time |
|---|---|---|---|---|
| 1,000 | 0.45 s | 0.21 MB | 28 MB | 0.36 s |
| 3,000 | 1.6 s | 0.49 MB | 245 MB | 1.2 s |
| 10,000 | 6.5 s | 1.6 MB | about 2.7 GB (refused by the size guard) | n/a |

At 3,000 residues: about 1/500 of the memory for about 1.4× the time. The score alone (`global_score`) took 0.83 s at 10,000.

## 6. Implementation map

| Function | Responsibility |
|---|---|
| `last_row` | One vectorised row sweep: all of $a$ against every prefix of $b$ |
| `global_score` | Optimal score in $O(\min(n, m))$ memory |
| `hirschberg` | The divide step, base cases, and assembly of the alignment |

## 7. Validation

All tests live in [`engines/python/tests/test_hirschberg.py`](../../engines/python/tests/test_hirschberg.py). The full-table engine serves as the reference, since it is itself verified against brute force and Biopython.

1. **Golden fixtures.** On every linear-gap fixture, the score equals the expected optimum and the alignment is one of the listed co-optimal alignments. Every affine-gap fixture is refused.
2. **Membership in the co-optimal set.** 300 random inputs with random scoring: the result must be one of the alignments the Needleman-Wunsch engine enumerates.
3. **Row correctness.** `last_row` must equal the bottom row of the full table exactly.
4. **Threshold independence.** The same score whether the base threshold is 0, 16, or 4,096 cells.
5. **Asymmetric matrices.** A deliberately asymmetric matrix checks that the shorter-row swap in `global_score` transposes correctly.
6. **Memory.** A 3,000 × 3,000 alignment must peak below 1/50 of the full table's size.
7. **The tests themselves.** Two deliberate bugs were planted in the divide step (ignoring the backward row; forgetting to re-reverse it). Both made the suite fail immediately.

## 8. Limitations

Linear gaps only, for a reason that is easy to miss. With affine gaps, an optimal alignment can cross the middle row in the middle of a gap. Splitting there and solving the halves independently would charge that one gap's opening score twice and could return a suboptimal alignment. Myers and Miller (1988) fix this by also tracking, at the split, whether the path is inside a gap.

One alignment per run. For the count of co-optimal alignments, or to list them, use the full-table engine.

The visualizer does not yet replay Hirschberg runs. The DPTrace schema already lists `hirschberg` as an algorithm; the natural picture is the table being carved into ever smaller rectangles around the path.

## 9. References

Hirschberg, D. S. (1975). A linear space algorithm for computing maximal common subsequences. *Communications of the ACM*, 18(6), 341–343. https://doi.org/10.1145/360825.360861

Myers, E. W., & Miller, W. (1988). Optimal alignments in linear space. *Computer Applications in the Biosciences*, 4(1), 11–17. https://doi.org/10.1093/bioinformatics/4.1.11

Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*, 48(3), 443–453. https://doi.org/10.1016/0022-2836(70)90057-4 (PMID: 5420325)

<!--
  Hirschberg's paper is posed for the longest common subsequence problem
  (see its abstract). Applying the same split to weighted alignment scores is
  the standard generalisation; cite the paper for the idea, not for the
  weighted version.
-->
