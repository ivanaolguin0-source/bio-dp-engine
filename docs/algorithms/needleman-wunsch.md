# Needleman-Wunsch Global Alignment

Finds the best end-to-end alignment of two DNA, RNA, or protein sequences, and reports exactly how many equally good alignments exist.

| | |
|---|---|
| **Problem** | Optimal global pairwise alignment with affine gap scoring |
| **Input** | Two sequences of lengths $n$ and $m$; a substitution matrix; gap open and extend scores |
| **Output** | Optimal score; a canonical optimal alignment; the exact count of co-optimal alignments, with lazy enumeration and uniform sampling over them |
| **Time / space** | $O(nm)$ / $O(nm)$ for full traceback; $O(n + m)$ with [Hirschberg](hirschberg.md) (linear gaps) |
| **Engine** | Python, [`engines/python/src/biodp/align/`](../../engines/python/src/biodp/align/) |
| **Status** | Implemented; verified against brute-force enumeration and Biopython 1.88 |

## 1. Biological significance

Sequences that descend from a common ancestor accumulate substitutions, insertions, and deletions. An alignment is a hypothesis about that history: each column claims that two residues descend from the same ancestral position, or that one sequence gained or lost a residue there. Inferring homology, measuring divergence between orthologous genes, and checking a full-length amplicon against its expected sequence all start from such an alignment.

*Global* alignment assumes the two sequences are related along their whole length, which suits genes or proteins of similar length from related species. When only a region is shared (one conserved domain inside two otherwise different proteins), global alignment forces unrelated flanks into the answer, and local alignment (Smith-Waterman, planned for the Java engine) is the right tool.

Dynamic programming is not an optimisation here; it is the only practical approach. Two sequences of length 100 admit about $2.05 \times 10^{75}$ distinct alignments (the central Delannoy number $D(100, 100)$; the count grows roughly as $(3 + 2\sqrt{2})^n \approx 5.8^n$). The DP finds the best of them in $100 \times 100$ table cells.

With a log-odds matrix such as BLOSUM62, the score is a sum of log-likelihood ratios: each column compares how often that residue pair appears in alignments of related proteins against how often it would appear by chance. Gap scores model indels. They are *affine* (an expensive opening plus a cheaper extension) because a single mutational event often inserts or deletes several residues at once, so one gap of length 4 is more plausible than four separate gaps of length 1.

## 2. Problem statement

Let $a = a_1 \dots a_n$ and $b = b_1 \dots b_m$ be sequences over an alphabet $\Sigma$. Let $s : \Sigma \times \Sigma \to \mathbb{Z}$ be a substitution score, and let a gap of length $L \ge 1$ score

```math
\gamma(L) = o + (L - 1)\,e, \qquad o \le 0,\; e \le 0.
```

An alignment is a sequence of columns, each holding $(a_i, b_j)$, $(a_i, -)$, or $(-, b_j)$, that spells out all of $a$ on top and all of $b$ underneath, in order. Its score is the sum of $s$ over residue-residue columns plus $\gamma$ over every maximal run of same-kind gap columns. The task is to find an alignment of maximum score.

Two conventions are worth stating because implementations disagree on them:

1. **The gap convention.** Here $o$ is the score of the *first* gap position. This matches Biopython's `PairwiseAligner`. Some tools phrase gap penalties as positive costs, or charge $o + L\,e$; convert before comparing numbers.
2. **Adjacent opposite gaps are allowed.** A column $(a_i, -)$ may be followed directly by $(-, b_j)$, and the second gap pays its own opening score. Some implementations forbid this transition and can therefore miss optimal alignments (Flouri et al., 2015). Fixture `nw-004-adjacent-opposite-gaps` exercises it directly.

## 3. DP formulation

### 3.1 Subproblems

With linear gaps ($o = e = g$) one table suffices: $H(i, j)$ is the best score of any alignment of the prefixes $a_1 \dots a_i$ and $b_1 \dots b_j$. Optimal substructure holds because removing the last column of an optimal alignment leaves an optimal alignment of shorter prefixes; if it did not, swapping in a better one would improve the whole.

Affine gaps break that argument. Whether the next gap column costs $o$ or $e$ depends on whether the previous column was already a gap of the same kind, and $H$ alone forgets that. Gotoh's fix (1982) is to keep three tables, one per *state*, recording how the partial alignment ends:

| State | Meaning of the last column | Move in the table |
|---|---|---|
| $M(i, j)$ | $a_i$ opposite $b_j$ (match or substitution) | diagonal |
| $X(i, j)$ | $a_i$ opposite a gap | vertical (consumes $a$) |
| $Y(i, j)$ | $b_j$ opposite a gap | horizontal (consumes $b$) |

Each alignment corresponds to exactly one path through these states. That bijection is what makes exact counting of co-optimal alignments possible (section 4.2).

### 3.2 Recurrence relation

```math
\begin{aligned}
M(i,j) &= s(a_i, b_j) + \max\{\, M(i{-}1, j{-}1),\; X(i{-}1, j{-}1),\; Y(i{-}1, j{-}1) \,\} \\[4pt]
X(i,j) &= \max\{\, M(i{-}1, j) + o,\; X(i{-}1, j) + e,\; Y(i{-}1, j) + o \,\} \\[4pt]
Y(i,j) &= \max\{\, M(i, j{-}1) + o,\; X(i, j{-}1) + o,\; Y(i, j{-}1) + e \,\}
\end{aligned}
```

$M$ adds one aligned pair to the best alignment of the two shorter prefixes, whatever state it ended in. $X$ either extends a vertical gap already in progress (paying $e$) or opens a new one after any other kind of column (paying $o$). $Y$ is the mirror image. The $Y \to X$ and $X \to Y$ terms are the adjacent-opposite-gap transitions discussed above.

For linear gaps the three tables collapse into the familiar single recurrence, and the engine's exported projection $H = \max(M, X, Y)$ reproduces it exactly (tested in `test_projection_equals_textbook_table_for_linear_gaps`):

```math
H(i,j) = \max \begin{cases}
  H(i{-}1, j{-}1) + s(a_i, b_j) & \text{align } a_i \text{ with } b_j \\
  H(i{-}1, j) + g & a_i \text{ against a gap} \\
  H(i, j{-}1) + g & b_j \text{ against a gap}
\end{cases}
```

### 3.3 Base cases

```math
\begin{aligned}
M(0, 0) &= 0 && \text{the empty alignment} \\
X(i, 0) &= o + (i - 1)\,e && i \ge 1 \text{: one leading gap in } b \\
Y(0, j) &= o + (j - 1)\,e && j \ge 1 \text{: one leading gap in } a \\
\text{every other boundary cell} &= -\infty && \text{impossible states}
\end{aligned}
```

Getting these wrong is a classic source of bugs. Gotoh's original initialisation can produce suboptimal global alignments (Flouri et al., 2015). The engine represents $-\infty$ as a large negative `int64` sentinel with headroom against overflow, and serialises it as `null` in DPTrace output.

### 3.4 Fill order and where the answer lives

Cell $(i, j)$ needs $(i{-}1, j{-}1)$, $(i{-}1, j)$, and $(i, j{-}1)$. Row-major order satisfies this. The optimal score is

```math
S^{*} = \max\{\, M(n, m),\; X(n, m),\; Y(n, m) \,\}.
```

The default `"vectorized"` kernel fills a whole row per NumPy operation. $M$ and $X$ depend only on the previous row, so they vectorise directly. $Y$ depends on its left neighbour, which looks inherently sequential, but the recurrence unrolls to a running maximum:

```math
Y(i, j) = o + (j - 1)\,e + \max_{0 \le k < j} \bigl( B(k) - k\,e \bigr), \qquad B(k) = \max\{ M(i,k),\, X(i,k) \}
```

and `numpy.maximum.accumulate` computes that prefix maximum in one pass. Because all scores are integers the rewrite is exact, and a property test asserts that both kernels produce bit-identical score and pointer tensors.

### 3.5 How this differs from the 1970 paper

The recurrence above is the modern form, not a transcription of Needleman and Wunsch (1970). Their paper differs in three ways:

1. **Each cell searched a whole subrow and subcolumn.** A cell added its own similarity value to "the maximum value from among all the cells which lie on a pathway to it" (p. 444; the search along a subrow and subcolumn is pictured in Fig. 1, p. 445), so filling the array took more than constant work per cell. The three-neighbour form used here came later.
2. **The gap penalty was "a number subtracted for every gap made"** (p. 444), with the option of making it "a function of the size and/or direction of the gap." Affine gaps with a constant-time update are Gotoh's (1982) contribution.
3. **The array was filled from the ends of the sequences toward the origins** (p. 444; Fig. 1, p. 445, shows the fill proceeding "row-by-row towards the first row"), the reverse of the direction shown here. The optimum is the same either way.

The name stuck to the problem and its DP solution rather than to the exact 1970 procedure, which is why this document cites both papers.

## 4. Matrix backtracking

### 4.1 From table to answer

Alongside each score, the fill stores a 3-bit pointer mask per state recording *which* predecessor states attained the maximum (bit 1 from $M$, bit 2 from $X$, bit 4 from $Y$). The move direction is implied by the current state, so the mask only needs to name the predecessor state.

Traceback starts at whichever end states attain $S^{*}$, follows pointers back to $M(0, 0)$, and emits one column per step: $(a_i, b_j)$ on leaving an $M$ cell, $(a_i, -)$ on leaving $X$, $(-, b_j)$ on leaving $Y$. Reversing the emitted columns gives the alignment. The walk is iterative, so sequences of any length that fit in memory are safe from Python's recursion limit.

### 4.2 Ties and multiple optimal solutions

Most textbook implementations keep one pointer per cell and silently discard ties. That is a scientific problem, because co-optimal alignments are often biologically different hypotheses: in a homopolymer run, an indel has many equally good placements, and choosing one arbitrarily can make two pipelines report the same variant at different coordinates.

The pointer masks form a directed acyclic graph, and every path from an optimal end node to the start node is a distinct optimal alignment. The engine exploits this in three ways:

1. **Exact counting.** A second, small DP counts paths over just the *optimal subgraph* (nodes on at least one optimal path), bucketed by anti-diagonal $i + j$. Counts use Python's arbitrary-precision integers; `A`×200 against `A`×100 has $\binom{200}{100} \approx 9.05 \times 10^{58}$ optimal alignments, and the test suite checks that value exactly.
2. **Lazy enumeration.** Iterating `result.alignments` performs a depth-first walk with shared path tails, yielding alignments one at a time. Asking for five never touches the other $10^{58}$.
3. **Uniform sampling.** `result.alignments.sample(rng)` walks back choosing each predecessor with probability proportional to its path count. The step probabilities telescope to $1 / \text{count}$, so every optimal alignment is equally likely. A chi-square test checks this.

**Canonical tie-break.** Walking backwards, the engine prefers $M$, then $X$, then $Y$. Preferring the diagonal near the end defers gaps toward the start, so the canonical alignment places gaps as far left (toward the 5′ end) as the ties allow: `AAA` against `AA` yields `-AA` under `AAA`.

### 4.3 Worked example

GATTACA against GCATGCU, match $+1$, mismatch $-1$, linear gap $-1$. The projected table $H$ (rows are prefixes of GATTACA, columns of GCATGCU), generated by the engine:

```
        -    G    C    A    T    G    C    U
   -    0   -1   -2   -3   -4   -5   -6   -7
   G   -1    1    0   -1   -2   -3   -4   -5
   A   -2    0    0    1    0   -1   -2   -3
   T   -3   -1   -1    0    2    1    0   -1
   T   -4   -2   -2   -1    1    1    0   -1
   A   -5   -3   -3   -1    0    0    0   -1
   C   -6   -4   -2   -2   -1   -1    1    0
   A   -7   -5   -3   -1   -2   -2    0    0
```

The optimal score is $H(7, 7) = 0$, and there are exactly three optimal alignments:

```
G-ATTACA      G-ATTACA      G-ATTACA
GCA-TGCU      GCAT-GCU      GCATG-CU
```

A larger example: the peptides HEAGAWGHEE and PAWHEAE from Durbin et al. (1998), with BLOSUM50 and linear gap $-8$, score 1. The engine finds **three** co-optimal alignments. The canonical one is

```
HEAGAWGHE-E
  . || || |
--P-AW-HEAE
```

which matches the alignment I recall the textbook presenting (verify against your copy before quoting it). The other two differ only in where P and A sit among the leading gaps, which shows how much a single reported alignment can hide.

## 5. Time and space

### 5.1 Complexity

Filling is $\Theta(nm)$ time. Full traceback needs every pointer, so space is $\Theta(nm)$: three `int64` score layers and three `uint8` pointer layers, 27 bytes per cell. A 10,000 × 10,000 problem would need about 2.7 GB, which is why `needleman_wunsch` refuses tables above 25 million cells by default (about 675 MB) and suggests a linear-space method instead.

Counting co-optimal alignments costs time proportional to the optimal subgraph. For typical inputs that is a thin band around the main path; for tie-dominated inputs (long homopolymers) it can approach the full table.

### 5.2 Optimisations

| Technique | Saves | Costs | In this repo |
|---|---|---|---|
| Row-vectorised fill with prefix-max for $Y$ | Python interpreter overhead | Slightly less obvious code; guarded by kernel-equivalence test | Yes (default) |
| Score-only, two rows | Space to $O(\min(n, m))$ | No traceback possible | Yes: `global_score` (linear gaps) |
| Hirschberg divide and conquer | Space to $O(n + m)$ for the alignment itself | About 2× the time; one alignment, not all co-optimal ones | Yes: [`hirschberg`](hirschberg.md) (linear gaps) |
| Myers-Miller | Hirschberg's saving with affine gaps | Considerably fiddlier bookkeeping at the split point | Stretch goal |
| Banded alignment | Time and space to $O(k \cdot \max(n, m))$ | Exact only if the optimum stays within the band | Not planned |
| SIMD (striped, anti-diagonal) | Large constant factors | Platform-dependent code | Discussed for the Java engine |

Hirschberg's idea is worth stating because it is the canonical example of DP space reduction. Compute the last row of the forward scores for the top half of $a$ against all of $b$, and the last row of the *reverse* scores for the bottom half. The column $j^{*}$ where their sum peaks is where an optimal alignment crosses the middle row. Recurse on the two sub-rectangles. Every level of recursion touches half the area of the one above, so total time stays $O(nm)$ while only two rows are ever held.

### 5.3 Measured performance

Single run, no warm-up, random DNA with 10 % substitutions, affine gaps, Python 3.13, NumPy 2.5.3, in a cloud container. Treat these as orders of magnitude, not a benchmark.

| $n = m$ | Vectorized fill | Fill + optimal-subgraph count | Reference (loop) fill |
|---|---|---|---|
| 300 | 0.015 s | 0.016 s | 0.24 s |
| 1,000 | 0.08 s | 0.08 s | n/a |
| 3,000 | 0.41 s | 0.47 s | n/a |

## 6. Implementation map

| File | Responsibility |
|---|---|
| [`scoring/substitution.py`](../../engines/python/src/biodp/scoring/substitution.py) | Integer substitution matrices, NCBI-format loader, residue validation with exact error positions |
| [`scoring/gaps.py`](../../engines/python/src/biodp/scoring/gaps.py) | Affine gap model and its convention |
| [`align/model.py`](../../engines/python/src/biodp/align/model.py) | State enum, pointer flags, sentinel, the `DPMatrices` container |
| [`align/needleman_wunsch.py`](../../engines/python/src/biodp/align/needleman_wunsch.py) | Forward fill only: reference and vectorised kernels, public API |
| [`align/traceback.py`](../../engines/python/src/biodp/align/traceback.py) | Optimal-subgraph counting, lazy enumeration, uniform sampling, `Alignment` (CIGAR, identity, formatting) |
| [`trace/dptrace.py`](../../engines/python/src/biodp/trace/dptrace.py) | Serialisation to the DPTrace contract for the visualiser |
| [`contracts/schemas/dptrace.schema.json`](../../contracts/schemas/dptrace.schema.json) | The language-neutral contract every engine emits |

## 7. Validation

Correctness rests on four sources of evidence that share no code with the engine.

1. **Golden fixtures** (`fixtures/golden/needleman_wunsch/`). Eight hand-chosen inputs, each targeting one behaviour: empty inputs, a single affine gap, homopolymer ties, adjacent opposite gaps, the GATTACA example, gap merging under affine scoring, and two protein cases. Expected values come from exhaustive enumeration of every possible alignment, and the generator aborts if Biopython disagrees. The Java engine will consume the same files.
2. **Exhaustive oracle, property-based.** Hypothesis generates 300 random sequence pairs and scoring schemes; for each, the engine's score, count, and full *set* of optimal alignments must equal brute force. Enumeration must also never yield a duplicate.
3. **Biopython at scale.** Sixty random DNA and BLOSUM62 protein pairs up to 150 residues, with random affine gaps; score and co-optimal count must match `PairwiseAligner` exactly.
4. **Internal consistency.** Both fill kernels must agree bit for bit; every enumerated alignment must rescore to the optimum and spell out both input sequences; the linear-gap projection must equal an independently written single-matrix fill.

That last check earned its place during development: it caught the exporter drawing a spurious arrow at the origin cell.

## 8. Limitations

Scores are integers only. That is deliberate (exact ties are what make counting meaningful, and standard matrices are integer), but users with fractional scoring schemes must scale them.

End gaps are penalised like interior gaps. Semi-global ("glocal") alignment, where overhangs are free, is common when aligning a read to a longer reference and is a natural extension: it changes only the base cases and the answer cell.

A single substitution matrix scores every position the same way, so it cannot express per-position profiles; aligning to a protein family calls for profile HMMs instead.

Full-matrix traceback limits practical input size to tens of thousands of residues per side. Genome-scale alignment uses seed-and-extend heuristics, not exhaustive DP.

## 9. References

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163 (PMID: 19304878; PMCID: PMC2682512)

Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). *Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids*. Cambridge University Press.

Flouri, T., Kobert, K., Rognes, T., & Stamatakis, A. (2015). Are all global alignment algorithms and implementations correct? *bioRxiv*. https://doi.org/10.1101/031500

Gotoh, O. (1982). An improved algorithm for matching biological sequences. *Journal of Molecular Biology*, 162(3), 705–708. https://doi.org/10.1016/0022-2836(82)90398-9 (PMID: 7166760)

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89(22), 10915–10919. https://doi.org/10.1073/pnas.89.22.10915

Hirschberg, D. S. (1975). A linear space algorithm for computing maximal common subsequences. *Communications of the ACM*, 18(6), 341–343. https://doi.org/10.1145/360825.360861

Myers, E. W., & Miller, W. (1988). Optimal alignments in linear space. *Computer Applications in the Biosciences*, 4(1), 11–17. https://doi.org/10.1093/bioinformatics/4.1.11

Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*, 48(3), 443–453. https://doi.org/10.1016/0022-2836(70)90057-4 (PMID: 5420325)
