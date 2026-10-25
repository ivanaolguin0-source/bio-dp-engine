# Nussinov RNA Secondary Structure Prediction

Folds a single RNA strand onto itself by choosing the set of base pairs that is largest, among all structures without crossing pairs, and reports exactly how many structures tie for that maximum.

| | |
|---|---|
| **Problem** | Maximum base-pairing secondary structure of one RNA sequence |
| **Input** | An RNA sequence of length $n$; a minimum hairpin loop length $h$ (default 3); whether G-U wobble pairs are allowed (default yes) |
| **Output** | Maximum number of pairs; a canonical optimal structure in dot-bracket notation; the exact count of optimal structures, with lazy enumeration and uniform sampling |
| **Time / space** | $O(n^3)$ / $O(n^2)$ |
| **Engine** | Python, [`engines/python/src/biodp/rna/nussinov.py`](../../engines/python/src/biodp/rna/nussinov.py) |
| **Status** | Implemented, tested, visualized |

## 1. Biological significance

An RNA molecule is a single strand, and it folds back on itself: complementary stretches pair up (G with C, A with U, and the weaker G-U "wobble"), forming double-stranded stems closed by single-stranded loops. That pattern of pairs is the *secondary structure*, and for many RNAs it is what matters. Transfer RNAs fold into a cloverleaf of four stems; ribosomal RNAs, riboswitches and microRNA precursors are all recognised by the shape their pairs make.

Predicting which bases pair is a search over an enormous space: the number of possible structures grows exponentially with length. The problem becomes tractable because of one restriction. If no two pairs cross (no $i < k < j < l$ with pairs $(i, j)$ and $(k, l)$), then a pair splits the sequence into an inside and an outside that fold independently. Crossing pairs, called pseudoknots, do occur in real RNAs, but excluding them is what makes dynamic programming possible.

**The honest caveat.** Nussinov maximises the *number* of pairs. Real RNA folds toward minimum free energy, where what matters is stacking between adjacent pairs and the cost of loops. Practical predictors use energy models (Zuker and Stiegler, 1981; implemented in the ViennaRNA package). Nussinov is the cleanest example of *interval* dynamic programming, and the stepping stone to those models; it is not a realistic predictor, and section 4.3 shows exactly how unrealistic.

## 2. Problem statement

Let $x = x_1 \dots x_n$ be an RNA sequence over $\{A, C, G, U\}$. A pair $(i, j)$ with $i < j$ is *allowed* if $x_i x_j$ is A-U, U-A, G-C, C-G (or G-U, U-G with wobble enabled) and it encloses a hairpin loop of at least $h$ bases, $j - i - 1 \ge h$. A secondary structure is a set of allowed pairs in which each base pairs at most once and no two pairs cross. Find a structure with the largest number of pairs.

DNA's T is rejected rather than silently converted to U: a silent conversion would hide a wrong input.

## 3. DP formulation

### 3.1 Subproblems

$N(i, j)$ is the maximum number of pairs in a structure on the interval $x_i \dots x_j$ alone. This is *interval DP*: the subproblems are substrings, not prefixes as in alignment. Optimal substructure holds because, with no crossings, a pair $(k, j)$ divides the interval into $x_i \dots x_{k-1}$ and $x_{k+1} \dots x_{j-1}$, which cannot interact.

### 3.2 Recurrence relation

The textbook form (Durbin et al., 1998) considers four ways the interval can be built:

```math
N(i, j) = \max \begin{cases}
  N(i{+}1, j) & i \text{ unpaired} \\
  N(i, j{-}1) & j \text{ unpaired} \\
  N(i{+}1, j{-}1) + 1 & i \text{ pairs with } j \\
  \max_{i < k < j} \bigl( N(i, k) + N(k{+}1, j) \bigr) & \text{bifurcation}
\end{cases}
```

This finds the right *score*, but most structures can be built through several of these cases (an unpaired base at both ends can be removed from either end first, and a bifurcation can usually be placed at several $k$). Counting derivations would then overcount structures.

The engine uses the standard *unambiguous* equivalent, which decides only the fate of the last base $j$:

```math
N(i, j) = \max \Bigl( N(i, j{-}1),\;\; \max_{\substack{i \le k \le j - h - 1 \\ (k, j) \text{ allowed}}} \bigl( N(i, k{-}1) + 1 + N(k{+}1, j{-}1) \bigr) \Bigr)
```

Either $j$ is unpaired, or it pairs with exactly one $k$, and that choice splits off the inside $(k+1, j-1)$ and the left part $(i, k-1)$. Every structure has exactly one derivation, so counting optimal derivations counts optimal structures exactly. The test suite checks that both recurrences give the same maximum on hundreds of random sequences.

### 3.3 Base cases

```math
N(i, j) = 0 \quad \text{whenever } j - i \le h
```

An interval of $h + 1$ bases or fewer cannot close a hairpin, so it has no pairs. Empty intervals ($j < i$) also score 0 and contain exactly one structure, the empty one.

### 3.4 Fill order and where the answer lives

$N(i, j)$ reads only shorter intervals, so the table is filled one diagonal at a time, by interval length $d = j - i$ from 0 up to $n - 1$. Only the upper triangle is used. In the visualizer the cell being filled reads an L-shaped set of cells: along row $i$ to its left (the $N(i, k-1)$ terms) and down column $j - 1$ (the $N(k+1, j-1)$ terms).

The answer is $N(1, n)$, the top-right corner.

## 4. Matrix backtracking

### 4.1 From table to answer

Traceback starts at the whole sequence $(1, n)$ and asks which choice achieved the maximum: $j$ unpaired, so continue with $(i, j-1)$; or $j$ paired with $k$, so record the pair and continue with both $(i, k-1)$ and $(k+1, j-1)$. Because a pair spawns two sub-intervals, the traceback is a *tree*, not a path; the engine walks it with an explicit stack. Each interval it visits has a different right end, so the visualizer can reveal a pair the moment its interval is reached.

### 4.2 Ties and multiple optimal solutions

Ties are common, because many different pairings reach the same count. The engine records the count of optimal structures for every interval alongside its score, with arbitrary-precision integers: a random 1,000-base sequence has a 30-digit count. It then offers:

- the **canonical** structure: at each step prefer "$j$ unpaired", then the smallest partner $k$;
- **lazy enumeration** of every optimal structure, each exactly once;
- **uniform sampling**, weighting each choice by the number of structures it leads to.

### 4.3 Worked example: a real tRNA

Yeast tRNA-Phe, 76 nucleotides, with its reference cloverleaf from the tRNA database at tpsic.igcz.poznan.pl (entry 520; unmodified sequence). Results from the engine:

| | Nussinov (canonical) | Reference cloverleaf |
|---|---|---|
| Base pairs | 28 | 21 |
| Optimal structures | 1,872 tie | (one known structure) |
| Pairs shared | 10 | 10 |
| Sensitivity (shared / reference) | 48% | |
| PPV (shared / predicted) | 36% | |

```
sequence   GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA
predicted  (((((((((((...)))(((((((((((((((....).)))))((...)).)))))..))).)))).)))))....
reference  (((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....
```

Maximising the pair count produces more pairs than the molecule really has, recovers fewer than half of the true ones, and cannot even choose among 1,872 equally good answers. The real cloverleaf has fewer pairs arranged in stacked stems, which an energy model rewards and a pair count does not. The visualizer's "Yeast tRNA-Phe" case draws the two structures against each other.

## 5. Time and space

### 5.1 Complexity

There are $O(n^2)$ intervals, and each considers $O(n)$ partners $k$: $O(n^3)$ time, $O(n^2)$ space. The engine vectorises the scan over $k$ with NumPy; the counts use Python integers, since they routinely exceed 64 bits.

### 5.2 Optimisations

| Technique | Status |
|---|---|
| Vectorised scan over partners $k$ | Implemented |
| Skip partners that cannot pair or would close a too-small loop | Implemented (only allowed pairs are scored) |
| Energy model (Zuker) | Out of scope; this is the natural next algorithm in this family |
| Restricting pair span (local folding) | Not implemented |

### 5.3 Measured performance

Random sequences, single run, Python 3.13, NumPy 2.5.3, cloud container. Times include the structure counts.

| $n$ | Time | Maximum pairs | Digits in the optimal-structure count |
|---|---|---|---|
| 76 | 0.03 s | 27 | 5 |
| 200 | 0.22 s | 75 | 11 |
| 500 | 1.4 s | 191 | 18 |
| 1,000 | 7.7 s | 396 | 30 |

## 6. Implementation map

| Item | Responsibility |
|---|---|
| `FoldingRules` | Minimum hairpin loop and wobble setting |
| `nussinov()` | Fill by interval length; scores and structure counts |
| `NussinovResult.best / structures / sample` | Canonical traceback, enumeration, uniform sampling |
| `NussinovResult.choices`, `optimal_cells`, `visited_for` | The decisions at each cell, for the trace export and visualizer |
| `parse_dot_bracket`, `base_pair_metrics` | Reading reference structures; sensitivity and PPV |
| `biodp.trace.to_dptrace` | DPTrace export (schema 1.1: interval-length fill, structures, move labels) |
| [`web/src/viz/ArcDiagram.tsx`](../../web/src/viz/ArcDiagram.tsx) | Arc diagram: predicted pairs above the line, reference below |

## 7. Validation

All in [`engines/python/tests/test_nussinov.py`](../../engines/python/tests/test_nussinov.py) unless noted.

1. **Golden fixtures** (`fixtures/golden/nussinov/`): nine hand-chosen cases (empty, too short to pair, a simple hairpin, no loop constraint, wobble on and off, two hairpins side by side, a homopolymer stem, and one with 70 tied structures). Expected values come from exhaustive enumeration and are cross-checked against the classic recurrence.
2. **An independent oracle.** The brute-force enumerator builds structures by deciding the *first* base's partner; the engine decides the *last* base's. On 300 random sequences with random rules, the score, the count, and the full set of optimal structures must match exactly, with no duplicates.
3. **The classic recurrence.** On sequences up to 60 bases, the four-case textbook recurrence must give the same score.
4. **Validity.** Every structure the engine returns is checked for allowed pairs, minimum loops, single pairing, and no crossings.
5. **Uniform sampling.** A chi-square test on 9,000 draws over 9 tied structures.
6. **The tests themselves.** Four deliberate bugs were planted (hairpin loop off by one, the "j unpaired" count dropped, the inner interval forgotten in the score, and the inner interval skipped in traceback). Each made the suite fail.
7. **Trace contract** (`test_dptrace.py`, and `web/src/trace/validate.test.ts`): folding traces validate against the shared schema; the table's lower triangle is null; every listed structure's traceback lies on optimal cells.

## 8. Limitations

No pseudoknots, by construction.

A pair count is not an energy, so predictions for real RNAs are poor (section 4.3).

Every allowed pair counts the same: G-C, with three hydrogen bonds, counts no more than G-U. Weighting pairs is a small change to the recurrence; an energy model is the real fix.

Memory grows as $n^2$ and time as $n^3$; the default limit is 2,000 bases.

## 9. References

Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). *Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids*. Cambridge University Press.

Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26. https://doi.org/10.1186/1748-7188-6-26 (PMID: 22115189; PMCID: PMC3319429)

Nussinov, R., & Jacobson, A. B. (1980). Fast algorithm for predicting the secondary structure of single-stranded RNA. *Proceedings of the National Academy of Sciences*, 77(11), 6309–6313. https://doi.org/10.1073/pnas.77.11.6309 (PMID: 6161375; PMCID: PMC350273)

Nussinov, R., Pieczenik, G., Griggs, J. R., & Kleitman, D. J. (1978). Algorithms for loop matchings. *SIAM Journal on Applied Mathematics*, 35(1), 68–82. https://doi.org/10.1137/0135006

Zuker, M., & Stiegler, P. (1981). Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. *Nucleic Acids Research*, 9(1), 133–148. https://doi.org/10.1093/nar/9.1.133 (PMID: 6163133; PMCID: PMC326673)

Reference tRNA structure: tRNA database at tpsic.igcz.poznan.pl, entry 520, *Saccharomyces cerevisiae* tRNA-Phe.
