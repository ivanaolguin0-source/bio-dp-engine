<!--
  ALGORITHM DOCUMENT TEMPLATE

  Copy this file to docs/algorithms/<algorithm-slug>.md and replace every
  <placeholder>. Guidance lives in HTML comments like this one; they do not
  render on GitHub, so you can leave them in while drafting and delete them
  before publishing.

  Audience: a hiring manager who reads the first screen, and a technical
  interviewer who reads all of it. Write sections 0 to 2 for the first reader
  and sections 3 onward for the second.

  House rules:
  * Every claim about correctness points at a test that demonstrates it.
  * Every claim about speed points at a committed benchmark result.
  * Every citation carries enough detail to find it (authors, year, venue,
    volume, pages, DOI when verified). Mark anything you have not checked.
  * Math uses ```math fenced blocks for display equations and $...$ inline;
    both render on GitHub.
-->

# <Algorithm Name>

<!-- One sentence, no jargon a biologist would not know. What it computes and why anyone would want it. -->
<One-sentence summary.>

| | |
|---|---|
| **Problem** | <e.g. optimal global alignment of two sequences> |
| **Input** | <types and sizes> |
| **Output** | <score, structure, path, plus anything extra such as the count of co-optimal solutions> |
| **Time / space** | <e.g. $O(nm)$ / $O(nm)$; linear-space variant $O(n + m)$> |
| **Engine** | <Python or Java, and the module path> |
| **Status** | <implemented / tested / benchmarked> |

## 1. Biological significance

<!--
  Two or three short paragraphs. Answer, in order:
  1. What biological question does this answer? Name a concrete use (homology
     search, variant calling, structure prediction, gene annotation).
  2. Why is a *dynamic programming* solution required? Usually because the
     search space is exponential; quantify it.
  3. What does the score mean biologically? (Log-odds, base pairs, probability.)
  Avoid overselling: if production tools use a refined method, say so here.
-->

## 2. Problem statement

<!--
  Formal inputs and outputs, with the notation every later section uses.
  Define indices (0- or 1-based), alphabets, and the scoring model. State the
  objective as one line: "find X maximising Y subject to Z".
-->

## 3. DP formulation

### 3.1 Subproblems

<!--
  Define the table entry in words first, then symbols. "F(i, j) is the best
  score of ..." A reader should be able to predict the recurrence from this
  sentence alone. State why the optimal-substructure property holds.
-->

### 3.2 Recurrence relation

<!-- Display math. One equation per state. Explain each term in a sentence underneath. -->

```math
F(i, j) = \max \begin{cases}
  <\text{case 1}> \\
  <\text{case 2}>
\end{cases}
```

### 3.3 Base cases

<!-- The boundary row/column or diagonal, and the value that means "impossible" (usually -infinity). -->

### 3.4 Fill order and where the answer lives

<!--
  Which cells must be ready before (i, j)? Row-major, anti-diagonal, or by
  interval length (Nussinov)? Where is the final answer read from? If the
  engine uses a non-obvious order (wavefront), explain why it is still valid.
-->

## 4. Matrix backtracking

### 4.1 From table to answer

<!-- What is stored per cell to make reconstruction possible, and how the walk proceeds. -->

### 4.2 Ties and multiple optimal solutions

<!--
  State plainly whether the engine returns one solution or all of them, how
  ties are broken for the canonical answer, and why that convention was
  chosen. If you count co-optimal solutions, explain how.
-->

### 4.3 Worked example

<!--
  A small input, the filled table (or the relevant part), the path, and the
  final answer. Generate this from the engine and say so; never hand-type a
  table you have not computed.
-->

## 5. Time and space

### 5.1 Complexity

<!-- Asymptotic time and space, and the constant that matters in practice (bytes per cell). -->

### 5.2 Optimisations

<!--
  For each technique: what it saves, what it costs, and whether this repo
  implements it. Typical entries: linear-space score (two rows), divide and
  conquer (Hirschberg, Myers-Miller), banding, SIMD, parallel wavefronts.
-->

### 5.3 Measured performance

<!-- Table of sizes and timings, with machine, language version, and commit. Link the benchmark script. -->

## 6. Implementation map

<!-- Table: file path, responsibility. A reviewer should know where to click. -->

| File | Responsibility |
|---|---|
| `<path>` | <what it does> |

## 7. Validation

<!--
  How do we know it is right? List each independent source of evidence and
  the test that runs it: golden fixtures, exhaustive oracle, property tests,
  comparison with an established tool, cross-language agreement.
-->

## 8. Limitations

<!-- What this does not model, inputs it handles badly, and what a production tool would do instead. -->

## 9. References

<!--
  Author(s), year, title, venue, volume(issue), pages, DOI.
  Append "(DOI not verified)" to anything you have not resolved yourself.
-->
